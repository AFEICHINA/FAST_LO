#ifndef _STATE_PROCESSING_HPP
#define _STATE_PROCESSING_HPP

#include <cmath>
#include <math.h>
#include <deque>
#include <mutex>
#include <thread>
#include <fstream>
#include <csignal>
#include <ros/ros.h>
#include <so3_math.h>
#include <Eigen/Eigen>
#include <common_lib.h>
#include <pcl/common/io.h>
#include <pcl/point_cloud.h>
#include <pcl/point_types.h>
#include <condition_variable>
#include <nav_msgs/Odometry.h>
#include <pcl/common/transforms.h>
#include <pcl/kdtree/kdtree_flann.h>
// #include <tf/transform_broadcaster.h>
#include <eigen_conversions/eigen_msg.h>
#include <pcl_conversions/pcl_conversions.h>
// #include <sensor_msgs/Imu.h>
#include <sensor_msgs/PointCloud2.h>
// #include <lidar_imu_init/States.h>
// #include <geometry_msgs/Vector3.h>

bool b_first_frame_ = true;
double time_last_scan;
V3D cov_acc_scale = V3D(0.1, 0.1, 0.1);
V3D cov_gyr_scale = V3D(0.1, 0.1, 0.1);

void set_gyr_cov(const V3D &scaler)
{
  cov_gyr_scale = scaler;
}

void set_acc_cov(const V3D &scaler)
{
  cov_acc_scale = scaler;
}

// void set_mean_acc_norm(const double &mean_acc_norm){
//     IMU_mean_acc_norm = mean_acc_norm;
// }

// void set_gyr_bias_cov(const V3D &b_g)
// {
//   cov_bias_gyr = b_g;
// }

// void set_acc_bias_cov(const V3D &b_a)
// {
//   cov_bias_acc = b_a;
// }

void Forward_propagation_without_imu(const MeasureGroup &meas, StatesGroup &state_inout,
                             PointCloudXYZI &pcl_out) {
    pcl_out = *(meas.lidar);
    /*** sort point clouds by offset time ***/
    const double &pcl_beg_time = meas.lidar_beg_time;
    sort(pcl_out.points.begin(), pcl_out.points.end(), time_list);
    const double &pcl_end_offset_time = pcl_out.points.back().curvature / double(1000);

    MD(DIM_STATE, DIM_STATE) F_x, cov_w;
    double dt = 0.0;

    if (b_first_frame_) {
        dt = 0.1;
        b_first_frame_ = false;
    } else {
        dt = pcl_beg_time - time_last_scan;
        time_last_scan = pcl_beg_time;
    }

    /* covariance propagation */
    F_x.setIdentity();
    cov_w.setZero();
    /** In CV model, bias_g represents angular velocity **/
    /** In CV model，bias_a represents linear acceleration **/
    M3D Exp_f = Exp(state_inout.bias_g, dt);
    F_x.block<3, 3>(0, 0) = Exp(state_inout.bias_g, -dt);
    F_x.block<3, 3>(0, 15) = Eye3d * dt;
    F_x.block<3, 3>(3, 12) = Eye3d * dt;


    cov_w.block<3, 3>(15, 15).diagonal() = cov_gyr_scale * dt * dt;
    cov_w.block<3, 3>(12, 12).diagonal() = cov_acc_scale * dt * dt;

    /** Forward propagation of covariance**/
    state_inout.cov = F_x * state_inout.cov * F_x.transpose() + cov_w;

    /** Forward propagation of attitude **/
    state_inout.rot_end = state_inout.rot_end * Exp_f;

    /** Position Propagation **/
    state_inout.pos_end += state_inout.vel_end * dt;

    /**CV model： un-distort pcl using linear interpolation **/
    // if(lidar_type != L515){
        auto it_pcl = pcl_out.points.end() - 1;
        double dt_j = 0.0;
        for(; it_pcl != pcl_out.points.begin(); it_pcl --)
        {
            dt_j= pcl_end_offset_time - it_pcl->curvature/double(1000);
            M3D R_jk(Exp(state_inout.bias_g, - dt_j));
            V3D P_j(it_pcl->x, it_pcl->y, it_pcl->z);
            // Using rotation and translation to un-distort points
            V3D p_jk;
            p_jk = - state_inout.rot_end.transpose() * state_inout.vel_end * dt_j;

            V3D P_compensate =  R_jk * P_j + p_jk;

            /// save Undistorted points and their rotation
            it_pcl->x = P_compensate(0);
            it_pcl->y = P_compensate(1);
            it_pcl->z = P_compensate(2);
        }
    // }
}

#endif