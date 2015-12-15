/*
Copyright 2011. All rights reserved.
Institute of Measurement and Control Systems
Karlsruhe Institute of Technology, Germany

This file is part of libviso2.
Authors: Andreas Geiger

libviso2 is free software; you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation; either version 2 of the License, or any later version.

libviso2 is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
libviso2; if not, write to the Free Software Foundation, Inc., 51 Franklin
Street, Fifth Floor, Boston, MA 02110-1301, USA 
*/

#ifndef VISUALODOMETRY_H
#define VISUALODOMETRY_H

#include <iostream>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "matcher.h"
#include "matrix.h"

class VisualOdometry {

public:

  // constructor
  VisualOdometry();
  
  // deconstructor
  ~VisualOdometry();

  // set calibration parameters. you have to call this function
  // at least once after instanciating the class and before starting
  // the visual odometry process. to estimate scale, the baseline
  // is used in the stereo case. for the monocular case, we assume
  // that the camera has no roll, and estimate the scale from the camera
  // height over ground by fitting the ground plane to SfM features.
  // input: f ........ focal length (assumes fu=fv)
  //        cu, cv ... principal point
  //        b ........ baseline (stereo) / height over ground (mono)
  //        p ........ camera pitch (only needed for mono version)
  void setCalibration(FLOAT f,FLOAT cu,FLOAT cv,FLOAT b,FLOAT p=0);
  
  // takes matched features as input and updates the visual odometry estimate.
  // input:  p_matched ............. matched features
  //         deltaT ................ time difference between the 2 frames
  //         stereo ................ use stereo for estimation
  //         record_raw_odometry ... debugging flag, not needed
  // output: true on success
  bool update (std::vector<Matcher::p_match> p_matched,
               FLOAT deltaT,bool stereo,bool record_raw_odometry=false);
  
  // returns the 4x4 transformation matrix computed by the last update
  // call. the transformation matrix projects points from the previous
  // camera coordinates to the current camera coordinates.
  Matrix getTransformation ();
  
  // get inlier feature match indices (those which are consistent with egomotion)
  std::vector<int32_t> getInliers () { return inliers; }

private:
  
  template<class T> struct idx_cmp {
    idx_cmp(const T arr) : arr(arr) {}
    bool operator()(const size_t a, const size_t b) const { return arr[a] < arr[b]; }
    const T arr;
  };

  enum result { UPDATED, FAILED, CONVERGED };

  bool    processMono(std::vector<Matcher::p_match> p_matched,FLOAT deltaT);
  bool    processStereo(std::vector<Matcher::p_match> p_matched,FLOAT deltaT);
  void    parcpy(FLOAT *param_dst,FLOAT *param_src);
  
  // monocular version
  Matrix  smallerThanMedian(Matrix &X,FLOAT &median);
  bool    monoNormalizeFeaturePoints(std::vector<Matcher::p_match> &p_matched,Matrix &T1,Matrix &T2);
  std::vector<int32_t> monoGetRandomSample(int32_t N,int32_t num);
  void    monoFundamentalMatrix (const std::vector<Matcher::p_match> &p_matched,const std::vector<int32_t> &active,Matrix &F);
  std::vector<int32_t> monoGetFundamentalInlier (std::vector<Matcher::p_match> &p_matched,Matrix &F,FLOAT tau);
  void    monoEtoRt(Matrix &E,Matrix &K,std::vector<Matcher::p_match> &p_matched,std::vector<int32_t> &inliers,Matrix &X,Matrix &R,Matrix &t);
  int32_t monoTriangulateChieral (std::vector<Matcher::p_match> &p_matched,Matrix &K,Matrix &R,Matrix &t,Matrix &X);
  
  // stereo version
  void    stereoPlotErrors(const std::vector<Matcher::p_match> p_matched,const std::vector<int32_t> &active,FLOAT *param);
  std::vector<int32_t> stereoGetInlier(const std::vector<Matcher::p_match> p_matched,const FLOAT *param,FLOAT tau);
  std::vector<int32_t> stereoGetRandomSample(const std::vector<Matcher::p_match> &p_matched,const float &min_dist);
  result  monoUpdateParameters(const std::vector<Matcher::p_match> &p_matched,const std::vector<int32_t> &active,FLOAT *param,const FLOAT &step_size,const FLOAT &eps);
  result  stereoUpdateParameters(const std::vector<Matcher::p_match> &p_matched,const std::vector<int32_t> &active,FLOAT *param,const FLOAT &step_size,const FLOAT &eps);
  void    monoComputeObservations(const std::vector<Matcher::p_match> &p_matched,const std::vector<int32_t> &active);
  void    stereoComputeObservations(const std::vector<Matcher::p_match> &p_matched,const std::vector<int32_t> &active);
  void    monoComputePredictionsAndJacobian(const FLOAT *param,int32_t num_pt);
  void    stereoComputePredictionsAndJacobian(const FLOAT *param,const std::vector<int32_t> &active);
  void    stereoTestJacobian();

  // focal length, principal point, baseline/camera height, pitch
  FLOAT calib_f,calib_cu,calib_cv,calib_b,calib_p;

  FLOAT *X,*Y,*Z;               // 3d points
  FLOAT *J;                     // jacobian
  FLOAT *p_observe,*p_predict;  // observed and predicted 2d points
  FLOAT param[6];               // parameter set
  std::vector<int32_t> inliers; // inlier set

  // kalman filter
  Matrix KF_x; // state vector
  Matrix KF_z; // observation vector
  Matrix KF_A; // state transition matrix
  Matrix KF_H; // observation matrix
  Matrix KF_P; // covariance of state vector
  Matrix KF_Q; // process noise
  Matrix KF_R; // measurement covariance

};

#endif // VISUALODOMETRY_H
