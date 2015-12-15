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

#include "mex.h"
#include <iostream>
#include <string.h>
#include "visualodometry.h"

using namespace std;

static VisualOdometry *VO;

void mexFunction (int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[]) {

  // read command
  char command[128];
  mxGetString(prhs[0],command,128);

  // init
  if (!strcmp(command,"init")) {
    
    // check arguments
    if(nrhs!=1+4 && nrhs!=1+5)
      mexErrMsgTxt("4 or 5 parameters required (f,cu,cv,base/cam_height,[pitch]).");
    if(!mxIsDouble(prhs[1]) || mxGetM(prhs[1])*mxGetN(prhs[1])!=1)
      mexErrMsgTxt("Input f must be a double scalar.");
    if(!mxIsDouble(prhs[2]) || mxGetM(prhs[2])*mxGetN(prhs[2])!=1)
      mexErrMsgTxt("Input cu must be a double scalar.");
    if(!mxIsDouble(prhs[3]) || mxGetM(prhs[3])*mxGetN(prhs[3])!=1)
      mexErrMsgTxt("Input cv must be a double scalar.");
    if(!mxIsDouble(prhs[4]) || mxGetM(prhs[4])*mxGetN(prhs[4])!=1)
      mexErrMsgTxt("Input base/cam_height must be a double scalar.");
    if(nrhs==1+5)
      if(!mxIsDouble(prhs[5]) || mxGetM(prhs[5])*mxGetN(prhs[5])!=1)
        mexErrMsgTxt("Input pitch must be a double scalar.");
    
    double f     = *((double*)mxGetPr(prhs[1]));
    double cu    = *((double*)mxGetPr(prhs[2]));
    double cv    = *((double*)mxGetPr(prhs[3]));
    double base  = *((double*)mxGetPr(prhs[4]));
    double pitch = 0;
    if(nrhs==1+5)
      pitch = *((double*)mxGetPr(prhs[5]));
    
    VO = new VisualOdometry();
    VO->setCalibration(f,cu,cv,base,pitch);

  // close
  } else if (!strcmp(command,"close")) {
    delete VO;

  // update via observations
  } else if (!strcmp(command,"update")) {
    
    // check arguments
    if(nrhs!=1+2)
      mexErrMsgTxt("2 inputs required (deltaT,p_matched).");
    if(!mxIsDouble(prhs[1]) || mxGetM(prhs[1])*mxGetN(prhs[1])!=1)
      mexErrMsgTxt("Input deltaT must be a double scalar.");
    if(!mxIsDouble(prhs[2]) || (mxGetM(prhs[2])!=4 && mxGetM(prhs[2])!=8))
      mexErrMsgTxt("Input p_matched must be a 4xN or 8xN double matrix.");
    if(nlhs!=1)
      mexErrMsgTxt("1 output required (success).");
    
    // get pointers
    double  deltaT    = *((double*)mxGetPr(prhs[1]));
    double* p_matched =   (double*)mxGetPr(prhs[2]);
    int     M         =             mxGetM(prhs[2]);
    int     N         =             mxGetN(prhs[2]);
    
    // copy matches
    vector<Matcher::p_match> matches;    
    Matcher::p_match match;
    int32_t k=0;
    
    // monocular version
    if (M==4) {
      for (int i=0; i<N; i++) {
        match.u1p = p_matched[k++];
        match.v1p = p_matched[k++];
        match.u1c = p_matched[k++];
        match.v1c = p_matched[k++];
        matches.push_back(match);
      }
      
    // stereo version
    } else {
      for (int i=0; i<N; i++) {
        match.u1p = p_matched[k++];
        match.v1p = p_matched[k++];
        match.u2p = p_matched[k++];
        match.v2p = p_matched[k++];
        match.u1c = p_matched[k++];
        match.v1c = p_matched[k++];
        match.u2c = p_matched[k++];
        match.v2c = p_matched[k++];
        matches.push_back(match);
      }
    }
    
    // process
    const int dims[] = {1,1};
    plhs[0] = mxCreateNumericArray(2,dims,mxDOUBLE_CLASS,mxREAL);
    *((double*)mxGetPr(plhs[0])) = VO->update (matches,deltaT,M==8);

  // read transformation
  } else if (!strcmp(command,"gettransformation")) {
    
    // check arguments
    if(nlhs!=1)
      mexErrMsgTxt("1 output required (transformation H).");
    
    // get state
    Matrix H = ~VO->getTransformation();
    
    // create outputs
    const int dims[] = {4,4};
    plhs[0] = mxCreateNumericArray(2,dims,mxDOUBLE_CLASS,mxREAL);
    H.getData((double*)mxGetPr(plhs[0]));

  // get inliers
  } else if (!strcmp(command,"getinliers")) {
    
    // check arguments
    if(nlhs!=1)
      mexErrMsgTxt("1 output required (inliers).");
    
    // get inliers
    std::vector<int32_t> inliers = VO->getInliers();
    const int dims[] = {1,inliers.size()};
    plhs[0] = mxCreateNumericArray(2,dims,mxDOUBLE_CLASS,mxREAL);
    double *inliers_mex = (double*)mxGetPr(plhs[0]);
    for (int i=0; i<inliers.size(); i++)
      inliers_mex[i] = inliers[i]+1; // convert to matlab index (start with 1)
  
  // unknown command
  } else {
    mexPrintf("Unknown command: %s\n",command);
  }
  
}
