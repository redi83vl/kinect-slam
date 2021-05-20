# kinect_slam
Localization of moving Kinect RGBD Camera

RGBD SLAM consist in localizing an RGBD(red-green-blue-depth) sensor in an unknown environment and simultaneously build the map. Recently, research has focused on this technique because of the low cost MS Kinect sensor. With the Kinect sensor is possible to obtain dense 3D point cloud of indoor environment with RGB info. The same information may be obtained with stereo cameras by the quality of the stereo is not comparable with the Kinect.
The presented work demonstrate the efficency of the RGBD localization. The sensor is localized in the environment by using epipolar geoemtry to match consecutive RGB frames and the 3D depth map to compute the absolute orientation.
Image matching is done with a KLT like image matcher implemente in c++. Once thye image features are computed and matched then the 3D location is calculated. At this point we have two 3D point clouds associated. The 3D rigid transformation then connects the two point clouds gives us also the rigid transfrom describing sensor motion. Sensor motion is also filtered using UKF with a constant velocity assumption.
3D reconstruction(Mapping) is not implemented, this because of the real time contraint.

