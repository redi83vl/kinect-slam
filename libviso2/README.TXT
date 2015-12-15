####################################################################################
# Copyright 2011. All rights reserved.                                             #
# Institute of Measurement and Control Systems                                     #
# Karlsruhe Institute of Technology, Germany                                       #
#                                                                                  #
# This file is part of libviso2.                                                   #
# Authors:  Andreas Geiger                                                         #
#           Please send any bugreports to geiger@kit.edu                           #
#                                                                                  #
# libviso2 is free software; you can redistribute it and/or modify it under the    #
# terms of the GNU General Public License as published by the Free Software        #
# Foundation; either version 2 of the License, or any later version.               #
#                                                                                  #
# libviso2 is distributed in the hope that it will be useful, but WITHOUT ANY      #
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A  #
# PARTICULAR PURPOSE. See the GNU General Public License for more details.         #
#                                                                                  #
# You should have received a copy of the GNU General Public License along with     #
# libviso2; if not, write to the Free Software Foundation, Inc., 51 Franklin       #
# Street, Fifth Floor, Boston, MA 02110-1301, USA                                  #
####################################################################################

+++++++++++++++++++++++++++++++++++
+          INTRODUCTION           +
+++++++++++++++++++++++++++++++++++

Libviso 2 (LIBrary for ViSual Odometry 2) is a cross-platfrom (Linux, Windows) C++
library with MATLAB wrappers for computing the 6 DOF motion of a moving stereo
camera. Input is a sequence of rectified stereo images. Output is a 4x4 matrix
which projects a point from the previous to the current camera coordinates.

Version 2 now also supports visual odometry from monocular sequences, and the
feature matching functions implement sparse stereo, scene flow and optical flow.
For more details, please have a look at the corresponding publication.

If you distribute a software that uses libviso, you have to distribute it under GPL
with the source code. Another option is to contact us to purchase a commercial license.

If you find this software useful or if you use this software for your research,
we would be happy if you cite the following related publication:

@INPROCEEDINGS{Geiger11,
 author = {Andreas Geiger and Julius Ziegler and Christoph Stiller},
 title = {StereoScan: Dense 3d Reconstruction in Real-time},
 booktitle = {IEEE Intelligent Vehicles Symposium},
 year = {2011},
 month = {June},
 address = {Baden-Baden, Germany}
}

+++++++++++++++++++++++++++++++++++
+    COMPILING MATLAB WRAPPERS    +
+++++++++++++++++++++++++++++++++++

Prerequisites needed for compiling libviso2:
- Matlab (if you want to compile the matlab wrappers)
- gcc or Visual Studio (if you want to include the code in your C++ project)

If you want to use libviso directly from MATLAB you can easily do this by using
the MATLAB wrappers provided. They also include some demo files for testing.

In the MATLAB directory of libviso, simply run 'make.m' to generate the mex wrappers.
(Run mex -setup before to choose your desired compiler)

Now try to run the demo*.m files!
For some of them you will need the Karlsruhe dataset from www.cvlibs.net.

+++++++++++++++++++++++++++++++++++
+     BUILDING A C++ LIBRARY      +
+++++++++++++++++++++++++++++++++++

You can directly use the headers from the 'src' directory in your C++ project.

For more information on the usage of this library, have a look into the MATLAB wrappers and
the documentation of the header files (especially matcher.h and visual visualodometry.h).

Please send any feedback and bugreports to geiger@kit.edu,
Andreas Geiger

