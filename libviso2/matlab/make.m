dbclear all;

% compile matlab wrappers
disp('Building wrappers ...');
mex('matcherMex.cpp','../src/matcher.cpp','../src/filter.cpp','../src/triangle.cpp','../src/matrix.cpp','-I../src','CXXFLAGS=-msse3 -fPIC');
mex('visualOdometryMex.cpp','../src/visualodometry.cpp','../src/matrix.cpp','-I../src');
mex('reconstructionMex.cpp','../src/reconstruction.cpp','../src/matrix.cpp','-I../src');
disp('...done!');

