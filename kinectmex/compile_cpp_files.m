% This function compile_cpp_files will compile the c++ code files
% which wraps OpenNI for the Kinect in Matlab.
%
% Please install first on your computer:
% - NITE-Bin-Win32-v1.3.0.18
% - OpenNI-Bin-Win32-v1.0.0.25
%
% Please change the Path inside this function if OpenNI is not 
% installed to C:\Program Files (x86)\OpenNI\
%
% Just execute by:
%
%   compile_c_files 
%

clear all; close all; clc

OpenNiPath='C:\Program Files (x86)\OpenNI\';
OpenNiPathLib='C:\Program Files (x86)\OpenNI\Lib';
OpenNiPathInclude='C:\Program Files (x86)\OpenNI\Include';

cd('Mex');
files=dir('*.cpp');
for i=1:length(files)
    Filename=files(i).name;
    clear(Filename); 
    mex('-v',['-L' OpenNiPathLib],'-lopenNI',['-I' OpenNiPathInclude],Filename);
end
cd('..');