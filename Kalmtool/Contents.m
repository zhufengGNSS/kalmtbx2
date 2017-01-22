%KALMTOOL version 2
%   A toolbox for state estimation for nonlinear systems.
%
%
%Filter functions
% dd1      - DD1-filter.
% dd1m     - DD1-filter for systems with multiple observation streams.
% dd2      - DD2-filter.
% dd2m     - DD2-filter for systems with multiple observation streams.
% ekf      - Extended Kalman filter.
% ekfm     - Extended Kalman filter for multiple observation streams.
%
%Mex files
% dd1c     - C-Mex version of 'dd1'.
% dd1mc    - C-Mex version of 'dd1mc'.
% dd2c     - C-Mex version of 'dd2'.
% dd2mc    - C-Mex version of 'dd2mc'.
% xytest   - Tests your own C-functions before the filtering is performed.
%
%Utilities
% kalmeval - Evaluate filter performance.
% mat2cov  - Extract covariance matrix from vector of upper triangular elements.
% mat2var  - Extract variance estimates from matrix of covariance estimates.
% smat2var - Calculate variance estimate for each state based on the Chol. fact..
% smat2cov - Make covariance matrix from vector of Cholesky factor elements.
% triag    - Triangularization with Householder transformation.
% kalmlb** - Object file that must be linked with the mex files.
% 
%Demonstrations (in Subdirectory 'Demo')
% agvdemo  - Position and orientation estimation of an AGV.
%            Simultaneous model calibration.
% falldemo - Falling body example (a continuous-discrete example).
% demomex  - Generate MEX-files to speed up the demonstration examples.
%            (Compiles the source code in dd1agv.c, dd1agvm.c, dd2agv.c,
%             dd2agvm.c, agvfct.c, agvfctm.c, dd2fall.c, dd2fall.c, fallfct.c).
%
% All functions are written by Magnus Norgaard
% November 2001
