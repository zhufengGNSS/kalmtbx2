% ========================== AGV pose example ============================
% Position and orientation ("pose") of an AGV are estimated along with
% the radius of each of the two driving wheels and the distance between
% the wheels. Pose and model parameters are estimated in two different ways:
% First approach: Combine the "encoder model" with a pose detection from
% the roof-mounted camera.
% Second approach: Use also observations of two simple guidemarks. To 
% show the flexibility, only x and y estimates are provided in this
% observation stream.
%
% The demo shows how to call the various filters with multiple observation
% streams. It is not a particularly good demonstration of their performance
% as the problem is relatively linear in between samples.
%
% Written by Magnus Norgaard
% LastEditDate: Nov. 20, 2001

%=========================================================================
Method = 5;           % Select filter method (1-5):
                      % 1=EKF, 2=DD1, 3=DD1 (mex-file), 4=DD2, 5=DD2 (mex)
%=========================================================================

close all; clear optpar;
xfunc = 'agvtu';      % File containing state equations
yfunc = 'agvobs';     % File containing output equations
yfunc2= 'agvobs2';    % File containing output equations (2nd stream)
linfunc = 'agvlin';   % File containing the linearization
b0   = 0.5;           % Nominal distance between wheels (m)
rr0  = 0.1;           % Nominal radius of left wheel (m)
rl0  = 0.1;           % Nominal radius of right wheel (m)
P0= diag([1e-4 1e-4 1e-4 1e-6 1e-6 1e-6]); % Initial covariance
Q = diag([1/12 1/12]);     % Cov. of process noise
R = diag([1e-7 1e-7 1e-8]);           % Cov. of measurement noise
[v,d] = eig(P0);      % Cholesky factor of initial state covariance
Sx0 = real(v*sqrt(d));
[v,d] = eig(Q);       % Cholesky factor of process noise covariance
Sv    = real(v*sqrt(d));
[v,d] = eig(R);       % Cholesky factor of measurement noise covariance
Sw    = real(v*sqrt(d));
opt.C = [eye(3) zeros(3)];
opt.G = eye(3);



%----- Initializations -----
load agvexp           % Load observations from AGV experiment
x0 = [y(1,:) rr0 rl0 b0]';   % Initial state = first observation

  
%----- Run the filters using only observations from the roof camera -----
fprintf('\nExperiment with one stream of observations (roof camera).\n')
switch Method
   case 1,
     [xhat,Pmat]=ekf(xfunc,yfunc,linfunc,x0,P0,Q,R,u,y,tidx);
   case 2,
     [xhat,Smat]=dd1(xfunc,yfunc,x0,P0,Q,R,u,y,tidx,opt);
   case 3
     [xhat,Smat]=dd1agv(x0,Sx0,Sv,Sw,u,y,tidx,opt);
   case 4,
     [xhat,Smat]=dd2(xfunc,yfunc,x0,P0,Q,R,u,y,tidx,opt);
   case 5,
     [xhat,Smat]=dd2agv(x0,Sx0,Sv,Sw,u,y,tidx,opt);
   otherwise,
      error('No valid filter method selected. Method= 1...5');
end

%----- Now add the x,y-observations from the guide mark detection -----
fprintf('\nExperiment with two observation streams (roof camera and guide marks).\n')
R2 = 10*R(1:2,1:2);        % Cov. of measurement noise
[v,d] = eig(R2);           % Cholesky factor of measurement noise cov.
Sw2   = real(v*sqrt(d));
optm.C = {[eye(3) zeros(3)],[eye(2) zeros(2,4)]};
optm.G = {eye(3),eye(2)};
switch Method
   case 1,
     [xhat2,Pmat2]=ekfm(xfunc,{yfunc,yfunc2},linfunc,x0,P0,Q,{R,R2},...
                    u,{y,y2},{tidx,tidx2});
   case 2,
     [xhat2,Smat2]=dd1m(xfunc,{yfunc,yfunc2},x0,P0,Q,{R,R2},u,{y,y2},{tidx,tidx2},optm);
   case 3,
     [xhat2,Smat2]=dd1magv(x0,Sx0,Sv,{Sw,Sw2},u,{y,y2},{tidx,tidx2},optm);
   case 4,
     [xhat2,Smat2]=dd2m(xfunc,{yfunc,yfunc2},x0,P0,Q,{R,R2},u,{y,y2},{tidx,tidx2},optm);
   case 5,
     [xhat2,Smat2]=dd2magv(x0,Sx0,Sv,{Sw,Sw2},u,{y,y2},{tidx,tidx2},optm);
   otherwise,
      error('No valid filter method selected. Method= 1...5');
end


%----- Display the results -----
samples=size(u,1);
minval = min([xhat(:,[1:3]);y]);
maxval = max([xhat(:,[1:3]);y]);
figure(1)
subplot(311)
plot(0:samples,xhat(:,1),'-',tidx,y(:,1),'+')
grid
axis([0 samples minval(1) maxval(1)])
ylabel('x-coordinate')

subplot(312)
plot(0:samples,xhat(:,2),'-',tidx,y(:,2),'+')
grid
axis([0 samples minval(2) maxval(2)])
ylabel('y-coordinate')

subplot(313)
plot(0:samples,xhat(:,3),'-',tidx,y(:,3),'+')
grid
axis([0 samples minval(3) maxval(3)])
ylabel('theta')
xlabel('Time (samples)')

figure(2)
clf
subplot(211)
plot(0:samples,xhat(:,4:5))
grid
axis([0 samples min(min(xhat(:,4:5))) max(max(xhat(:,4:5)))])
ylabel('Right and left wheel radius')

subplot(212)
plot(0:samples,xhat(:,6))
grid
axis([0 samples min(xhat(:,6)) max(xhat(:,6))])
ylabel('Wheel distance')
xlabel('Time (samples)')
