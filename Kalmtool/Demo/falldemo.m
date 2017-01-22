% ====================== Falling body example ===========================
% Details about the benchmark example can be found in the paper
% "Suboptimal state estimation for continuous-time nonlinear systems from
% discrete noisy measurements", M. Athans, R.P. Wishner, A. Bertolini
% IEEE Trans. Automatic Control, vol. 13, no. 5, 1968 (pp. 504-514)
%
% Written by Magnus Norgaard
% LastEditDate: Dec. 11, 2001
Method = 5;                % Select filter method (1-5):
                           % 1=EKF, 2=DD1, 3=DD1 (mex-file), 4=DD2, 5=DD2 (mex)
xfunc = 'body1';           % File containing state equations
yfunc = 'body2';           % File containing output equations
linfunc = 'bodylin';       % File containing the linearization
x0 = [3e5;2e4;1e-3];       % Initial state vector
Q  = [zeros(3)];           % Covariance of process noise
r  = 1e4;                  % Covariance of measurement noise
P0 = diag([1e6 4e6 1e-4]); % Initial covariance on state estimate
gamma = 5e-5;              % Model parameter
M     = 1e5;               % Horizontal radar position (ft)
H     = 1e5;               % Vertical radar position (ft)
rksteps = 64;              % RK steps / sampling period
delta = 1/rksteps;         % Fast "Sampling period"
runs  = 50;                % Monte carlo repetitions


% ---- Generate test data ----
randn('seed',0);           % Set seed for random noise
Tfinal = 60;               % Simulate 'Tfinal' seconds
ysim   = zeros(Tfinal,1);  % Store true y sequence
xtrue  = [x0';zeros(Tfinal,3)];
xhatmat= zeros(Tfinal+1,3,runs);
v      = zeros(3,1);       % No process noise
w0     = 0;                % Mean of measurement noise
clear optpar
optpar.init = [delta M H gamma]; % Prepare initialization parameters
optpar.F = eye(3);
optpar.G=1;

% Run the simulation
x      = x0;
body1(optpar.init);              % Initialze state update
body2(optpar.init);              % Initialize observation equation
for k=1:Tfinal,
  for kk=1:1/delta,
     x=body1(x,[],v);
  end
  xtrue(k+1,:) = x';
  ysim(k)=body2(x,w0);
end

% Generate 'run' different data sets
ytrue = repmat(ysim,1,runs) + sqrt(r)*randn(Tfinal,runs);


%----- Do the Monte Carlo shit -----
x0hat = [x0(1:2);3e-5];        % Initial state estimate
idx   = [1:Tfinal]'*rksteps;   % Measurement time stamps (in rk-periods)

[v,d] = eig(P0);           % Cholesky factor of initial state covariance
Sx0 = real(v*sqrt(d));
[v,d] = eig(Q);            % Cholesky factor of process noise covariance
Sv    = real(v*sqrt(d));
[v,d] = eig(r);            % Cholesky factor of measurement noise covariance
Sw    = real(v*sqrt(d));

for k=1:runs,
fprintf('\nExperiment no. %d\n',k);

  %----- Estimate state trajectory -----
  switch Method
     case 1,
       [xhat,Pmat]=ekf(xfunc,yfunc,linfunc,x0hat,P0,Q,r,[],ytrue(:,k),idx,optpar);
     case 2,
       [xhat,Smat]=dd1(xfunc,yfunc,x0hat,P0,Q,r,[],ytrue(:,k),idx,optpar);
     case 3
       [xhat,Smat]=dd1fall(x0hat,Sx0,Sv,Sw,[],ytrue(:,k),idx,optpar);
     case 4,
       [xhat,Smat]=dd2(xfunc,yfunc,x0hat,P0,Q,r,[],ytrue(:,k),idx,optpar);
     case 5,
       [xhat,Smat]=dd2fall(x0hat,Sx0,Sv,Sw,[],ytrue(:,k),idx,optpar);
     otherwise
        error('No valid filter method selected. Method=1...5')
  end
  
  % ----- Store data -----
  xhatmat(:,:,k) = xhat([1;(idx+1)],:);
end
xhatmean = mean(xhatmat,3);    % Calculate mean values


%----- Display Results -----
close all
figure(1)
subplot(111)
plot(0:Tfinal,abs(xtrue(:,1)-xhatmean(:,1))); axis([0 60 0 300])
ylabel('Absolute value of average altitude error (ft)');
xlabel('Time (sec)');

figure(2)
plot(0:Tfinal,abs(xtrue(:,2)-xhatmean(:,2))); axis([0 60 0 350])
ylabel('Absolute value of average velocity error (ft)');
xlabel('Time (sec)');

figure(3)
semilogy(0:Tfinal,abs(xtrue(:,3)-xhatmean(:,3))); axis([0 60 1e-9 1e-2])
ylabel('Absolute value of average error in ballistic coefficient');
xlabel('Time (sec)');

