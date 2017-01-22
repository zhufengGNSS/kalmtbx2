function [xhat_data,Pmat]=ekf(kalmfilex,kalmfiley,linfile,xbar,...
                P0,q,r,u,y,timeidx,optpar)
% EKF
%  This function is an implementation of the conventional
%  extended Kalman filter (EKF).
%  The filter estimates the states for nonlinear systems written in
%  the general form:
%               x(k+1) = f[x(k),u(k),v(k)]
%               y(k)   = g[x(k),w(k)]
%
% Call: [xhat,Pmat]=ekf(xfunc,yfunc,linfunc,x0,P0,q,r,u,y,tidx,optpar) 
%
% Input:
%   xfunc   - Function containing the state equations.
%   yfunc   - Function containing the output equations.
%   linfunc - Function containing linearization procedure.
%   x0      - Initial state vector.
%   P0      - Initial covariance matrix (symmetric, nonnegative definite).
%   q,r     - Covariance matrices for v and w, respectively.
%   u       - Input signal. Dimension is [samples x inputs].
%             Use [] if there are no inputs.
%   y       - Output signal. Dimension is [observations x outputs].
%   tidx    - Vector containing time stamps (in samples) for the 
%             observations in y.
%   optpar  - Data structure containing optional parameters:
%             .init : Initial parameters for 'xfile', 'yfile', and
%                     'linfile' (use an arbitrary format).
%
% Output:
%   xhat    - State estimates. Dimension is [samples+1 x states].
%   Pmat    - Matrix where each row contains the upper triangular elements
%             of the covariance matrix estimates. The dimension is 
%             [samples+1 x 0.5*states*(states+1)]. The individual covariance 
%             matrices can later be extracted with MAT2COV.
%
% The user must write the three functions 'xfunc', 'yfunc', and 'linfunc' 
% containing state update, output equation, and linearization. The 
% function containing the state update should have the header 
% (the function name is arbitrary):
%       function x=my_xfile(x,u,v)
%
% the function containing the output equation must have the header
%       function y=my_yfile(x,w)
%
% while the function containing the linearization must have the header
%      function [M,N]=my_linfile(x,u,vw,flag)
% flag=0: Linearization of the state equation
% flag=1: Linerization of the output equation.
%  
% In all three cases, an initialization of constant parameters can be 
% made using the parameter 'optpar.init'. This parameter is passed through
% x if the functions are called with only one parameter.
%
% Written by Magnus Norgaard
% LastEditDate: Nov. 9, 2001 

% >>>>>>>>>>>>>>>>>>>>>>>>>>> INITIALIZATIONS <<<<<<<<<<<<<<<<<<<<<<<<<<
if isempty(u),             % No inputs
  nu = 0; samples = timeidx(end); uk1 = [];
else
  [samples,nu] = size(u);  % # of samples and inputs
end
ny           = size(y,2);  % # of outputs
nx           = size(P0,1); % # of states
nv           = size(q,1);  % # of process noise sources
nw           = size(r,1);  % # of measurement noise sources
if isempty(xbar),          % Set to x0=0 if not specified 
  xbar=zeros(nx,1);
elseif length(xbar)~=nx,
  error('Dimension mismatch between x0 and P0');
end
if size(y,1)~=size(timeidx,1)
  error('Dimension mismatch between y and timeidx');
end
Pxbar = P0;                % A priori estimate = initial covariance
xhat_data = zeros(samples+1,nx); % Matrix for storing state estimates
Pmat      = zeros(samples+1,0.5*nx*(nx+1)); % Matrix for storing cov. matrices
pidx      = find(tril(reshape(1:nx*nx,nx,nx))); % Index in P
yidx  = 1;                 % Index into y-vector 


% ----- Initialize state+output equations and linearization -----
if nargin<11,              % No optional parameters passed
   optpar = [];
end
if isfield(optpar,'init')  % Parameters for m-functions
   initpar = optpar.init;
else
   initpar = [];
end
vmean = zeros(nv,1);       % Mean of process noise
wmean = zeros(nw,1);       % Mean of measurement noise
feval(kalmfilex,initpar);  % Initialize state equation
feval(kalmfiley,initpar);  % Initialize output equation
feval(linfile,initpar);    % Initialize linearization
counter = 0;               % Counts the progress of the filtering
waithandle=waitbar(0,'Filtering in progress');  % Initialize waitbar

% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> FILTERING <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
for k=0:samples,

  % --- Measurement update (a posteriori update) ---
  ybar = feval(kalmfiley,xbar,wmean);
  if (k<=timeidx(end) & timeidx(yidx)==k),
    [C,G] = feval(linfile,xbar,[],wmean,1); % Linearization
    if isempty(G),                      % Kalman gain
       K = Pxbar*C'/(C*Pxbar*C'+r);     % Noise enters directly
    else
       K = Pxbar*C'/(C*Pxbar*C'+G*r*G');% General update
    end
    Px   = Pxbar-K*C*Pxbar;             % A posteriori covariance
    xhat = xbar + K*(y(yidx,:)'-ybar);  % State estimate
    yidx = yidx + 1;                    % Update index in time vector
  
  % no observations available at this sampling time
  else
    xhat = xbar;                        % Copy a priori state estimate
    Px   = Pxbar;                       % Copy a priori covariance factor
  end

  % --- Time update (a'priori update) of state and covariance ---
  if k<samples,
    if nu>0 uk1 = u(k+1,:)'; end
    xbar=feval(kalmfilex,xhat,uk1,vmean);    % State update
    [A,F] = feval(linfile,xhat,uk1,vmean,0); % Linearization
    if isempty(F),                     % Covariance update
       Pxbar = A*Px*A' + q;            % Noise enters directly
    else
       Pxbar = A*Px*A' + F*q*F';       % General update
    end
  end
  
  % --- Store results ---
  xhat_data(k+1,:) = xhat';
  Pmat(k+1,:)      = Px(pidx)';
  
  % --- How much longer? ---
  if (counter+0.01<= k/samples),
     counter = k/samples;
     waitbar(k/samples,waithandle);
  end
end
close(waithandle);
