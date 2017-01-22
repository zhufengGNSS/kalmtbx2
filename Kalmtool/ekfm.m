function [xhat_data,Pmat]=ekfm(kalmfilex,kalmfiley,linfile,xbar,...
                P0,q,r,u,y,timeidx,optpar)
% EKFM
%  This function is an implementation of the conventional
%  extended Kalman filter (EKF).
%  It is implemented to handle multiple observation streams.
%  The filter estimates the states for nonlinear systems written in
%  the general form:
%               x(k+1) = f[x(k),u(k),v(k)]
%               y1(k)  = g1[x(k),w1(k)]
%                        :
%               yn(k)  = gn[x(k),wn(k)]
%   where 'x' is the state vector, 'u' is a possible input, and 'v' and 'w'
%   are (white) noise sources.
%
% Call:
%   [xhat,Pmat]=ekfm(xfunc,yfunc,linfunc,x0,P0,q,r,u,y,tidx,optpar) 
%
% Input:
%   xfunc   - Name of function containing the state equations.
%   yfunc   - Cell array specifying the names of the functions 
%             containing the output equations.
%   linfunc - Function containing linearization procedure.
%   x0      - Initial state vector.
%   P0      - Initial covariance matrix (symmetric, nonnegative definite).
%   q       - Covariance matrices for process noise.
%   r       - Cell array containing the measurement noise cov. matrices. 
%   u       - Input signal. Dimension is [samples x inputs].
%             Use [] if there are no inputs.
%   y       - Cell array containing the output signals. 
%             Dimension of each stream is [observations x outputs-in-stream].
%   tidx    - Cell array containing vector with time stamps (in samples) 
%             for the observations in y.
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
% the functions containing output equations must have the header
%       function y=my_yfile(x,w)
%
% while the function containing the linearization must have the header
%      function [M,N]=my_linfile(x,u,vw,flag)
% flag=0: Linearization of the state equation
% flag=i: Linerization of the output equation no. i (i=1...n).
%  
% In all three cases, an initialization of constant parameters can be 
% made using the parameter 'optpar.init'. This parameter is passed through
% x if the functions are called with only one parameter.
% 
% Written by Magnus Norgaard
% LastEditDate: Nov. 9, 2001

% >>>>>>>>>>>>>>>>>>>>>>>>>>> INITIALIZATIONS <<<<<<<<<<<<<<<<<<<<<<<<<<
nx           = size(P0,1); % # of states
nv           = size(q,1);  % # of process noise sources
if isempty(xbar),          % Set to x0=0 if not specified 
  xbar = zeros(nx,1);
elseif length(xbar)~=nx,
  error('Dimension mismatch between x0 and P0');
end
streams    = length(y);
if ~(iscell(kalmfiley) & iscell(r) & iscell(timeidx) & iscell(y))
  error('"yfunc", "r", "tidx", and "y" must be cell array');
elseif (streams~=length(r) | streams~=length(timeidx) | ...
                                 streams~=length(kalmfiley))
  error('"yfunc", "r", "tidx", and "y" must have same number of cells');
end
ny         = 0;                % Total number of observations
lastsample = 0;                % Number of sample containing last observation
idx1 = zeros(streams,1);       % Index to start of each stream in ybar
idx2 = zeros(streams,1);       % Index to end of each stream in ybar
for n=1:streams,               % Wrap information about observation stream 
  obs(n).yfunc = kalmfiley{n}; % into data structure
  obs(n).y     = y{n};
  obs(n).tidx  = timeidx{n};
  obs(n).ny    = size(obs(n).y,2);
  obs(n).nobs  = size(obs(n).y,1);
  obs(n).r     = r{n};
  obs(n).nw    = size(obs(n).r,1);
  if (obs(n).nobs~=size(obs(n).tidx,1)),
    error('Dimension mismatch between y and tidx');
  end
  ny = ny + obs(n).ny;
  if obs(n).tidx(end)>lastsample,
     lastsample=obs(n).tidx(end);
  end
  idx1(n) = ny - obs(n).ny + 1;
  idx2(n) = ny;
end
if isempty(u),             % No inputs
  nu = 0; samples = lastsample; uk1 = [];
else
  [samples,nu] = size(u);  % # of samples and inputs
end

Pxbar = P0;                % A priori estimate = initial covariance
xhat_data = zeros(samples+1,nx); % Matrix for storing state estimates
Pmat      = zeros(samples+1,0.5*nx*(nx+1)); % Matrix for storing cov. matrices
pidx      = find(tril(reshape(1:nx*nx,nx,nx))); % Index in P
ybar      = zeros(ny,1);
yidx  = ones(streams,1);   % Index into y-vectors 


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
for n=1:streams,           % Mean of measurement noise
   obs(n).wmean = zeros(obs(n).nw,1);
end

feval(kalmfilex,initpar);      % Initialize state equation
for n=1:streams,
   feval(obs(n).yfunc,initpar);% Initialize output equations
end
feval(linfile,initpar);        % Initialize linearization

counter = 0;                   % Counts the progress of the filtering session
waithandle=waitbar(0,'Filtering in progress'); % Initialize waitbar

% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> FILTERING <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
for k=0:samples,

  % --- Measurement update (a posteriori update) ---
  for n=1:streams,
    ybar(idx1(n):idx2(n)) = feval(obs(n).yfunc,xbar,obs(n).wmean);
    if (k<=obs(n).tidx(end) & obs(n).tidx(yidx(n))==k),
    
      % Linearization
      [C,G] = feval(linfile,xbar,[],obs(n).wmean,n);

      % Kalman gain
      if isempty(G),                      
         K = Pxbar*C'/(C*Pxbar*C'+obs(n).r); % Noise enters directly
      else
         K = Pxbar*C'/(C*Pxbar*C'+G*obs(n).r*G');% General update
      end

      % A posteriori covariance
      Pxbar = Pxbar-K*C*Pxbar;
      
      % State estimate
      xbar = xbar + K*[obs(n).y(yidx(n),:)'-ybar(idx1(n):idx2(n))];      
      yidx(n) = yidx(n) + 1;              % Update index in time vector
    end
  end
  xhat = xbar;
  Px   = Pxbar;
  
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
     waitbar(k/samples,waithandle);  % Update waitbar
  end
end
close(waithandle);                   % Close waitbar window
