function [xhat_data,Smat]=dd2m(kalmfilex,kalmfiley,xbar,P0,q,r,u,y,timeidx,optpar)
% DD2M
%   This function performs a DD2-filtering; a state estimation for
%   nonlinear systems that is based on second-order polynomial approximations
%   of the nonlinear mappings. The approximations are derived by using a 
%   multidimensional extension of Stirling's interpolation formula. 
%   The function is implemented to handle multiple observation streams.
%   The model of the nonlinear system must be specified in the form:
%               x(k+1) = f[x(k),u(k),v(k)]
%               y1(k)  = g1[x(k),w1(k)]
%                        :
%               yn(k)  = gn[x(k),wn(k)]
%   where 'x' is the state vector, 'u' is a possible input, and 'v' and 'w'
%   are (white) noise sources.
%
% Call
%   [xhat,Smat]=dd2m(xfile,yfile,x0,P0,q,r,u,y,tidx,optpar) 
%
% Input
%   xfile   - File containing the state equations.
%   yfunc   - Cell array specifying the names of the functions 
%             containing the output equations.
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
%             .A:     State transition matrix.
%             .C:     Output sensitivity matrix (cell array).
%             .F:     Process noise coupling matrix.
%             .G:     Measurement noise coupling matrix (cell array).
%             .init : Initial parameters for 'xfile', 'yfile'
%                     (use an arbitrary format).
%
% Output
%   xhat    - State estimates. The dimension is [samples+1 x states].
%   Smat    - Matrix where each row contains elements of (the upper triangular
%             part of) the Cholesky factor of the covariance matrix. The 
%             dimension is [samples+1 x 0.5*states*(states+1)]. The individual
%             covariance matrices can later be extracted with SMAT2COV.
%
%  The user must write the two m-functions 'xfile' and 'yfile' containing the
%  state update and the output equation. The function containing the state
%  update should take three arguments:
%       function x=my_xfile(x,u,v)
%
%  while the function containing the output equation should take two
%  arguments:
%       function y=my_yfile(x,w)
%
%  In both cases, an initialization of constant parameters can be 
%  made using the parameter 'optpar.init'. This parameter is passed through
%  x if the functions are called with only one parameter.
%
%  Literature:
%     M. Norgaard, N.K. Poulsen, O. Ravn: "New Developments in State
%     Estimation for Nonlinear Systems", Automatica, (36:11), Nov. 2000,
%     pp. 1627-1638.
%
% Written by Magnus Norgaard
% LastEditDate: Nov. 9, 2001 

% >>>>>>>>>>>>>>>>>>>>>>>>>>>> INITIALIZATIONS <<<<<<<<<<<<<<<<<<<<<<<<<<<
h2    = 3;                 % Squared divided difference step
h     = sqrt(h2);          % Divided difference step
scal1 = 0.5/h;             % Scaling factor
scal2 = sqrt((h2-1)/(4*h2*h2)); % Another scaling factor
nx    = size(P0,1); % # of states
nv    = size(q,1);  % # of process noise sources
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
  [v,d]        = eig(r{n});
  obs(n).nw    = size(r{n},1);
  obs(n).nxnw2 = nx+obs(n).nw;
  obs(n).Sw    = real(v*sqrt(d));  % Cholesky factor of measurement noise cov.
  obs(n).hSw   = h*obs(n).Sw;
  if (obs(n).nobs~=length(obs(n).tidx)),
    error('Dimension mismatch between y and tidx');
  end
  ny = ny + obs(n).ny;
  if obs(n).tidx(end)>lastsample,
     lastsample=obs(n).tidx(end);
  end
  idx1(n) = ny - obs(n).ny + 1;
  idx2(n) = ny;
  obs(n).Cflag=0;
  obs(n).Gflag=0;
end
if isempty(u),             % No inputs
  nu = 0; samples = lastsample; uk1 = [];
else
  [samples,nu] = size(u);  % # of samples and inputs
end

xhat_data = zeros(samples+1,nx); % Matrix for storing state estimates
Smat      = zeros(samples+1,0.5*nx*(nx+1)); % Matrix for storing cov. matrices
[I,J]     = find(triu(reshape(1:nx*nx,nx,nx))'); % Index to elem. in Sx
sidx      = sub2ind([nx nx],J,I); 
ybar      = zeros(ny,1);
y0        = zeros(ny,1);
yidx      = ones(streams,1);% Index into y-vectors 


% ----- Initialize state+output equations and linearization -----
if nargin<10,              % No optional parameters passed
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

[v,d] = eig(P0);           % Cholesky factor of initial state covariance
Sxbar = triag(real(v*sqrt(d)));
[v,d] = eig(q);            % Cholesky factor of process noise covariance
Sv    = real(v*sqrt(d));
hSv   = h*Sv;
Aflag = 0; Fflag = 0;
nxnv2 = nx+nv;
if isfield(optpar,'A'),    % Deterministic dynamic model is linear
   A = optpar.A;
   if(size(A,1)~=nx | size(A,2)~=nx)
      error('"optpar.A" has the wrong dimension');
   end
	nxnv2 = nxnv2-nx;
   Aflag = 1;
end
if isfield(optpar,'F'),    % Linear process noise model in state equation
   F = optpar.F;
   if(size(F,1)~=nx | size(F,2)~=nv)
      error('"optpar.F" has the wrong dimension');
   end
	nxnv2 = nxnv2-nv;
   Fflag = 1;
end
SxxSxv = zeros(nx,nx+nv+nxnv2);  % Allocate compund matrix consisting of
                                 % Sxx, Sxv and Sxx2, Sxv2 
if Fflag,
   SxxSxv(:,nx+1:nx+nv) = F*Sv;
end

% Index to location of Sxv2 in SxxSxv matrix
if Aflag,
   idx_sxv2 = nx+nv;
else
   idx_sxv2 = 2*nx+nv;
end


if isfield(optpar,'C'),    % Deterministic observation model is linear
   if ~iscell(optpar.C),
      error('"optpar.C" must be a cell array');
   elseif streams~=length(optpar.C),
      error('Number of cells in "optpar.C" is wrong');
   end
   for n=1:streams,
		if ~isempty(optpar.C{n})  % Observation model "n" is linear
         obs(n).Cflag = 1;
			obs(n).nxnw2 = obs(n).nxnw2-nx;
			obs(n).C = optpar.C{n};
		   if(size(obs(n).C,1)~=obs(n).ny | size(obs(n).C,2)~=nx)
				errstr=sprintf('optpar.C{%d} has the wrong dimension',n);
            error(errstr);
			end
      end
   end
end
if isfield(optpar,'G'),    % Linear observation noise model
   if ~iscell(optpar.G),
      error('"optpar.G" must be a cell array');
   elseif streams~=length(optpar.G),
      error('Number of cells in "optpar.G" is wrong');
   end
   for n=1:streams,
		if ~isempty(optpar.G{n})  % Observation model "n" is linear
         obs(n).Gflag = 1;
			obs(n).nxnw2 = obs(n).nxnw2-obs(n).nw;
			obs(n).G = optpar.G{n};
		   if(size(obs(n).G,1)~=obs(n).ny | size(obs(n).G,2)~=obs(n).nw)
				errstr=sprintf('optpar.G{%d} has the wrong dimension',n);
            error(errstr);
			end
		end
	end
end

% Allocate space for Syx, Syw and Syx2, Syw2 and insert a possible
% constant term.
for n=1:streams,
   obs(n).SyxSyw  = zeros(obs(n).ny,nx+obs(n).nw+obs(n).nxnw2);
   if obs(n).Gflag,
      obs(n).SyxSyw(:,nx+1:nx+obs(n).nw) = obs(n).G*obs(n).Sw;
	end
	
	% Index to location of Syw2 in SyxSyw matrix
   if obs(n).Cflag,
	   obs(n).idx_syw2 = nx+obs(n).nw;
   else
	   obs(n).idx_syw2 = 2*nx+obs(n).nw;
   end
end

feval(kalmfilex,initpar);       % Initialize state equation
for n=1:streams,
   feval(obs(n).yfunc,initpar); % Initialize output equations
end
counter = 0;               % Counts the progress of the filtering
waithandle=waitbar(0,'Filtering in progress');  % Initialize waitbar


% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> FILTERING <<<<<<<<<<<<<<<<<<<<<<<<<<<<<
for k=0:samples,

  % --- Measurement update (a posteriori update) ---
  for n=1:streams,
    y0(idx1(n):idx2(n)) = feval(obs(n).yfunc,xbar,obs(n).wmean);
    ybar(idx1(n):idx2(n)) = ((h2-obs(n).nxnw2)/h2)*y0(idx1(n):idx2(n));
    if (k<=obs(n).tidx(end) & obs(n).tidx(yidx(n))==k),
		 if obs(n).Cflag,
          obs(n).SyxSyw(:,1:nx) = obs(n).C*Sxbar;
       else
          kx2 = nx+obs(n).nw;
          for kx=1:nx,
             syp = feval(obs(n).yfunc,xbar+h*Sxbar(:,kx),obs(n).wmean);
             sym = feval(obs(n).yfunc,xbar-h*Sxbar(:,kx),obs(n).wmean);
             obs(n).SyxSyw(:,kx)    = scal1*(syp-sym);
	          obs(n).SyxSyw(:,kx2+kx)= scal2*(syp+sym-2*y0(idx1(n):idx2(n)));
	          ybar(idx1(n):idx2(n))  = ybar(idx1(n):idx2(n)) + (syp+sym)/(2*h2);    
          end
		 end
       if ~obs(n).Gflag,
			 kw2 = obs(n).idx_syw2;
          for kw=1:obs(n).nw,
             swp = feval(obs(n).yfunc,xbar,obs(n).hSw(:,kw));
             swm = feval(obs(n).yfunc,xbar,-obs(n).hSw(:,kw));
             obs(n).SyxSyw(:,nx+kw)  = scal1*(swp-swm);
	          obs(n).SyxSyw(:,kw2+kw) = scal2*(swp+swm-2*y0(idx1(n):idx2(n)));
	          ybar(idx1(n):idx2(n))   = ybar(idx1(n):idx2(n)) + (swp+swm)/(2*h2);
			 end
      end
      
      % Cholesky factor of a'posteriori output estimation error covariance
      Sy   = triag(obs(n).SyxSyw);

      % Kalman gain
      K    = (Sxbar*obs(n).SyxSyw(:,1:nx)')/(Sy*Sy');
      
      % State estimate
      xbar = xbar + K*[obs(n).y(yidx(n),:)'-ybar(idx1(n):idx2(n))];
      
      % Cholesky factor of a'posteriori estimation error covariance
      Sxbar= triag([Sxbar-K*obs(n).SyxSyw(:,1:nx) K*obs(n).SyxSyw(:,nx+1:end)]);

      yidx(n) = yidx(n) + 1;              % Update index in time vector
    end
  end
  xhat = xbar;
  Sx   = Sxbar;


  % --- Time update (a'priori update) of state and covariance ---
  if k<samples, 
    if nu>0 uk1 = u(k+1,:)'; end
    fxbar = feval(kalmfilex,xhat,uk1,vmean);
    xbar = ((h2-nxnv2)/h2)*fxbar;
    if Aflag,
        SxxSxv(:,1:nx) = A*Sx;
    else
       kx2 = nx+nv;
       for kx=1:nx,
          sxp = feval(kalmfilex,xhat+h*Sx(:,kx),uk1,vmean);
          sxm = feval(kalmfilex,xhat-h*Sx(:,kx),uk1,vmean);
          SxxSxv(:,kx)     = scal1*(sxp-sxm);
          SxxSxv(:,kx2+kx) = scal2*(sxp+sxm-2*fxbar);
          xbar             = xbar + (sxp+sxm)/(2*h2);
       end
    end
    if ~Fflag,
       for kv=1:nv,
          svp = feval(kalmfilex,xhat,uk1,hSv(:,kv));
          svm = feval(kalmfilex,xhat,uk1,-hSv(:,kv));
          SxxSxv(:,nx+kv)       = scal1*(svp-svm);
          SxxSxv(:,idx_sxv2+kv) = scal2*(svp+svm-2*fxbar);
          xbar                  = xbar + (svp+svm)/(2*h2);
       end
    end
    
    % Cholesky factor of a'priori estimation error covariance
    Sxbar = triag(SxxSxv);
  end
  

  % --- Store results ---
  xhat_data(k+1,:) = xhat';
  Smat(k+1,:)      = Sx(sidx)';
    
  % --- How much longer? ---
  if (counter+0.01<= k/samples),
     counter = k/samples;
     waitbar(k/samples,waithandle);
  end
end
close(waithandle);
