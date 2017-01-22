function [xhat_data,Smat]=dd2(kalmfilex,kalmfiley,xbar,P0,q,r,u,y,timeidx,optpar)
% DD2
%   This function performs a DD2-filtering; a state estimation for nonlinear 
%   systems that is based on second-order polynomial approximations of the 
%   nonlinear mappings. The approximations are derived by using a 
%   multidimensional extension of Stirling's interpolation formula. 
%   The model of the nonlinear system must be specified in the form:
%               x(k+1) = f[x(k),u(k),v(k)]
%               y(k)   = g[x(k),w(k)]
%   where 'x' is the state vector, 'u' is a possible input, and 'v' and 'w'
%   are (white) noise sources.
%
% Call
%   [xhat,Smat]=dd2(xfile,yfile,x0,P0,q,r,u,y,timeidx,optpar) 
%
% Input
%   xfile   - File containing the state equations.
%   yfile   - File containing the output equations.
%   x0      - Initial state vector.
%   P0      - Initial covariance matrix (symmetric, nonnegative definite).
%   q,r     - Covariance matrices for v and w, respectively.
%   u       - Input signal. Dimension is [samples x inputs].
%             Use [] if there are no inputs.
%   y       - Output signal. Dimension is [observations x outputs].
%   timeidx - Vector containing sample numbers for the availability of
%             the observations in 'y'. The vector has same length as 'y'.
%   optpar  - Data structure containing optional parameters:
%             .A:     State transition matrix.
%             .C:     Output sensitivity matrix.
%             .F:     Process noise coupling matrix.
%             .G:     Measurement noise coupling matrix.
%             .init : Initial parameters for 'xfile' and 'yfile'
%                     (use an arbitrary format).
%
% Output
%   xhat    - State estimates. Dimension is [samples+1 x states].
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
% Written by: Magnus Norgaard
% LastEditDate: Nov. 9, 2001 

% >>>>>>>>>>>>>>>>>>>>>>>>>>> INITIALIZATIONS <<<<<<<<<<<<<<<<<<<<<<<<<<
h2    = 3;                 % Squared divided-difference step size
h     = sqrt(h2);          % Divided difference step-size
scal1 = 0.5/h;             % A convenient scaling factor
scal2 = sqrt((h2-1)/(4*h2*h2)); % Another scaling factor
if isempty(u),             % No inputs
  nu = 0; samples = timeidx(end); uk1 = [];
else
  [samples,nu] = size(u);  % # of samples and inputs
end
nx           = size(P0,1); % # of states
if isempty(xbar),          % Set to x0=0 if not specified 
  xbar = zeros(nx,1);
elseif length(xbar)~=nx,
  error('Dimension mismatch between x0 and P0');
end
ny    = size(y,2);         % # of outputs
nv    = size(q,1);         % # of process noise sources
nw    = size(r,1);         % # of measurement noise sources
[v,d] = eig(P0);           % Square root of initial state covariance
Sxbar = triag(real(v*sqrt(d)));
[v,d] = eig(q);            % Square root of process noise covariance
Sv    = real(v*sqrt(d));
hSv   = h*Sv;
[v,d] = eig(r);
Sw    = real(v*sqrt(d));   % Square root of measurement noise cov.
hSw   = h*Sw;
SxxSxv = zeros(nx,2*(nx+nv));  % Allocate compund matrix consisting of Sxx and Syv 
SyxSyw = zeros(ny,2*(nx+nw));  % Allocate compund matrix consisting of Syx and Syw
xhat_data = zeros(samples+1,nx); % Matrix for storing state estimates
Smat      = zeros(samples+1,0.5*nx*(nx+1)); % Matrix for storing cov. matrices
[I,J]     = find(triu(reshape(1:nx*nx,nx,nx))'); % Index to elem. in Sx
sidx      = sub2ind([nx nx],J,I); 
yidx  = 1;                 % Index into y-vector 
vmean = zeros(nv,1);       % Mean of process noise
wmean = zeros(nw,1);       % Mean of measurement noise

% ----- Initialize state+output equations and linearization -----
if nargin<10,              % No optional parameters passed
   optpar = [];
end
if isfield(optpar,'init')  % Parameters for m-functions
   initpar = optpar.init;
else
   initpar = [];
end

Aflag = 0; Cflag = 0; Fflag = 0; Gflag = 0;
nxnv2 = nx+nv;
nxnw2 = nx+nw;
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
   SxxSxv(:,nx+1:nx+nv) = F*Sv;
   nxnv2 = nxnv2-nv;
   Fflag = 1;
end
if isfield(optpar,'C'),    % Deterministic observation model linear
   C = optpar.C;
   if(size(C,1)~=ny | size(C,2)~=nx)
      error('"optpar.C" has the wrong dimension');
   end
   nxnw2 = nxnw2-nx;
   Cflag = 1;
end
if isfield(optpar,'G'),    % Linear observation noise model
   G = optpar.G;
   if(size(G,1)~=ny | size(G,2)~=nw)
      error('"optpar.G" has the wrong dimension');
   end
   SyxSyw(:,nx+1:nx+nw) = G*Sw;
   nxnw2 = nxnw2-nw;
   Gflag = 1;
end

% Index to location of Sxv2 and Syw2 in SxxSxv and SyxSyw matrices
if Cflag,
	idx_syw2 = nx+nw;
else
	idx_syw2 = 2*nx+nw;
end
if Aflag,
	idx_sxv2 = nx+nv;
else
	idx_sxv2 = 2*nx+nv;
end
SxxSxv = [SxxSxv zeros(nx,nxnv2)];
SyxSyw = [SyxSyw zeros(ny,nxnw2)];
feval(kalmfilex,initpar);  % State equation
feval(kalmfiley,initpar);  % Output equation
counter = 0;               % Counts the progress of the filtering
waithandle=waitbar(0,'Filtering in progress');  % Initialize waitbar


% >>>>>>>>>>>>>>>>>>>>>>>>>>>> FILTERING <<<<<<<<<<<<<<<<<<<<<<<<<<<
for k=0:samples,

  % --- Measurement update (a posteriori update) ---
  y0 = feval(kalmfiley,xbar,wmean);
  if (k<=timeidx(end) & timeidx(yidx)==k),
    ybar = ((h2-nxnw2)/h2)*y0;
    if Cflag,
       SyxSyw(:,1:nx) = C*Sxbar;
    else
       kx2 = nx+nw;
       for kx=1:nx,
          syp = feval(kalmfiley,xbar+h*Sxbar(:,kx),wmean);
          sym = feval(kalmfiley,xbar-h*Sxbar(:,kx),wmean);
          SyxSyw(:,kx)  = scal1*(syp-sym);
          SyxSyw(:,kx2+kx) = scal2*(syp+sym-2*y0);
          ybar = ybar + (syp+sym)/(2*h2);    
       end
    end
    if ~Gflag,
       for kw=1:nw,
          swp = feval(kalmfiley,xbar,hSw(:,kw));
          swm = feval(kalmfiley,xbar,-hSw(:,kw));
          SyxSyw(:,nx+kw)       = scal1*(swp-swm);
          SyxSyw(:,idx_syw2+kw) = scal2*(swp+swm-2*y0);
          ybar = ybar + (swp+swm)/(2*h2);          
       end
    end
    
    % Cholesky factor of a'posteriori output estimation error covariance
    Sy   = triag(SyxSyw);
    K    = (Sxbar*SyxSyw(:,1:nx)')/(Sy*Sy');
    xhat = xbar + K*(y(yidx,:)'-ybar);  % State estimate

    % Cholesky factor of a'posteriori estimation error covariance
    Sx   = triag([Sxbar-K*SyxSyw(:,1:nx) K*SyxSyw(:,nx+1:end)]);
    yidx = yidx + 1; 

  % No observations available at this sampling instant
  else
    xhat = xbar;                       % Copy a priori state estimate
    Sx   = Sxbar;                      % Copy a priori covariance factor
  end

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
          SxxSxv(:,kx) = scal1*(sxp-sxm);
          SxxSxv(:,kx2+kx) = scal2*(sxp+sxm-2*fxbar);
          xbar            =  xbar + (sxp+sxm)/(2*h2);
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
