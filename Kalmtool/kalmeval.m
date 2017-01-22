function [yhat_data,RMS]=kalmeval(method,kalmfiley,r,x_data,PS,y,tidx,optpar);
% KALMEVAL    Evaluate filter performance.
%
%    Assuming a (nonlinear) observation equation
%               y(k) = g[x(k),w(k)],
%    where 'x' is the state vector and 'w' is (white) observation noise,
%    KALMEVAL estimates the output (y) from the state estimates obtained
%    by filtering. The function plots the observed and estimated outputs
%    as well as the state estimates along with 3 times their standard
%    deviations.
%    
% Call
%   [yhat,RMS]=kalmeval('method','yfile',R,xhat,PS,y,tidx,optpar) 
%
% Input
%   method  - Filter method (ekf, dd1, dd2, ekfm, dd1m, dd2m).
%   yfile   - Name of file containing the observation equations.
%   R       - Covariance matrix for the measurement noise.
%             only used if method='dd2' or 'dd2m'.
%   xhat    - State estimates. Dimension is [samples+1 x states].
%   PS      - Matrix where each row contains elements of (the upper triangular
%             part of) the Cholesky factor of the covariance matrix
%             (dd1, dd2, dd1m, dd2m) or the covariance matrix (ekf, ekfm).
%             The dimension is [samples+1 x 0.5*states*(states+1)].
%   y       - Signal of observed outputs. Dimension: [observations x outputs].
%   tidx    - Vector containing time stamps (in samples) for the 
%             observations in y.
%   optpar  - (Optional) Data structure containing optional parameters:
%	           .G:     Measurement noise coupling matrix (/cell array).
%             .init : Initial parameters for 'xfile', 'yfile'
%                     (use an arbitrary format).
%
%   In case of multiple observation streams (dd1m, dd2m, ekfm), the
%   arguments yfile, R, y, and tidx must be cell arrays. 
%
% Output
%   yhat    - Output estimate. Dimension is [samples+1 x outputs].
%   RMS     - RMS error for each output.
%
%  The user must provide the m-functions 'yfile' containing the
%  output equation. The function should take two arguments:
%       function y=my_yfile(x,w)
%
%  Initialization of constant parameters can be made using the parameter 
%  'optpar.init'. This parameter is passed through x if the functions 
%  are called with only one parameter.
%
%
% Written by Magnus Norgaard
% LastEditDate: Nov. 22, 2001 

% >>>>>>>>>>>>>>>>>>>>>>>>>>> INITIALIZATIONS <<<<<<<<<<<<<<<<<<<<<<<<<<
h2    = 3;                 % Squared divided difference step
h     = sqrt(h2);          % Divided difference step
if nargin<8,               % No optional parameters passed
   optpar = [];
end
if isfield(optpar,'init')  % Parameters for m-functions
   initpar = optpar.init;
else
   initpar = [];
end


% ----- Check if there are muliple streams -----
if strcmpi(method,'dd1m') | strcmpi(method,'dd2m') | strcmpi(method,'ekfm'),
   streamflag = 1;
   streams    = length(y);
   if ~(iscell(kalmfiley) & iscell(r) & iscell(tidx) & iscell(y))
     error('"yfunc", "R", "tidx", and "y" must be cell structures');
   elseif (streams~=length(r) | streams~=length(tidx) | ...
                                 streams~=length(kalmfiley))
     error('"yfunc", "r", "tidx", and "y" must have same number of cells');
   end
   idx1 = zeros(streams,1);       % Index to start of each stream in yhat/ybar
   idx2 = zeros(streams,1);       % Index to end of each stream in yhat/ybar
   ny   = 0;
   for n=1:streams,               % Wrap information about observation stream 
     obs(n).yfunc = kalmfiley{n}; % into data structure
     obs(n).y     = y{n};
     obs(n).tidx  = tidx{n};
     obs(n).ny    = size(obs(n).y,2);
     obs(n).nobs  = size(obs(n).y,1);
     [v,d]        = eig(r{n});
     obs(n).nw    = size(r{n},1);  
     obs(n).hSw   = h*real(v*sqrt(d)); % Square root of measurement noise cov.
     if (obs(n).nobs~=size(obs(n).tidx,1)),
       error('Dimension mismatch between y and tidx');
     end
     ny = ny + obs(n).ny;
     idx1(n) = ny - obs(n).ny + 1;
     idx2(n) = ny;
	  obs(n).wmean = zeros(obs(n).nw,1); % Mean of measurement noise
	  obs(n).nxnw2 = nx+obs(n).nw;
	  obs(n).Gflag = 0;
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
			end
		end
	end
else
   streamflag = 0;
   ny    = size(y,2);         % # of outputs
   nw    = size(r,1);         % # of measurement noise sources
   [v,d] = eig(r);            % Square root of process noise covariance
   hSw   = h*real(v*sqrt(d)); % ..times h
   feval(kalmfiley,initpar);  % Output equation
   wmean = zeros(nw,1);       % Mean of measurement noise
	Gflag = 0;
   nxnw2 = nx+nw;
	if isfield(optpar,'G'),    % Linear observation noise model
      nxnw2 = nxnw2-nw;
      Gflag = 1;
   end
end

yhat      = zeros(ny,1);      % Allocate yhat 
[samples,nx]= size(x_data);   % # of samples and states
samples     = samples -1;
yhat_data = zeros(samples+1,ny); % Matrix for storing estimates of output
[I,J]     = find(triu(reshape(1:nx*nx,nx,nx))'); % Index to elem. in Sx
sidx      = sub2ind([nx nx],J,I);
Sx        = zeros(nx,nx);  % Allocate memory for covariance 'root'
counter   = 0;             % Counts the progress of the evaluation
waithandle=waitbar(0,'Evaluating filter');  % Initialize wait bar


% >>>>>>>>>>>>>>>>>>>>>>>>>>>> EVALUATE OUTPUT <<<<<<<<<<<<<<<<<<<<<<<<<<<
% ---------- DD1 and EKF ----------
if strcmpi(method,'dd1') | strcmpi(method,'ekf')
   for k=0:samples,

     % --- Output prediction ---  
     yhat = feval(kalmfiley,x_data(k+1,:)',wmean); % output estimate  
     yhat_data(k+1,:) = yhat';

     % --- How much longer? ---
     if (counter+0.01<= k/samples),
		  counter = k/samples;
		  waitbar(k/samples,waithandle);
     end
   end
   
% ---------- DD2 ----------
elseif strcmpi(method,'dd2')
  for k=0:samples,

     % --- Output prediction ---  
     Sx(sidx) = PS(k+1,:);                % Extract covariance estimate
     xhat = x_data(k+1,:)';               % Current state estimate
     yhat = feval(kalmfiley,x_data(k+1,:)',wmean); % "Nominal" output estimate
     yhat = ((h2-nxnw2)/h2)*yhat;
     for kx=1:nx,
       syp = feval(kalmfiley,xhat+h*Sx(:,kx),wmean);
       sym = feval(kalmfiley,xhat-h*Sx(:,kx),wmean);
       yhat = yhat + (syp+sym)/(2*h2);    
     end
	  if ~Gflag,
        for kw=1:nw,
           swp = feval(kalmfiley,xhat,hSw(:,kw));
           swm = feval(kalmfiley,xhat,-hSw(:,kw));
           yhat = yhat + (swp+swm)/(2*h2);          
        end
	  end
     yhat_data(k+1,:) = yhat';

     % --- How much longer? ---
     if (counter+0.01<= k/samples),
		  counter = k/samples;
		  waitbar(k/samples,waithandle);
     end
   end

% ---------- DD1M and EKFM ----------
elseif strcmpi(method,'dd1m') | strcmpi(method,'ekfm'),
   for k=0:samples,
      for n=1:streams,
			% --- Output prediction ---
			xhat =   x_data(k+1,:)';
			yhat(idx1(n):idx2(n)) = feval(obs(n).yfunc,xhat,obs(n).wmean);  
		end
      yhat_data(k+1,:) = yhat';
      
      % --- How much longer? ---
		if (counter+0.01<= k/samples),
			counter = k/samples;
			waitbar(k/samples,waithandle);
     end
   end


% ---------- DD2M ----------
elseif strcmpi(method,'dd2m'),
   for k=0:samples,
      Sx(sidx) = h*PS(k+1,:);              % Extract covariance estimate
      xhat = x_data(k+1,:)';               % Current state estimate

      for n=1:streams,
	      % --- Output prediction ---
			yhat(idx1(n):idx2(n)) = feval(obs(n).yfunc,xhat,obs(n).wmean);
         yhat(idx1(n):idx2(n)) = ((h2-obs(n).nxnw2)/h2)*yhat(idx1(n):idx2(n));
         for kx=1:nx,
            syp = feval(obs(n).yfunc,xhat+Sx(:,kx),obs(n).wmean);
            sym = feval(obs(n).yfunc,xhat-Sx(:,kx),obs(n).wmean);
			   yhat(idx1(n):idx2(n)) = yhat(idx1(n):idx2(n)) + (syp+sym)/(2*h2);
			end
			if ~obs(n).Gflag,
			   for kw=1:obs(n).nw,
				   swp = feval(obs(n).yfunc,xhat,obs(n).hSw(:,kw));
				   swm = feval(obs(n).yfunc,xhat,-obs(n).hSw(:,kw));
				   yhat(idx1(n):idx2(n)) = yhat(idx1(n):idx2(n)) + (swp+swm)/(2*h2);
			   end
			end
      end
      yhat_data(k+1,:) = yhat';
      
     % --- How much longer? ---
     if (counter+0.01<= k/samples),
		  counter = k/samples;
		  waitbar(k/samples,waithandle);
     end
   end

% ---------- Invalid Filter Specified ----------
else
   close(waithandle);                  % Close wait bar
   error(['"' method '" is not a valid filter method']);
end

% Calculate RMS error
if streamflag==0,
   RMS = sqrt(sum((y-yhat_data(tidx,:)).^2,1)/size(y,1));
else
   RMS = zeros(1,ny);
   for n=1:streams,
      RMS(idx1(n):idx2(n))= ...
        sqrt(sum((obs(n).y-yhat_data(obs(n).tidx,idx1(n):idx2(n))).^2,1)/obs(n).nobs);
   end
end

% Close wait bar
close(waithandle);


% >>>>>>>>>>>>>>>>>>>>>>>>>>> PLOT RESULTS <<<<<<<<<<<<<<<<<<<<<<<<<<
Maxsubplots = 3;             % Maximum # of plots in one figure window

% --- Plot observed and estimated output ---
clf
% ---------- Multiple observation streams ------
if streamflag==1,
   % ----- Plot observations -----
   for n=1:streams
     subplot(streams,1,n)
     plot(0:samples,yhat_data(:,idx1(n):idx2(n)),'-',obs(n).tidx,obs(n).y,'+')
     grid
     set(gca,'Xlim',[0 samples])
     ylabel(['Stream #' num2str(n)]);
     if n==1,
       title('Observed output (+) and estimated output (line)')
     end
     if n==streams,
       xlabel('Time (samples)')
     end
   end

   % ----- Plot residuals -----
   figure
   for n=1:streams
     subplot(streams,1,n)
     plot(obs(n).tidx,[obs(n).y-yhat_data(obs(n).tidx,idx1(n):idx2(n))])
     grid
     set(gca,'Xlim',[0 samples])
     ylabel(['Stream #' num2str(n)]);
     if n==1,
       title('Error between observed and estimated output')
     end
     if n==streams,
       xlabel('Time (samples)')
     end
   end

% ---------- A single observation stream ------
else
   subplot(211)
   plot(0:samples,yhat_data,'-',tidx,y,'+')
   grid
   set(gca,'Xlim',[0 samples])
   title('Observed output (x) and estimated output (line)')
   xlabel('Time (samples)')

   % --- Plot residuals ---
   subplot(212)
   plot(tidx,[y-yhat_data(tidx,:)])
   grid
   set(gca,'Xlim',[0 samples])
   title('Error between observed and estimated outputs')
   xlabel('Time (samples)')
end

% --- Plot state estimates ---
if strcmpi(method,'ekf'),
   statestd = 3*sqrt(mat2var(PS));
else
   statestd = 3*sqrt(smat2var(PS));
end
noplots = ceil(nx/Maxsubplots);
s = 1;    % State counter

for p=1:noplots,
   figure
   for i=1:Maxsubplots,
      subplot(Maxsubplots,1,i);
      plot(0:samples,x_data(:,s),'b-',0:samples,x_data(:,s)+statestd(:,s),'r:',...
                 0:samples,x_data(:,s)-statestd(:,s),'r:');
      set(gca,'Xlim',[0 samples]);
      ylabel(['State #' num2str(s)]);
      s = s+1;
      if s>nx, break; end
   end
   xlabel('Time (samples)');
end		 
