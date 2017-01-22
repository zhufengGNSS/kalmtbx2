function [M,N]=bodylin(x,u,vw,flag)
% Linearization of AGV model
persistent A0 F0 C0 G0 gamma H Mp delta I % Make variables static

% Check if variables should be initialized
if nargin==1
   delta = x(1);
   gamma = x(4);
   Mp    = x(2);
   H     = x(3);
   
   A0 = eye(3);
   F0 = [];
   C0 = zeros(3,5); C0(1,1) = 1; C0(2,2) = 1; C0(3,3) = 1;
   G0 = [];
   return
end

% Linearize state equation 
if flag==0,
   ak  = exp(-gamma*x(1))*delta;
   a12 = -ak*exp(gamma*x(1));
   a21 = ak*gamma*x(2)*x(2)*x(3);
   a22 = -ak*2*x(2)*x(3);
   a23 = -ak*x(2)*x(2);
   
   M = [1+0.5*a12*a21  a12+0.5*a12*a22  0.5*a12*a23;
        a21+0.5*a22*a21  1+a22+0.5*(a12*a21+a22*a22)  a23+0.5*a22*a23;
	0 0 1]; 
   N  = F0;

% Linearize output equation 
elseif flag==1,
   tmp = x(1)-H;
   M = [tmp/sqrt(Mp*Mp+tmp*tmp) 0 0];
   N = G0;
end
