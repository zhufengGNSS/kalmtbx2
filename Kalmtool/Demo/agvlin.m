function [M,N]=agvlin(x,u,vw,flag)
% Linearization of AGV model
persistent kappa xout A0 F0 C0 G0  % Make variables static

% Check if variables should be initailized
if nargin==1
   kenc  = 800/(2*pi);                   % Encoder gain
   Ng=36;                                % Gear factor
   kappa = 0.5/(Ng*kenc);
   xout = zeros(6,1);
   
   A0 = diag([1 1 1 1 1 1]);
   F0 = zeros(6,2);
   C0 = zeros(3,6); C0(1,1) = 1; C0(2,2) = 1; C0(3,3) = 1;
   G0 = [];
   return
end

% Linearize state equation 
if flag==0,
   x1 = x(1); x2 = x(2); x3 = x(3); x4 = x(4); x5 = x(5); x6=x(6);
   u      = u + vw(1:2);
   s   = kappa*(x4*u(1) + x5*u(2));
   t   = kappa*(x4*u(1) - x5*u(2))/x6;
   phi = x3+t;
   a   = kappa*u(1);
   b   = kappa*u(2);
   alpha = kappa*x4;
   beta  = kappa*x5;
   sinp  = sin(phi);
   cosp  = cos(phi);
   ssinp = s*sinp;
   scosp = s*cosp;
   ax6   = a/x6;
   bx6   = b/x6;
   stx6  = s*t/x6;
   alpx6 = alpha/x6;
   betx6 = beta/x6;
   M=A0; N=F0;

   M(1,3) = -ssinp;
   M(1,4) = a*cosp - ax6*ssinp;
   M(1,5) = b*cosp + bx6*ssinp;
   M(1,6) = stx6*sinp;
   M(2,3) = scosp;
   M(2,4) = a*sinp + ax6*scosp;
   M(2,5) = b*sinp - bx6*scosp;
   M(2,6) = -stx6*cosp;
   M(3,4) = 2*a/x6;
   M(3,5) = -2*b/x6;
   M(3,6) = -2*t/x6;

   N(1,1) = alpha*cosp - alpx6*ssinp;
   N(1,2) = beta*cosp + betx6*ssinp;
   N(2,1) = alpha*sinp + alpx6*scosp;
   N(2,2) = beta*sinp - betx6*scosp;
   N(3,1) = 2*alpx6;
   N(3,2) = -2*betx6;


% Linearize output equation 
elseif flag==1,
   M=C0; N=G0;
elseif flag==2,
   M=C0(1:2,:); N=eye(2);
end
