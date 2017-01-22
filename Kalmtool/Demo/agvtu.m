function xout=agvtu(x,u,v)
% A priori update of the states in the AGV model
% 
persistent k1;                         % Make variables static

% Check if variables should be initailized
if nargin==1,
   kenc = 800/(2*pi);                  % Encoder gain
   N    = 36;                          % Gear factor
   k1 = 1/(kenc*N);                    % a useful parameter
   return
end

xout   = zeros(6,1);
u      = u + v(1:2);                   % Add process noise

s = 0.5*k1*(x(4)*u(1)+x(5)*u(2));
t = 0.5*k1*(x(4)*u(1)-x(5)*u(2))/x(6);

xout(1)   = x(1) + s*cos(x(3)+t);
xout(2)   = x(2) + s*sin(x(3)+t);
xout(3)   = x(3) + 2*t;
xout(4:6) = x(4:6);                    % Maintain parameters
