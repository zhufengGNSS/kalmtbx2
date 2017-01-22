function xdot=bodydot(x,u,v)
% Falling body example. State equation.
persistent gamma
if nargin==1,
  gamma = x(4);
else
  xdot = [-x(2)+v(1);
          -exp(-gamma*x(1))*x(2)*x(2)*x(3)+v(2);
          v(3)];
end
