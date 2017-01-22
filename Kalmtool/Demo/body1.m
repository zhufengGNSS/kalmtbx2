function x=body1(x,u,v)
% 4th order Runge-Kutta solution of differential equation in 'bodydot'
persistent delta
if nargin==1,
  delta = x(1);
  bodydot(x);
else
  k1 = delta*bodydot(x,u,v);
  k2 = delta*bodydot(x+0.5*k1,u,v);
  k3 = delta*bodydot(x+0.5*k2,u,v);
  k4 = delta*bodydot(x+k3,u,v);

  x  = x + (k1+2*k2+2*k3+k4)/6;
end
