function y=body2(x,w)
% Falling body example. Output equation.
persistent M H
if nargin==1,
  M = x(2);
  H = x(3);
else
  tmp = x(1)-H;
  y = sqrt(M*M+tmp*tmp)+w;
end
