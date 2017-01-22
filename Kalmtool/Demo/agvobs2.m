function y=agvobs(x,w)
% Prediction of observation
% 
% Check if variables should be initailized
if nargin==1
   return
end
y = x(1:2)+w;
