function y=agvobs(x,w)
% Prediction of observation
% 
% Check if variables should be initailized
persistent C         % Make variables static
if nargin==1
   return
end
y = x(1:3)+w;
