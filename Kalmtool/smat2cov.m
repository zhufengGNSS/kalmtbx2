function P = smat2cov(Smat,idx)
% SMAT2COV   Make covariance matrix from vector of Cholesky factor elements.
%    P = smat2cov(Svec) returns the (quadratic) covariance matrix when given a 
%    vector containing the (upper triangular) Cholesky factor elements.
%    P = smat2cov(Smat,k) extracts the k'th row from the matrix of vectors
%    Smat. The vectors of Cholesky factor elements must be organized
%    row wise in Smat.
% 
% Written by: Magnus Norgaard, IMM/IAU, Technical University of Denmark
% LastEditDate: Apr. 15, 2000 
if nargin==1, 
  idx = 1;
end
states = round(-0.5 + sqrt(0.25 + 2*size(Smat,2)));
P      = triu(reshape(1:states*states,states,states));
[I,J]  = find(P');
ii     = sub2ind([states states],J,I);
P(ii)  = Smat(idx,:);
P      = P*P';
