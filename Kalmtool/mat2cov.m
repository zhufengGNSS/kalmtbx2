function P = mat2cov(Pmat,idx)
% MAT2COV   Extract covariance matrix from vector of upper triangular elements.
%
%    P = mat2cov(Pvec) returns the (quadratic) covariance matrix when given a 
%    vector containing the upper triangular elements.
%    P = mat2cov(Pmat,k) extracts the k'th row from the matrix of vectors
%    Pmat. The vectors of upper triangular elements must be organized
%    row wise in Pmat.
 
% Written by: Magnus Norgaard, IMM/IAU, Technical University of Denmark
% LastEditDate: Apr. 15, 2000 
if nargin==1, 
  idx = 1;
end
states = round(-0.5 + sqrt(0.25 + 2*size(Pmat,2)));
P      = tril(reshape(1:states*states,states,states));
ii     = find(P);
P(ii)  = Pmat(idx,:);
P      = P + P' - diag(diag(P));
