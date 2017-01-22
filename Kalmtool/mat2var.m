function vars = mat2var(Pmat)
% MAT2VAR   Extract variance estimates from matrix of covariance estimates.
%    The variance estimates (corresponding to the diagonal of a 
%    covariance matrix) are extracted from a matrix for which each
%    row contains the upper triangular elements of a covariance matrix.
%
%           varmat = mat2var(Pmat); 
% 
% Written by: Magnus Norgaard, IMM/IAU, Technical University of Denmark
% LastEditDate: Apr. 15, 2000 
states = round(-0.5 + sqrt(0.25 + 2*size(Pmat,2)));
idx    = cumsum ([1 states:-1:2]);
vars   = Pmat(:,idx);
