function vars = smat2var(Smat)
% SMAT2VAR   Calculate variance estimate for each state.
%    SMAT2VAR(Smat) returns a matrix where each column is the variance
%    of a state estimate. 'Smat' is a matrix where each row contains elements 
%    of (the upper triangular part of) the Cholesky factor of a covariance 
%    matrix.
%
%           varmat = smat2var(Smat); 
% 
% Written by: Magnus Norgaard, IMM/IAU, Technical University of Denmark
% LastEditDate: Apr. 15, 2000 
states = round(-0.5 + sqrt(0.25 + 2*size(Smat,2))); % Number of states
idx    = cumsum ([1 states:-1:2]);     % Index to first element in each row
idx2   = idx+[states-1:-1:0];          % Index to last element in each row
vars   = zeros(size(Smat,1),states);   % Allocate space for variance matrix
for k=1:states,
   vars(:,k) = sum(Smat(:,idx(k):idx2(k)).*Smat(:,idx(k):idx2(k)),2);
end
