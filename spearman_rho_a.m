function r = spearman_rho_a(X,Y)

% function r = spearman_rho_a(X,Y)
% 
% Function to compute Spearman's rho a, a version of Spearman rho that breaks
% ties randomly (akin to Kendall's tau a). Overall this leads to accurate
% inference for comparing representational similarities and thus obviates
% the need for Kendall's tau a. Without tied ranks, the result is identical
% to Spearman's rho.
% Adapted from compare_rho_a in https://github.com/rsagroup/rsatoolbox
%
% r = spearman_rho_a(X) returns a P-by-P matrix containing the pairwise
% Spearman rank correlation coefficient between each pair of columns in the
% N-by-P matrix X.
%
% For vectors,
% r = spearman_rho_a(X,Y) returns a single rank correlation coefficient.
% For matrices,
% r = spearman_rho_a(X,Y) returns a P1-by-P2 matrix containing the pairwise
% Spearman rank correlation coefficient between each pair of columns in the
% N-by-P1 and N-by-P2 matrices X and Y.
%
% Martin Hebart 2022/10/11
%
% See also corr, tiedrank

[n,m] = size(X);
assert(all(isreal(X(:))),'Input contains complex numbers. Cannot compute rank correlations.')
Xranked = tiedrank(X);
if nargin == 1
    Xranked = Xranked-mean(Xranked);
    r = 12*(Xranked'*Xranked)/(n^3-n);
    r(r>1) = 1; r(r<-1) = -1; % make sure result does not exceed 1 or -1 (this version deals better with NaN)
    r = r - diag(diag(r)) + eye(m); % make sure diagonal is 1 (unless all is NaN)
    return
end

% below code only runs if we are comparing two matrices X and Y
assert(size(Y,1)==n,'Both inputs need to be of the same length in the first dimension.')
assert(all(isreal(Y(:))),'Input contains complex numbers. Cannot compute rank correlations.')
Yranked = tiedrank(Y);
Yranked = Yranked-mean(Yranked);
r = (12*Xranked'*Yranked)/(n^3-n);
r(r>1) = 1; r(r<-1) = -1; % make sure result does not exceed 1 or -1 (this version deals better with NaN)
