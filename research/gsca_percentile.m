function [p,lower]=gsca_percentile(draws,alpha)
%GSCA_PERCENTILE One-sided empirical lower bound and its inverse p-value.
% Exact numerical convention: reject iff p<alpha, equivalent to lower>0.
% Correct the candidate ceil(alpha*B) rank using k/B comparisons to avoid
% floating multiplication putting an exact grid probability in the next bin.
validateattributes(draws,{'numeric'},{'2d','real','finite','nonempty'});
validateattributes(alpha,{'numeric'},{'scalar','real','>',0,'<',1});
B=size(draws,1); rank=max(1,ceil(alpha*B));
while rank>1 && (rank-1)/B>=alpha, rank=rank-1; end
while rank<B && rank/B<alpha, rank=rank+1; end
p=mean(draws<=0,1);
s=sort(draws,1); lower=s(rank,:);
end
