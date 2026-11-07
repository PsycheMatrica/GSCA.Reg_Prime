function risk = gsca_exact_risk(maps,popMu,popSigma)
%GSCA_EXACT_RISK Integrate affine squared residuals over known X moments.
%   This computes conditional risk for ONE fitted model, not the finite-N
%   training-averaged truth. Average independent training fits for that.
%   It remains exact for non-Gaussian X with these finite second moments.
%   The maps include all centering and outcome-scaling terms.
J = size(maps.full,1)-1;
if numel(popMu)~=J || ~isequal(size(popSigma),[J,J]) || ...
        any(~isfinite(popMu(:))) || any(~isfinite(popSigma(:))) || ...
        ~isreal(popMu) || ~isreal(popSigma)
    error('gsca:exactRiskInput','Population moments have invalid dimensions or values.');
end
popMu = reshape(popMu,1,J);
symTol = 100*eps*max(norm(popSigma,'fro'),realmin);
if norm(popSigma-popSigma','fro')>symTol
    error('gsca:exactRiskInput','Population covariance must be symmetric.');
end
popSigma = (popSigma+popSigma')/2;
if min(eig(popSigma)) < -100*eps*max(norm(popSigma,2),realmin)
    error('gsca:exactRiskInput','Population covariance must be positive semidefinite.');
end
risk.full = affineLoss(maps.full,popMu,popSigma);
risk.removed = affineLoss(maps.removed,popMu,popSigma);
risk.zeroRemoved = affineLoss(maps.zeroRemoved,popMu,popSigma);
risk.deltaEPE = risk.removed-risk.full;
risk.zeroDeltaEPE = risk.zeroRemoved-risk.full;
end

function value = affineLoss(A,mu,Sigma)
slope = A(2:end,:);
residualMean = A(1,:)+mu*slope;
value = (sum(slope.*(Sigma*slope),1)+residualMean.^2)';
end
