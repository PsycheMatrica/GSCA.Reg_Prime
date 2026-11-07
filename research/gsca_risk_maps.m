function maps = gsca_risk_maps(fit,C0,B0,refMu,refSigma)
%GSCA_RISK_MAPS Affine full/removal residuals for exploratory GSCA inference.
%   [ones(N,1), X] * maps.full gives one full residual per specified edge.
%   maps.removed uses linear-projection replacement of the source score.
%   maps.zeroRemoved uses coefficient deletion on the same outcome scale.
%
%   Columns follow find(C0), then find(B0), in MATLAB column-major order.
%   C-edge losses are in raw observed-indicator units. B-edge losses concern
%   the FITTED endogenous component, with variance one under refSigma.
%   The latter remains an algorithm-dependent outcome, not a fixed latent
%   population component. Paired RESIDUAL variances are needed across fits.
%
%   Projection includes an intercept and only the other specified parents
%   of the same equation. Population refMu/refSigma give the population
%   projection conditional on the fitted weights. Full-data empirical
%   moments instead define a plug-in nuisance estimate, not an oracle rule.
%   Structural normalization fixes an outcome unit; no EPE denominator or
%   ratio is used. Component signs must be aligned before calling.
%
%   A failed fit or singular projection returns invalid/NaN columns and a
%   reason. Malformed inputs throw an error. No pseudoinverse is used.
required = {'mu','scale','W','C','B','converged'};
if ~isstruct(fit) || ~all(isfield(fit,required))
    error('gsca:riskMapsInput','fit lacks required fields.');
end
[J,P] = size(fit.W);
if ~isequal(size(C0),[P,J]) || ~isequal(size(B0),[P,P]) || ...
        ~isequal(size(fit.C),[P,J]) || ~isequal(size(fit.B),[P,P])
    error('gsca:riskMapsInput','Incompatible model dimensions.');
end
if ~isnumeric(refMu) || ~isreal(refMu) || numel(refMu)~=J || ...
        ~isnumeric(refSigma) || ~isreal(refSigma) || ...
        ~isequal(size(refSigma),[J,J]) || ...
        any(~isfinite(refMu(:))) || any(~isfinite(refSigma(:)))
    error('gsca:riskMapsInput','Reference moments must be finite and real.');
end
if any(~isfinite(C0(:))) || any(~isfinite(B0(:)))
    error('gsca:riskMapsInput','Masks must be finite.');
end
refMu = reshape(refMu,1,J);
symTol = 100*eps*max(norm(refSigma,'fro'),realmin);
if norm(refSigma-refSigma','fro')>symTol
    error('gsca:riskMapsInput','refSigma must be symmetric.');
end
refSigma = (refSigma+refSigma')/2;
if min(eig(refSigma)) < -100*eps*max(norm(refSigma,2),realmin)
    error('gsca:riskMapsInput','refSigma must be positive semidefinite.');
end
Cmask = C0~=0;
Bmask = B0~=0;
[cs,ct] = find(Cmask);
[bs,bt] = find(Bmask);
nC = numel(cs);
E = nC+numel(bs);
maps.full = nan(J+1,E);
maps.removed = nan(J+1,E);
maps.zeroRemoved = nan(J+1,E);
maps.kind = [repmat({'C'},nC,1);repmat({'B'},numel(bs),1)];
maps.source = [cs;bs];
maps.target = [ct;bt];
maps.index = [find(Cmask);find(Bmask)];
maps.edgeValid = false(E,1);
maps.edgeReason = repmat({''},E,1);
maps.targetScale = nan(E,1);
maps.valid = false;
maps.reason = '';
if ~isscalar(fit.converged) || ~fit.converged
    maps.reason = 'fit_not_converged';
    maps.edgeReason(:) = {maps.reason};
    return
end
if numel(fit.mu)~=J || numel(fit.scale)~=J || ...
        any(~isfinite(fit.mu(:))) || any(~isfinite(fit.scale(:))) || ...
        any(fit.scale(:)<=0) || any(~isfinite(fit.W(:))) || ...
        any(~isfinite(fit.C(:))) || any(~isfinite(fit.B(:))) || ...
        ~isreal(fit.W) || ~isreal(fit.C) || ~isreal(fit.B)
    maps.reason = 'nonfinite_or_invalid_fit';
    maps.edgeReason(:) = {maps.reason};
    return
end
mu = reshape(fit.mu,1,J);
sigma = reshape(fit.scale,1,J);
rawW = bsxfun(@rdivide,fit.W,sigma');
covNorm = norm(refSigma,2);
for e=1:E
    p = maps.source(e);
    q = maps.target(e);
    if e<=nC
        parents = find(Cmask(:,q));
        unit = zeros(J,1);
        unit(q) = 1;
        fullSlope = unit-rawW*fit.C(:,q)*sigma(q);
        multiplier = fit.C(p,q)*sigma(q);
        maps.targetScale(e) = 1;
    else
        parents = find(Bmask(:,q));
        targetVar = rawW(:,q)'*refSigma*rawW(:,q);
        varianceScale = (norm(rawW(:,q))^2)*covNorm;
        if ~isfinite(targetVar) || targetVar<=eps*max(varianceScale,realmin)
            maps.edgeReason{e} = 'zero_reference_target_variance';
            continue
        end
        targetSD = sqrt(targetVar);
        fullSlope = (rawW(:,q)-rawW*fit.B(:,q))/targetSD;
        multiplier = fit.B(p,q)/targetSD;
        maps.targetScale(e) = targetSD;
    end
    other = parents(parents~=p);
    residualWeight = rawW(:,p);
    if ~isempty(other)
        U = rawW(:,other);
        G = U'*refSigma*U;
        G = (G+G')/2;
        if any(~isfinite(G(:))) || rcond(G)<1e-10
            maps.edgeReason{e} = 'singular_parent_projection';
            continue
        end
        residualWeight = residualWeight-U*(G\(U'*refSigma*residualWeight));
    end
    full = [-mu*fullSlope;fullSlope];
    removedTerm = multiplier*[-refMu*residualWeight;residualWeight];
    zeroTerm = multiplier*[-mu*rawW(:,p);rawW(:,p)];
    if any(~isfinite([full;removedTerm;zeroTerm]))
        maps.edgeReason{e} = 'nonfinite_residual_map';
        continue
    end
    maps.full(:,e) = full;
    maps.removed(:,e) = full+removedTerm;
    maps.zeroRemoved(:,e) = full+zeroTerm;
    maps.edgeValid(e) = true;
end
maps.valid = all(maps.edgeValid);
if ~maps.valid
    bad = find(~maps.edgeValid,1);
    maps.reason = sprintf('edge_%d:%s',bad,maps.edgeReason{bad});
end
end
