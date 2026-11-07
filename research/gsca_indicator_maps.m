function maps = gsca_indicator_maps(fit,C0,B0,refMu,refSigma)
%GSCA_INDICATOR_MAPS Fixed observed-indicator risks for structural GSCA paths.
%   maps = gsca_indicator_maps(fit,C0,B0,refMu,refSigma)
%
%   For each free B(p,q), predict every observed indicator j in block q:
%     yhat_j = mu_j + scale_j*C(q,j)*sum_{s in parents(q)} CV_s*B(s,q).
%   CV_s uses ONLY the fitted parent block's observed indicators. The target
%   block's own indicators are never inputs to this prediction. Outcomes
%   are the fixed observed X(:,j), not re-estimated endogenous scores.
%
%   [ones(N,1),X]*maps.full/removed/zeroRemoved returns affine residuals,
%   one column per edge-indicator coordinate. Projection removal replaces
%   CV_p by its affine projection onto the OTHER specified parent scores.
%   It uses refMu/refSigma conditional on the fitted parent weights.
%   Population moments give the oracle projection; full-X moments give a
%   plug-in approximation. No EPE denominator or target-score normalization
%   is used. Raw indicator-unit changes therefore rescale the losses.
%
%   Groups are B-only edges in find(B0) order. Within each group, coordinates
%   follow find(C0(q,:)). maps.aggregation is E-by-K, with 1/J_q entries:
%     edgeRisk = maps.aggregation * coordinateRisk
%   Here E is the number of B edges and K the number of coordinates.
%
%   Supported measurement structure: simple reflective, disjoint blocks;
%   fitted W support within the corresponding C0 block; no self paths.
%   Unsupported crossloadings, overlapping score support, and target-block
%   leakage throw errors. Numerical failures return valid=false with NaNs.
required = {'mu','scale','W','C','B','converged'};
if ~isstruct(fit) || ~isscalar(fit) || ~all(isfield(fit,required))
    error('gsca:indicatorMapsInput','fit lacks required fields.');
end
if ~isnumeric(fit.W) || ~isreal(fit.W) || ndims(fit.W)>2
    error('gsca:indicatorMapsInput','fit.W must be a real numeric matrix.');
end
[J,P] = size(fit.W);
if J<1 || P<1 || ~isequal(size(C0),[P,J]) || ...
        ~isequal(size(B0),[P,P]) || ~isequal(size(fit.C),[P,J]) || ...
        ~isequal(size(fit.B),[P,P])
    error('gsca:indicatorMapsInput','Incompatible model dimensions.');
end
if (~isnumeric(C0) && ~islogical(C0)) || ...
        (~isnumeric(B0) && ~islogical(B0)) || ...
        ~isreal(C0) || ~isreal(B0) || ...
        any(~isfinite(C0(:))) || any(~isfinite(B0(:)))
    error('gsca:indicatorMapsInput','Masks must be finite and real.');
end
Cmask = C0~=0;
Bmask = B0~=0;
if any(sum(Cmask,1)>1)
    error('gsca:indicatorMapsUnsupported','Crossloadings or overlapping measurement blocks are unsupported.');
end
if any(sum(Cmask,2)==0)
    error('gsca:indicatorMapsUnsupported','Each component requires a specified indicator block.');
end
if any(diag(Bmask))
    error('gsca:indicatorMapsUnsupported','Self paths would use target indicators as prediction inputs.');
end
if ~isnumeric(refMu) || ~isreal(refMu) || numel(refMu)~=J || ...
        ~isnumeric(refSigma) || ~isreal(refSigma) || ...
        ~isequal(size(refSigma),[J,J]) || ...
        any(~isfinite(refMu(:))) || any(~isfinite(refSigma(:)))
    error('gsca:indicatorMapsInput','Reference moments must be finite and real.');
end
refMu = reshape(refMu,1,J);
symTol = 100*eps*max(norm(refSigma,'fro'),realmin);
if norm(refSigma-refSigma','fro')>symTol
    error('gsca:indicatorMapsInput','refSigma must be symmetric.');
end
refSigma = (refSigma+refSigma')/2;
if min(eig(refSigma)) < -100*eps*max(norm(refSigma,2),realmin)
    error('gsca:indicatorMapsInput','refSigma must be positive semidefinite.');
end

[edgeSource,edgeTarget] = find(Bmask);
E = numel(edgeSource);
indicatorCount = sum(Cmask(edgeTarget,:),2);
K = sum(indicatorCount);
maps.full = nan(J+1,K);
maps.removed = nan(J+1,K);
maps.zeroRemoved = nan(J+1,K);
maps.kind = repmat({'B'},K,1);
maps.group = zeros(K,1);
maps.indicator = zeros(K,1);
maps.source = zeros(K,1);
maps.target = zeros(K,1);
maps.index = zeros(K,1);
maps.edgeKind = repmat({'B'},E,1);
maps.edgeSource = edgeSource;
maps.edgeTarget = edgeTarget;
maps.edgeIndex = find(Bmask);
maps.indicatorCount = indicatorCount;
maps.nCoordinates = indicatorCount;
maps.aggregation = zeros(E,K);
maps.coordinateValid = false(K,1);
maps.coordinateReason = repmat({''},K,1);
maps.edgeValid = false(E,1);
maps.edgeReason = repmat({''},E,1);
maps.valid = false;
maps.reason = '';
start = 0;
for e=1:E
    columns = start+(1:indicatorCount(e));
    indicators = find(Cmask(edgeTarget(e),:));
    maps.group(columns) = e;
    maps.indicator(columns) = indicators(:);
    maps.source(columns) = edgeSource(e);
    maps.target(columns) = edgeTarget(e);
    maps.index(columns) = maps.edgeIndex(e);
    maps.aggregation(e,columns) = 1/indicatorCount(e);
    start = start+indicatorCount(e);
end
if ~isscalar(fit.converged) || ~isfinite(fit.converged) || ~fit.converged
    maps = fail_all(maps,'fit_not_converged');
    return
end
if ~isnumeric(fit.mu) || ~isnumeric(fit.scale) || ...
        ~isreal(fit.mu) || ~isreal(fit.scale) || ...
        numel(fit.mu)~=J || numel(fit.scale)~=J || ...
        any(~isfinite(fit.mu(:))) || any(~isfinite(fit.scale(:))) || ...
        any(fit.scale(:)<=0) || any(~isfinite(fit.W(:))) || ...
        ~isnumeric(fit.C) || ~isreal(fit.C) || any(~isfinite(fit.C(:))) || ...
        ~isnumeric(fit.B) || ~isreal(fit.B) || any(~isfinite(fit.B(:)))
    maps = fail_all(maps,'nonfinite_or_invalid_fit');
    return
end
% Check exact support, rather than silently zeroing illegal fitted entries.
if any(fit.W(~Cmask')~=0)
    error('gsca:indicatorMapsUnsupported','Fitted weights use indicators outside their own simple measurement blocks.');
end
if any(fit.C(~Cmask)~=0) || any(fit.B(~Bmask)~=0)
    error('gsca:indicatorMapsUnsupported','Fitted coefficients violate the supplied model masks.');
end
if any(sum(abs(fit.W),1)==0)
    maps = fail_all(maps,'zero_component_weights');
    return
end
mu = reshape(fit.mu,1,J);
sigma = reshape(fit.scale,1,J);
rawW = bsxfun(@rdivide,fit.W,sigma');
for e=1:E
    p = edgeSource(e);
    q = edgeTarget(e);
    columns = find(maps.group==e);
    targetIndicators = maps.indicator(columns);
    parents = find(Bmask(:,q));
    if any(any(rawW(targetIndicators,parents)~=0))
        error('gsca:indicatorMapsUnsupported','Target indicators occur in a parent score used for prediction.');
    end
    other = parents(parents~=p);
    residualWeight = rawW(:,p);
    if ~isempty(other)
        U = rawW(:,other);
        G = U'*refSigma*U;
        G = (G+G')/2;
        if any(~isfinite(G(:))) || rcond(G)<1e-10
            maps.edgeReason{e} = 'singular_parent_projection';
            maps.coordinateReason(columns) = {maps.edgeReason{e}};
            continue
        end
        residualWeight = residualWeight-U*(G\(U'*refSigma*residualWeight));
    end
    parentPrediction = rawW(:,parents)*fit.B(parents,q);
    for k=reshape(columns,1,[])
        j = maps.indicator(k);
        loadingScale = fit.C(q,j)*sigma(j);
        predictionSlope = parentPrediction*loadingScale;
        unit = zeros(J,1);
        unit(j) = 1;
        fullSlope = unit-predictionSlope;
        full = [-mu*fullSlope;fullSlope];
        multiplier = fit.B(p,q)*loadingScale;
        removedTerm = multiplier*[-refMu*residualWeight;residualWeight];
        zeroTerm = multiplier*[-mu*rawW(:,p);rawW(:,p)];
        if any(~isfinite([full;removedTerm;zeroTerm]))
            maps.coordinateReason{k} = 'nonfinite_residual_map';
            continue
        end
        maps.full(:,k) = full;
        maps.removed(:,k) = full+removedTerm;
        maps.zeroRemoved(:,k) = full+zeroTerm;
        maps.coordinateValid(k) = true;
    end
    maps.edgeValid(e) = all(maps.coordinateValid(columns));
    if ~maps.edgeValid(e) && isempty(maps.edgeReason{e})
        bad = columns(find(~maps.coordinateValid(columns),1));
        maps.edgeReason{e} = maps.coordinateReason{bad};
    end
end
maps.valid = all(maps.edgeValid);
if ~maps.valid
    bad = find(~maps.edgeValid,1);
    maps.reason = sprintf('edge_%d:%s',bad,maps.edgeReason{bad});
end
end

function maps = fail_all(maps,reason)
maps.reason = reason;
maps.edgeReason(:) = {reason};
maps.coordinateReason(:) = {reason};
end

