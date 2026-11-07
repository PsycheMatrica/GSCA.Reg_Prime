function test_gsca_risk_maps
%TEST_GSCA_RISK_MAPS Independent algebra checks (no fitted-model simulation).
fit.mu = [1,-2,.5];
fit.scale = [2,3,4];
fit.W = eye(3);
fit.C = diag([.5,.7,.8]);
fit.B = [0,0,.4;0,0,-.2;0,0,0];
fit.converged = true;
C0 = eye(3);
B0 = fit.B~=0;
refMu = [-.2,.3,.4];
S = [4,1.2,.4;1.2,9,1;.4,1,16];
X = [1,2,3;-2,1,4;5,-1,2;0,3,-2;2,4,1;-1,0,5];
A = [ones(size(X,1),1),X];
maps = gsca_risk_maps(fit,C0,B0,refMu,S);
assert(maps.valid);
assert(isequal(maps.source,[1;2;3;1;2]));
assert(isequal(maps.target,[1;2;3;3;3]));
rawW = bsxfun(@rdivide,fit.W,fit.scale');
CV = bsxfun(@minus,X,fit.mu)*rawW;
for e=1:numel(maps.source)
    p = maps.source(e);
    q = maps.target(e);
    if strcmp(maps.kind{e},'C')
        directFull = X(:,q)-(fit.mu(q)+CV*fit.C(:,q)*fit.scale(q));
        directZero = directFull+CV(:,p)*fit.C(p,q)*fit.scale(q);
        multiplier = fit.C(p,q)*fit.scale(q);
        others = find(C0(:,q));
    else
        targetSD = sqrt(rawW(:,q)'*S*rawW(:,q));
        directFull = (CV(:,q)-CV*fit.B(:,q))/targetSD;
        directZero = directFull+CV(:,p)*fit.B(p,q)/targetSD;
        multiplier = fit.B(p,q)/targetSD;
        others = find(B0(:,q));
    end
    assert(max(abs(A*maps.full(:,e)-directFull))<1e-12);
    assert(max(abs(A*maps.zeroRemoved(:,e)-directZero))<1e-12);
    others = others(others~=p);
    projectionSlope = (maps.removed(2:end,e)-maps.full(2:end,e))/multiplier;
    assert(norm(rawW(:,others)'*S*projectionSlope)<1e-12);
    % The replacement residual has mean zero under reference moments.
    h = maps.removed(:,e)-maps.full(:,e);
    assert(abs([1,refMu]*h)<1e-12);
end
% Exact integrated risk must equal a direct average for empirical moments.
risk = gsca_exact_risk(maps,mean(X,1),cov(X,1));
assert(max(abs(risk.full-mean((A*maps.full).^2,1)'))<1e-12);
assert(max(abs(risk.deltaEPE-mean((A*maps.removed).^2-(A*maps.full).^2,1)'))<1e-12);
assert(max(abs(risk.zeroDeltaEPE-mean((A*maps.zeroRemoved).^2-(A*maps.full).^2,1)'))<1e-12);
% Common location shifts preserve all three residuals.
shift = [7,-5,2];
shiftFit = fit;
shiftFit.mu = fit.mu+shift;
shiftMaps = gsca_risk_maps(shiftFit,C0,B0,refMu+shift,S);
shiftA = [ones(size(X,1),1),bsxfun(@plus,X,shift)];
assert(max(abs(shiftA*shiftMaps.removed-A*maps.removed),[],'all')<1e-12);
% Changes of raw indicator units rescale C residuals, preserve B residuals.
units = [2,.5,3];
unitFit = fit;
unitFit.mu = fit.mu.*units;
unitFit.scale = fit.scale.*units;
unitMaps = gsca_risk_maps(unitFit,C0,B0,refMu.*units,S.*(units'*units));
unitA = [ones(size(X,1),1),bsxfun(@times,X,units)];
expected = A*maps.removed;
expected(:,1:3) = bsxfun(@times,expected(:,1:3),units);
assert(max(abs(unitA*unitMaps.removed-expected),[],'all')<1e-12);
% Joint sign changes preserve C risks and B squared risks.
D = diag([-1,1,-1]);
signFit = fit;
signFit.W = fit.W*D;
signFit.C = D*fit.C;
signFit.B = D*fit.B*D;
signMaps = gsca_risk_maps(signFit,C0,B0,refMu,S);
signRisk = gsca_exact_risk(signMaps,refMu,S);
risk = gsca_exact_risk(maps,refMu,S);
assert(max(abs(signRisk.deltaEPE-risk.deltaEPE))<1e-12);
% Singular parent projection is explicit, never silently pseudoinverted.
badFit = fit;
badFit.W(:,2) = 0;
badMaps = gsca_risk_maps(badFit,C0,B0,refMu,S);
assert(~badMaps.valid);
assert(~badMaps.edgeValid(4));
assert(strcmp(badMaps.edgeReason{4},'singular_parent_projection'));
% Nonconvergence propagates and cannot produce a seemingly valid risk.
badFit = fit;
badFit.converged = false;
badMaps = gsca_risk_maps(badFit,C0,B0,refMu,S);
assert(~badMaps.valid && all(isnan(badMaps.full(:))));
fprintf('PASS: affine residuals, projection, exact risk, units, signs, failure handling.\n');
end
