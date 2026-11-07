function test_gsca_indicator_maps
%TEST_GSCA_INDICATOR_MAPS Independent checks for fixed indicator prediction.
% Run with research/ on the MATLAB path.
test_multiblock;
test_singleton_ols;
fprintf('PASS: fixed-indicator maps, aggregation, OLS bridge, leakage, projection, units, signs.\n');
end

function test_multiblock
C0 = kron(eye(3),ones(1,2));
B0 = [0,0,1;0,0,1;0,0,0];
fit.mu = [1,-2,.5,2,-1,.3];
fit.scale = [2,3,4,1.5,.8,2.5];
fit.W = [ .8,0,0; .4,0,0; 0,.7,0; 0,.5,0; 0,0,.9; 0,0,.3];
fit.C = [ .7,.9,0,0,0,0; 0,0,.8,.6,0,0; 0,0,0,0,.85,.75];
fit.B = [0,0,.4;0,0,-.2;0,0,0];
fit.converged = true;
refMu = [-.2,.3,.4,-.1,.5,.7];
L = [1,0,0,0,0,0;.2,1,0,0,0,0;.3,.1,1,0,0,0; ...
     0,.2,.1,1,0,0;.1,0,.1,.2,1,0;0,.1,0,.1,.3,1];
S = L*L';
X = reshape(sin((1:60)*sqrt(2))+cos((1:60)*sqrt(3)),10,6);
A = [ones(size(X,1),1),X];
maps = gsca_indicator_maps(fit,C0,B0,refMu,S);
assert(maps.valid);
assert(isequal(maps.edgeSource,[1;2]));
assert(isequal(maps.edgeTarget,[3;3]));
assert(isequal(maps.group,[1;1;2;2]));
assert(isequal(maps.indicator,[5;6;5;6]));
assert(isequal(maps.indicatorCount,[2;2]));
assert(max(abs(sum(maps.aggregation,2)-1))<1e-15);
rawW = bsxfun(@rdivide,fit.W,fit.scale');
CV = bsxfun(@minus,X,fit.mu)*rawW;
for k=1:numel(maps.indicator)
    p = maps.source(k);
    q = maps.target(k);
    j = maps.indicator(k);
    parents = find(B0(:,q));
    directPrediction = fit.mu(j)+fit.scale(j)*fit.C(q,j)*(CV(:,parents)*fit.B(parents,q));
    directFull = X(:,j)-directPrediction;
    multiplier = fit.scale(j)*fit.C(q,j)*fit.B(p,q);
    assert(max(abs(A*maps.full(:,k)-directFull))<1e-12);
    assert(max(abs(A*maps.zeroRemoved(:,k)-(directFull+CV(:,p)*multiplier)))<1e-12);
    other = parents(parents~=p);
    projectionSlope = (maps.removed(2:end,k)-maps.full(2:end,k))/multiplier;
    assert(norm(rawW(:,other)'*S*projectionSlope)<1e-12);
    assert(abs([1,refMu]*(maps.removed(:,k)-maps.full(:,k)))<1e-12);
    % Every target-block prediction coefficient is exactly zero.
    unit = zeros(6,1);
    unit(j) = 1;
    assert(all(unit([5,6])-maps.full([6,7],k)==0));
    assert(all(unit([5,6])-maps.removed([6,7],k)==0));
end
% Change target indicators while holding predictor blocks and the fit fixed.
changedX = X;
changedX(:,5:6) = changedX(:,5:6)+reshape(1:20,10,2)*10;
changedA = [ones(size(X,1),1),changedX];
for name = {'full','removed','zeroRemoved'}
    R = maps.(name{1});
    originalPrediction = X(:,maps.indicator)-A*R;
    changedPrediction = changedX(:,maps.indicator)-changedA*R;
    assert(max(abs(originalPrediction-changedPrediction),[],'all')<1e-11);
end
% Exact integration equals a direct finite average, before and after grouping.
risk = gsca_exact_risk(maps,mean(X,1),cov(X,1));
directD = (A*maps.removed).^2-(A*maps.full).^2;
assert(max(abs(risk.deltaEPE-mean(directD,1)'))<1e-12);
assert(max(abs(maps.aggregation*risk.deltaEPE-mean(directD*maps.aggregation',1)'))<1e-12);
% Location changes cancel in residuals.
shift = [7,-5,2,3,-1,4];
shiftFit = fit;
shiftFit.mu = fit.mu+shift;
shiftMaps = gsca_indicator_maps(shiftFit,C0,B0,refMu+shift,S);
shiftA = [ones(size(X,1),1),bsxfun(@plus,X,shift)];
assert(max(abs(shiftA*shiftMaps.removed-A*maps.removed),[],'all')<1e-12);
% Raw-unit equivariance: each coordinate's residual scales with its outcome.
units = [2,.5,3,1.2,.7,4];
unitFit = fit;
unitFit.mu = fit.mu.*units;
unitFit.scale = fit.scale.*units;
unitMaps = gsca_indicator_maps(unitFit,C0,B0,refMu.*units,S.*(units'*units));
unitA = [ones(size(X,1),1),bsxfun(@times,X,units)];
expected = bsxfun(@times,A*maps.removed,units(maps.indicator));
assert(max(abs(unitA*unitMaps.removed-expected),[],'all')<1e-12);
% Joint component sign changes leave the observed-indicator maps invariant.
D = diag([-1,1,-1]);
signFit = fit;
signFit.W = fit.W*D;
signFit.C = D*fit.C;
signFit.B = D*fit.B*D;
signMaps = gsca_indicator_maps(signFit,C0,B0,refMu,S);
assert(max(abs(signMaps.full-maps.full),[],'all')<1e-12);
assert(max(abs(signMaps.removed-maps.removed),[],'all')<1e-12);
% Unsupported structures fail explicitly instead of leaking target outcomes.
crossC = C0;
crossC(1,5) = 1;
expect_error(@()gsca_indicator_maps(fit,crossC,B0,refMu,S),'gsca:indicatorMapsUnsupported');
leakFit = fit;
leakFit.W(5,1) = .1;
expect_error(@()gsca_indicator_maps(leakFit,C0,B0,refMu,S),'gsca:indicatorMapsUnsupported');
selfB = B0;
selfB(3,3) = 1;
expect_error(@()gsca_indicator_maps(fit,C0,selfB,refMu,S),'gsca:indicatorMapsUnsupported');
badFit = fit;
badFit.converged = false;
badMaps = gsca_indicator_maps(badFit,C0,B0,refMu,S);
assert(~badMaps.valid && all(isnan(badMaps.full(:))));
end

function test_singleton_ols
% An actual singleton GSCA fit must reproduce OLS predictions and removal.
i = (1:160)';
x1 = sin(i*sqrt(2));
x2 = .3*x1+cos(i*sqrt(3));
y = .2+.35*x1-.4*x2+.6*sin(i*sqrt(11));
X = [x1,x2,y];
W0 = eye(3);
C0 = eye(3);
B0 = [0,0,1;0,0,1;0,0,0];
fit = gsca_fit(X,W0,C0,B0,0,struct('Tol',1e-10,'MaxIter',1000));
assert(fit.converged,fit.reason);
maps = gsca_indicator_maps(fit,C0,B0,mean(X,1),cov(X,1));
assert(maps.valid);
A = [ones(size(X,1),1),X];
beta = [ones(size(X,1),1),X(:,1:2)]\y;
directFull = y-[ones(size(X,1),1),X(:,1:2)]*beta;
assert(max(abs(A*maps.full(:,1)-directFull))<1e-10);
assert(max(abs(A*maps.full(:,2)-directFull))<1e-10);
% With in-sample reference projection, FWL equals reduced OLS.
red1 = [ones(size(X,1),1),X(:,2)];
red2 = [ones(size(X,1),1),X(:,1)];
assert(max(abs(A*maps.removed(:,1)-(y-red1*(red1\y))))<1e-10);
assert(max(abs(A*maps.removed(:,2)-(y-red2*(red2\y))))<1e-10);
assert(isequal(maps.aggregation,eye(2)));
end

function expect_error(f,id)
try
    f();
catch err
    assert(strcmp(err.identifier,id),['Unexpected error: ',err.identifier]);
    return
end
error('test:MissingError','Expected error was not raised: %s',id);
end

