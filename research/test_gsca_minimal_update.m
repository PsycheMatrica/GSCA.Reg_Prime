function test_gsca_minimal_update
%TEST_GSCA_MINIMAL_UPDATE Algebra tests for the single-fit OOB helper only.
% No GSCA refitting, bootstrap, correction, or nested inference is introduced.
rootDirectory = fileparts(fileparts(mfilename('fullpath')));
addpath(rootDirectory);
test_supplied_fit;
test_singleton_ols;
fprintf('PASS: minimal single-fit auxiliary replacement and attached-indicator losses.\n');
end

function test_supplied_fit
i = (1:30)';
t = (31:38)';
ib = [sin(i),cos(i),sin(i*.7),cos(i*.4),sin(i*.3),cos(i*.8)];
oob = [sin(t),cos(t),sin(t*.7),cos(t*.4),sin(t*.3),cos(t*.8)];
location = [2,-1,.5,4,-3,1];
scale = [1,2,.8,1.5,3,.7];
ib = bsxfun(@plus,bsxfun(@times,ib,scale),location);
oob = bsxfun(@plus,bsxfun(@times,oob,scale),location);
W = [.8,0,0;.4,0,0;0,.7,0;0,.5,0;0,0,.9;0,0,.3];
C = [.7,.9,0,0,0,0;0,0,.8,.6,0,0;0,0,0,0,.85,.75];
B = [0,.25,.4;0,0,-.2;0,0,0];
C0 = C~=0;
B0 = B~=0;
Wcopy = W; Ccopy = C; Bcopy = B;
s = GSCA_Indicator_OOB(ib,oob,W,C,B,C0,B0);
assert(isequal(W,Wcopy) && isequal(C,Ccopy) && isequal(B,Bcopy));
assert(isequal(s.edgeIndex,find(B0)));
assert(isequal(s.CedgeIndex,find(C0)));
mu = mean(ib,1);
sd = std(ib,1,1);
Zi = bsxfun(@rdivide,bsxfun(@minus,ib,mu),sd);
Zo = bsxfun(@rdivide,bsxfun(@minus,oob,mu),sd);
Ci = Zi*W;
Co = Zo*W;
% Full B prediction uses only parent scores and fixed supplied coefficients.
fullScore = Co(:,1:2)*B(1:2,3);
directFull = bsxfun(@plus,bsxfun(@times,fullScore*C(3,5:6),sd(5:6)),mu(5:6));
assert(max(abs(s.fullPrediction_by_target{3}-directFull),[],'all')<1e-12);
fullError = oob(:,5:6)-directFull;
assert(abs(s.epeFull_by_target(3)-sum(fullError(:).^2)/numel(fullError))<1e-12);
% Compute the one-other-parent OLS replacement independently by covariances.
x = Ci(:,2); y = Ci(:,1);
slope = sum((x-mean(x)).*(y-mean(y)))/sum((x-mean(x)).^2);
intercept = mean(y)-slope*mean(x);
replacement = intercept+Co(:,2)*slope;
e = find(s.source==1 & s.target==3);
assert(max(abs(s.replacementScore{e}-replacement))<1e-12);
replacedScore = fullScore+(replacement-Co(:,1))*B(1,3);
directRemoved = bsxfun(@plus,bsxfun(@times,replacedScore*C(3,5:6),sd(5:6)),mu(5:6));
removedError = oob(:,5:6)-directRemoved;
assert(max(abs(s.removedPrediction_by_edge{e}-directRemoved),[],'all')<1e-12);
assert(abs(s.delOPE_b(e)-(mean(removedError(:).^2)-mean(fullError(:).^2)))<1e-12);
% A sole predictor has an intercept-only auxiliary regression.
single = find(s.source==1 & s.target==2);
assert(isempty(s.auxiliaryParents{single}));
assert(abs(s.auxiliaryCoefficients{single}-mean(Ci(:,1)))<1e-12);
% Changing OOB target indicators cannot alter their full/removed predictions.
changed = oob;
changed(:,5:6) = changed(:,5:6)+[100*ones(size(oob,1),1),-50*ones(size(oob,1),1)];
sc = GSCA_Indicator_OOB(ib,changed,W,C,B,C0,B0);
assert(max(abs(sc.fullPrediction_by_target{3}-s.fullPrediction_by_target{3}),[],'all')<1e-12);
assert(max(abs(sc.removedPrediction_by_edge{e}-s.removedPrediction_by_edge{e}),[],'all')<1e-12);
assert(isequal(sc.auxiliaryCoefficients{e},s.auxiliaryCoefficients{e}));
% C losses retain their legacy standardized-indicator reconstruction scale.
for k=1:numel(s.Csource)
    p = s.Csource(k); j = s.Ctarget(k);
    fullC = Co(:,p)*C(p,j);
    removedC = mean(Ci(:,p))*C(p,j);
    directDelta = mean((Zo(:,j)-removedC).^2)-mean((Zo(:,j)-fullC).^2);
    assert(abs(s.delOPE_c(k)-directDelta)<1e-12);
end
% Joint sign changes preserve all observed-indicator prediction differences.
D = diag([-1,1,-1]);
ss = GSCA_Indicator_OOB(ib,oob,W*D,D*C,D*B*D,C0,B0);
assert(max(abs(ss.delOPE_b-s.delOPE_b))<1e-12);
assert(max(abs(ss.delOPE_c-s.delOPE_c))<1e-12);
% Illegal use of a target indicator in a source score is rejected.
badW = W;
badW(5,1) = .1;
expect_error(@()GSCA_Indicator_OOB(ib,oob,badW,C,B,C0,B0), ...
    'GSCA_Indicator_OOB:UnsupportedModel');
end

function test_singleton_ols
i = (1:80)';
x1 = sin(i*sqrt(2));
x2 = .3*x1+cos(i*sqrt(3));
y = .2+.35*x1-.4*x2+.6*sin(i*sqrt(11));
ib = [x1,x2,y];
oob = [.2,.4,1;-.5,.3,-.2;.7,-.1,.6;-.2,-.4,.1];
mu = mean(ib,1);
sd = std(ib,1,1);
Zi = bsxfun(@rdivide,bsxfun(@minus,ib,mu),sd);
B0 = [0,0,1;0,0,1;0,0,0];
B = zeros(3);
B(1:2,3) = Zi(:,1:2)\Zi(:,3);
s = GSCA_Indicator_OOB(ib,oob,eye(3),eye(3),B,eye(3),B0);
beta = [ones(size(ib,1),1),ib(:,1:2)]\ib(:,3);
fullPrediction = [ones(size(oob,1),1),oob(:,1:2)]*beta;
assert(max(abs(s.fullPrediction_by_target{3}-fullPrediction))<1e-11);
% FWL: replacing a predictor by its in-bag projection reproduces reduced OLS.
red1 = [ones(size(ib,1),1),ib(:,2)]\ib(:,3);
red2 = [ones(size(ib,1),1),ib(:,1)]\ib(:,3);
assert(max(abs(s.removedPrediction_by_edge{1}-[ones(size(oob,1),1),oob(:,2)]*red1))<1e-11);
assert(max(abs(s.removedPrediction_by_edge{2}-[ones(size(oob,1),1),oob(:,1)]*red2))<1e-11);
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

