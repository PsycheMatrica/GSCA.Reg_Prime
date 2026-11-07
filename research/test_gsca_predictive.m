function test_gsca_predictive()
% Meaningful regression checks for the research estimator, no toolboxes.
rng(5723,'twister');
N=100; X=randn(N,2); y=.2*X(:,1)+.5*X(:,2)+randn(N,1);
Data=[X,y]; W0=eye(3); C0=zeros(3); B0=zeros(3);
B0(1:2,3)=1;
o=struct('Boot',199,'Seed',731,'ReferenceMean',zeros(1,3),...
    'ReferenceCov',[1,0,.2;0,1,.5;.2,.5,1.29]);
r=gsca_predictive(Data,W0,C0,B0,o);
assert(r.valid,r.reason); assert(r.minOOBCount>=2);
assert(r.identityError<1e-10,'Paired residual/bagging identity failed.');
assert(all(abs(r.deltaCorrected-r.deltaRaw-r.correction)<1e-12));
assert(isequal(r.pRaw,mean(r.draws.raw<=0,1)));
s=sort(r.draws.raw,1); expected=s(ceil(.05*o.Boot),:);
assert(isequal(r.lowerRaw,expected));
assert(isequal(r.lowerRaw>0,r.pRaw<.05));
% An exact grid probability must not cross bins through roundoff.
[pgrid,lgrid]=gsca_percentile([zeros(7,1);ones(93,1)],.07);
assert(pgrid==.07 && lgrid==0 && ~(pgrid<.07));
[pgrid,lgrid]=gsca_percentile([zeros(7,1);ones(93,1)],.07000000000000002);
assert(lgrid==1 && pgrid<.07000000000000002);
threw=false;
try, gsca_percentile([0;Inf],.05); catch, threw=true; end
assert(threw,'Nonfinite draws must not produce a p-value.');
% Seed controls every resample and preserves caller RNG.
st=rng; r2=gsca_predictive(Data,W0,C0,B0,o); after=rng;
assert(isequal(st,after)); assert(isequal(r.draws.raw,r2.draws.raw));
% Confirm the moving-target covariance term cannot be omitted.
t=[1;3;2;6]; f=[.2;1;2;3]; h=[.1;2;1;2]; g=f-h;
ef=t-f; eg=t-g;
cmat=cov([t,h]);
residualCorrection=var(ef,0)-var(eg,0);
predictionOnly=var(f,0)-var(g,0);
assert(abs(residualCorrection-(predictionOnly-2*cmat(1,2)))<1e-12);
assert(abs(residualCorrection-predictionOnly)>.1);
assert(abs(mean(eg.^2-ef.^2)+residualCorrection-...
    (mean(eg)^2-mean(ef)^2+residualCorrection/numel(t)))<1e-12);
% Exercise whole-estimator resampling; this checks execution, not calibration.
o.Boot=99; o.Outer=9; o.InnerBoot=99;
nested=gsca_predictive(Data,W0,C0,B0,o);
assert(nested.valid && nested.outerValid,'Nested estimator failed.');
assert(isequal(nested.pOuterCorrected,mean(nested.outer.corrected<=0,1)));
s=sort(nested.outer.corrected,1);
assert(isequal(nested.lowerOuterCorrected,s(1,:)));
fprintf('test_gsca_predictive: passed (including nested smoke; no calibration claim).\n');
end
