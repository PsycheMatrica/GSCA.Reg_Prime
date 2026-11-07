function test_gsca_indicator_predictive()
% Independent arithmetic tests for fixed observed-indicator OOB inference.
previous=rng; cleanup=onCleanup(@()rng(previous)); %#ok<NASGU>
rng(4921,'twister');
N=100;
x1=2+1.2*randn(N,1);
x2=-3+.4*(x1-2)+.8*randn(N,1);
y=4+.35*x1-.2*x2+1.7*randn(N,1);
Data=[x1,x2,y];
W0=eye(3); C0=eye(3); B0=zeros(3); B0(1:2,3)=1;
o=struct('Boot',119,'Seed',7132,'Alpha',.05);
before=rng;
r=gsca_indicator_predictive(Data,W0,C0,B0,o);
assert(r.valid,r.reason);
assert(isequal(before,rng),'Seeded call changed caller RNG.');
assert(size(r.draws.raw,2)==2 && numel(r.deltaRaw)==2,'B-edge indexing changed.');

% Reconstruct every bootstrap with raw-data OLS and ordinary auxiliary OLS.
% No GSCA fit/map helper is used in this independent bridge calculation.
rng(o.Seed,'twister');
ef=nan(N,o.Boot,2); eg=ef; ez=ef;
draws=nan(o.Boot,2); zeroDraws=draws;
for b=1:o.Boot
    ib=randi(N,N,1); used=false(N,1); used(ib)=true; ix=find(~used);
    beta=[ones(N,1),Data(ib,1:2)]\Data(ib,3);
    full=Data(ix,3)-[ones(numel(ix),1),Data(ix,1:2)]*beta;
    for p=1:2
        other=3-p;
        aux=[ones(N,1),Data(:,other)]\Data(:,p);
        residual=Data(ix,p)-[ones(numel(ix),1),Data(ix,other)]*aux;
        removed=full+beta(p+1)*residual;
        deleted=full+beta(p+1)*(Data(ix,p)-mean(Data(ib,p)));
        ef(ix,b,p)=full; eg(ix,b,p)=removed; ez(ix,b,p)=deleted;
        draws(b,p)=mean(removed.^2-full.^2);
        zeroDraws(b,p)=mean(deleted.^2-full.^2);
    end
end
rowRaw=nan(N,2); rowCorrection=rowRaw; counts=zeros(N,1);
predictionCorrection=rowRaw; bagged=rowRaw; finiteB=rowRaw;
for i=1:N
    keep=isfinite(ef(i,:,1)); counts(i)=sum(keep);
    for p=1:2
        a=reshape(ef(i,keep,p),[],1); b=reshape(eg(i,keep,p),[],1);
        rowRaw(i,p)=mean(b.^2-a.^2);
        rowCorrection(i,p)=var(a,0)-var(b,0);
        predictionCorrection(i,p)=var(Data(i,3)-a,0)-var(Data(i,3)-b,0);
        bagged(i,p)=mean(b)^2-mean(a)^2;
        finiteB(i,p)=rowCorrection(i,p)/numel(a);
    end
end
close(r.draws.raw,draws,2e-10,'raw OLS bootstrap draws');
close(r.draws.zeroDeletion,zeroDraws,2e-10,'zero-deletion OLS draws');
close(r.counts,counts,0,'observation OOB counts');
close(r.deltaRaw,mean(rowRaw,1),2e-10,'equal-observation weighting');
close(r.correction,mean(rowCorrection,1),2e-10,'sample residual variance correction');
close(r.correction,mean(predictionCorrection,1),2e-10,'fixed-outcome prediction variance');
close(r.deltaCorrected,mean(bagged+finiteB,1),2e-10,'finite-bootstrap bagging identity');
assert(norm(mean(rowRaw,1)-mean(draws,1),Inf)>1e-8, ...
    'Test dataset did not distinguish observation from replicate weighting.');
assert(isequal(r.pRaw,mean(draws<=0,1)));
assert(isequal(r.lowerRaw>0,r.pRaw<o.Alpha));
assert(r.identityError<1e-10);

% Multi-indicator outcomes: mean squared coordinates, never squared mean.
rng(8193,'twister');
u=randn(N,1); v=.6*u+randn(N,1);
M=[1+u+.4*randn(N,1), -2+.8*u+.5*randn(N,1), ...
    5+1.4*v+.6*randn(N,1), -1+.7*v+.5*randn(N,1)];
MW=[1,0;1,0;0,1;0,1]; MC=MW'; MB=[0,1;0,0];
mo=struct('Boot',79,'Seed',9321);
mr=gsca_indicator_predictive(M,MW,MC,MB,mo);
assert(mr.valid,mr.reason);
assert(isequal(size(mr.aggregation),[1,2]));
close(mr.aggregation,[.5,.5],0,'equal indicator aggregation');
close(mr.deltaRaw,mean(mr.indicator.deltaRaw),1e-12,'coordinate point aggregation');
close(mr.correction,mean(mr.indicator.correction),1e-12,'coordinate variance aggregation');
rng(mo.Seed,'twister');
ib=randi(N,N,1); used=false(N,1); used(ib)=true; ix=find(~used);
f=gsca_fit(M(ib,:),MW,MC,MB,0);
mp=gsca_indicator_maps(f,MC,MB,mean(M,1),cov(M,1));
F=[ones(numel(ix),1),M(ix,:)]*mp.full;
G=[ones(numel(ix),1),M(ix,:)]*mp.removed;
right=mean(mean(G.^2-F.^2,2),1);
wrong=mean(mean(G,2).^2-mean(F,2).^2);
close(mr.draws.raw(1),right,2e-10,'mean squared block loss');
assert(abs(right-wrong)>1e-7,'Block test failed to distinguish squared means.');

% Changing all outcome units preserves decisions and scales squared loss.
scaled=M; scaled(:,3:4)=3*scaled(:,3:4);
sr=gsca_indicator_predictive(scaled,MW,MC,MB,mo);
assert(sr.valid,sr.reason);
close(sr.deltaRaw,9*mr.deltaRaw,1e-8,'raw outcome-unit change');
close(sr.correction,9*mr.correction,1e-8,'correction outcome-unit change');
close(sr.deltaCorrected,9*mr.deltaCorrected,1e-8,'corrected outcome-unit change');
assert(isequal(sr.pRaw,mr.pRaw) && isequal(sr.pShifted,mr.pShifted), ...
    'Positive common outcome scaling changed percentile decisions.');

% Replay nested sampling. The default reference moments must be recomputed
% for each outer sample, rather than frozen at the original sample moments.
no=struct('Boot',79,'Outer',3,'InnerBoot',79,'Seed',6324);
nr=gsca_indicator_predictive(Data,W0,C0,B0,no);
assert(nr.valid && nr.outerValid,'Nested indicator estimator failed.');
rng(no.Seed,'twister');
for b=1:no.Boot, randi(N,N,1); end
outerIndices=randi(N,N,no.Outer);
innerSeeds=randi(2^31-1,no.Outer,1);
assert(isequal(nr.outer.bootstrapIndices,outerIndices),'Outer sample stream changed.');
assert(isequal(nr.outer.innerSeeds,innerSeeds),'Inner seed stream changed.');
inner=struct('Boot',no.InnerBoot,'Seed',[]);
frozenDifference=0;
for a=1:no.Outer
    ib=outerIndices(:,a); sample=Data(ib,:);
    inner.Seed=innerSeeds(a);
    expected=gsca_indicator_predictive(sample,W0,C0,B0,inner);
    assert(expected.valid,expected.reason);
    close(expected.referenceMean,mean(sample,1),0,'outer reference mean');
    close(expected.referenceCov,cov(sample,1),0,'outer reference covariance');
    close(nr.outer.raw(a,:),expected.deltaRaw,2e-10,'nested raw replay');
    close(nr.outer.corrected(a,:),expected.deltaCorrected,2e-10,'nested corrected replay');
    frozen=inner; frozen.ReferenceMean=mean(Data,1); frozen.ReferenceCov=cov(Data,1);
    wrong=gsca_indicator_predictive(sample,W0,C0,B0,frozen);
    assert(wrong.valid,wrong.reason);
    frozenDifference=max(frozenDifference,norm(expected.deltaRaw-wrong.deltaRaw,Inf));
end
assert(frozenDifference>1e-8,'Nested test did not distinguish frozen reference moments.');
assert(isequal(nr.pOuterCorrected,mean(nr.outer.corrected<=0,1)));
assert(isequal(nr.lowerOuterCorrected>0,nr.pOuterCorrected<.05));
% Inner-budget sensitivity uses identical outer rows and prefix inner draws.
larger=no; larger.InnerBoot=99;
large=gsca_indicator_predictive(Data,W0,C0,B0,larger);
assert(large.valid && large.outerValid,'Larger inner budget failed.');
assert(isequal(large.outer.bootstrapIndices,nr.outer.bootstrapIndices));
assert(isequal(large.outer.innerSeeds,nr.outer.innerSeeds));
sample=Data(outerIndices(:,1),:); inner.Seed=innerSeeds(1);
short=gsca_indicator_predictive(sample,W0,C0,B0,inner);
inner.Boot=larger.InnerBoot;
long=gsca_indicator_predictive(sample,W0,C0,B0,inner);
assert(isequal(short.draws.raw,long.draws.raw(1:no.InnerBoot,:)), ...
    'Increasing inner budget changed its existing bootstrap prefix.');
fprintf('test_gsca_indicator_predictive: independent OLS, row weighting, block loss, units, and nested checks passed.\n');
end

function close(actual,expected,tolerance,label)
assert(isequal(size(actual),size(expected)),['Shape mismatch: ',label]);
assert(all(isfinite(actual(:))) && norm(actual(:)-expected(:),Inf)<=tolerance, ...
    ['Numerical mismatch: ',label]);
end

