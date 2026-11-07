function out = gsca_predictive(Data,W0,C0,B0,opts)
%GSCA_PREDICTIVE Experimental paired-residual OOB correction for Delta EPE_p.
% Fixed prespecified scalar Lambda only. This is a research implementation:
% neither variance correction nor nested percentile guarantees calibrated
% finite-N predictive inference. Production GSCA_Reg is unchanged.
% Edges: find(C0), then find(B0), in MATLAB column-major order.
% ReferenceMean/ReferenceCov: externally specified population X moments, or
% empty (default) for full-sample X plug-in moments. These are recomputed
% inside each outer bootstrap when not supplied.
% All intervals below are ONE-SIDED lower bounds. Empirical inverse-CDF
% quantiles use order ceil(alpha*B); reject iff lower>0, equivalently p<alpha.
% Any failed inner fit makes the point estimator invalid (no silent redraw).
if nargin<5, opts=struct(); end
opts=defaults(opts);
validateattributes(Data,{'numeric'},{'2d','real','finite','nonempty'});
N=size(Data,1); E=nnz(C0)+nnz(B0);
if N<4 || E==0, error('gsca:input','Need N>=4 and a specified edge.'); end
if ~isempty(opts.Seed)
    previous=rng; cleanup=onCleanup(@()rng(previous)); %#ok<NASGU>
    rng(opts.Seed,'twister');
end
out=struct('valid',false,'reason','','failures',0,'outerValid',false);
names={'deltaRaw','correction','deltaCorrected','legacyRaw','pLegacy',...
    'pRaw','pShifted','lowerLegacy','lowerRaw','lowerShifted',...
    'pOuterRaw','pOuterCorrected','lowerOuterRaw','lowerOuterCorrected'};
for j=1:numel(names), out.(names{j})=nan(1,E); end
out.options=opts;
refMu=opts.ReferenceMean; refSigma=opts.ReferenceCov;
if isempty(refMu), refMu=mean(Data,1); end
if isempty(refSigma), refSigma=cov(Data,1); end
out.referenceMean=refMu; out.referenceCov=refSigma;
if isempty(opts.ReferenceCov)
    out.referenceMode='full-data indicator-moment plug-in (including target scale)';
else
    out.referenceMode='external covariance; reference mean follows ReferenceMean option';
end
initial=gsca_fit(Data,W0,C0,B0,opts.Lambda,opts);
out.fit=initial;
if ~initial.converged, out.reason=['Original fit: ' initial.reason]; return; end
base=gsca_risk_maps(initial,C0,B0,refMu,refSigma);
out.edges=struct('kind',{base.kind},'source',base.source,...
    'target',base.target,'index',base.index);
if ~base.valid, out.reason='Original residual map is undefined.'; return; end
sumF=zeros(N,E); sumG=sumF; sumF2=sumF; sumG2=sumF;
counts=zeros(N,1); draws=nan(opts.Boot,E); zeroDraws=draws;
fitReasons=cell(opts.Boot,1);
for b=1:opts.Boot
    idx=randi(N,N,1); used=false(N,1); used(idx)=true; oob=find(~used);
    if isempty(oob)
        out.failures=out.failures+1; fitReasons{b}='empty OOB'; continue
    end
    model=gsca_fit(Data(idx,:),W0,C0,B0,opts.Lambda,opts);
    if ~model.converged
        out.failures=out.failures+1; fitReasons{b}=model.reason; continue
    end
    mp=gsca_risk_maps(model,C0,B0,refMu,refSigma);
    if ~mp.valid
        out.failures=out.failures+1; fitReasons{b}='undefined residual map'; continue
    end
    X=[ones(numel(oob),1),Data(oob,:)];
    ef=X*mp.full; eg=X*mp.removed; ez=X*mp.zeroRemoved;
    if any(~isfinite([ef(:);eg(:);ez(:)]))
        out.failures=out.failures+1; fitReasons{b}='nonfinite loss'; continue
    end
    losses=[ef.^2,eg.^2,ez.^2];
    if any(~isfinite(losses(:)))
        out.failures=out.failures+1; fitReasons{b}='nonfinite squared loss'; continue
    end
    counts(oob)=counts(oob)+1;
    sumF(oob,:)=sumF(oob,:)+ef; sumG(oob,:)=sumG(oob,:)+eg;
    sumF2(oob,:)=sumF2(oob,:)+ef.^2; sumG2(oob,:)=sumG2(oob,:)+eg.^2;
    draws(b,:)=mean(eg.^2-ef.^2,1);
    zeroDraws(b,:)=mean(ez.^2-ef.^2,1);
end
out.counts=counts; out.minOOBCount=min(counts);
out.failureReasons=fitReasons;
out.draws=struct('raw',draws,'zeroDeletion',zeroDraws);
if out.failures>0 || any(counts<2)
    out.reason='At least one failed draw or fewer than two OOB evaluations; no inference returned.';
    return
end
if any(~isfinite([sumF(:);sumG(:);sumF2(:);sumG2(:)]))
    out.reason='Nonfinite accumulated residual moments.'; return
end
centerF=sumF2-sumF.^2./counts; centerG=sumG2-sumG.^2./counts;
if any(~isfinite([centerF(:);centerG(:)]))
    out.reason='Nonfinite centered residual moments.'; return
end
rawByRow=(sumG2-sumF2)./counts;
sF2=max(0,centerF./(counts-1));
sG2=max(0,centerG./(counts-1));
correctionByRow=sF2-sG2;
out.deltaRaw=mean(rawByRow,1);
out.correction=mean(correctionByRow,1);
out.deltaCorrected=out.deltaRaw+out.correction;
out.legacyRaw=mean(zeroDraws,1);
out.rawReplicateMean=mean(draws,1);
out.baggedLossDifference=mean((sumG./counts).^2-(sumF./counts).^2,1);
out.finiteBootAdjustment=mean(correctionByRow./counts,1);
out.identityError=max(abs(out.deltaCorrected-out.baggedLossDifference-out.finiteBootAdjustment));
if any(~isfinite([out.deltaRaw,out.correction,out.deltaCorrected,...
        out.legacyRaw,out.rawReplicateMean,out.baggedLossDifference,...
        out.finiteBootAdjustment,out.identityError]))
    out.reason='Nonfinite risk summary.'; return
end
[out.pRaw,out.lowerRaw]=gsca_percentile(draws,opts.Alpha);
[out.pLegacy,out.lowerLegacy]=gsca_percentile(zeroDraws,opts.Alpha);
[out.pShifted,out.lowerShifted]=gsca_percentile(draws+out.correction,opts.Alpha);
out.valid=true; out.reason='ok';
out.rejectionRule='lower > 0, equivalent to empirical p < Alpha (strict at ties)';
out.interpretation='Residual-variance/bagging correction candidate; not proven finite-N GSCA bias correction.';
if opts.Outer>0
    outerRaw=nan(opts.Outer,E); outerCorrected=outerRaw;
    outerReasons=cell(opts.Outer,1);
    inner=opts; inner.Outer=0; inner.Boot=opts.InnerBoot; inner.Seed=[];
    % Positions, including duplicates, are the units of ordinary empirical
    % pairs bootstrap. Source-ID exclusion would define a different procedure.
    for a=1:opts.Outer
        idx=randi(N,N,1);
        oa=gsca_predictive(Data(idx,:),W0,C0,B0,inner);
        if oa.valid
            outerRaw(a,:)=oa.deltaRaw; outerCorrected(a,:)=oa.deltaCorrected;
        else
            outerReasons{a}=oa.reason;
        end
    end
    out.outer=struct('raw',outerRaw,'corrected',outerCorrected,...
        'failureReasons',{outerReasons},'failed',sum(any(~isfinite(outerCorrected),2)));
    out.outerValid=out.outer.failed==0;
    if out.outerValid
        [out.pOuterRaw,out.lowerOuterRaw]=gsca_percentile(outerRaw,opts.Alpha);
        [out.pOuterCorrected,out.lowerOuterCorrected]=gsca_percentile(outerCorrected,opts.Alpha);
    end
end
end

function o=defaults(o)
if ~isstruct(o)||~isscalar(o), error('gsca:input','opts must be a scalar struct.'); end
d=struct('Boot',499,'Outer',0,'InnerBoot',[],'Lambda',0,'Alpha',.05,...
    'Seed',[],'ReferenceMean',[],'ReferenceCov',[],...
    'MaxIter',1000,'Tol',1e-8,'MinRcond',1e-10);
f=fieldnames(o);
for i=1:numel(f)
    if ~isfield(d,f{i}), error('gsca:input','Unknown option %s.',f{i}); end
end
f=fieldnames(d);
for i=1:numel(f)
    if ~isfield(o,f{i}), o.(f{i})=d.(f{i}); end
end
if isempty(o.InnerBoot), o.InnerBoot=o.Boot; end
validateattributes(o.Boot,{'numeric'},{'scalar','integer','>=',2});
validateattributes(o.Outer,{'numeric'},{'scalar','integer','>=',0});
validateattributes(o.InnerBoot,{'numeric'},{'scalar','integer','>=',2});
validateattributes(o.Lambda,{'numeric'},{'scalar','real','finite','nonnegative'});
validateattributes(o.Alpha,{'numeric'},{'scalar','>',0,'<',1});
validateattributes(o.MaxIter,{'numeric'},{'scalar','integer','positive'});
validateattributes(o.Tol,{'numeric'},{'scalar','real','finite','positive'});
validateattributes(o.MinRcond,{'numeric'},{'scalar','>',0,'<',1});
if ~isempty(o.Seed), validateattributes(o.Seed,{'numeric'},{'scalar','integer','nonnegative','finite'}); end
end
