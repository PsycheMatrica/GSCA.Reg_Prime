function Results = Run_Simulation_Correction_v3(mode,overrides)
%RUN_SIMULATION_CORRECTION_V3 Finite-N GSCA predictive-contribution study.
% Replaces the old coefficient-zero/nonnull classification with independently
% estimated risk truth for the specified fitting and replacement algorithms.
%
%   R=Run_Simulation_Correction_v3('smoke');
%   R=Run_Simulation_Correction_v3('ci');       % 12 conditions, R=100
%   R=Run_Simulation_Correction_v3('pilot');    % larger two-DGP comparison
%   R=Run_Simulation_Correction_v3('sensitivity');
%   R=Run_Simulation_Correction_v3('nested_smoke'); % execution check only
%   R=Run_Simulation_Correction_v3('ci',struct('Replications',20));
%
% All intervals evaluated here are ONE-SIDED lower confidence bounds.
% p<Alpha matches empirical inverse-CDF lower-bound rejection.
% "boundary" is a Monte Carlo calibrated candidate, not an exact null.
% Truth is calculated from independent training samples and analytic test-
% risk integration. Failed truth fits invalidate truth and calibration;
% completer-only diagnostics are retained but never relabeled as truth.
%
% Requires gsca_fit, gsca_risk_maps, gsca_exact_risk, gsca_predictive,
% and gsca_correction_dgp from this research directory. No private DGP_SEM
% package, gramm, parallel pool, or Statistics Toolbox is required.
if nargin<1 || isempty(mode), mode='smoke'; end
if nargin<2, overrides=struct; end
cfg=defaults(char(mode));
names=fieldnames(overrides);
for k=1:numel(names)
    if ~isfield(cfg,names{k})
        error('gsca:UnknownSimulationOption','Unknown option: %s',names{k});
    end
    cfg.(names{k})=overrides.(names{k});
end
validate_config(cfg);
if ~exist(cfg.OutputDir,'dir'), mkdir(cfg.OutputDir); end
records=repmat(blank_record(),0,1);
summaries=repmat(blank_summary(),0,1);
truthRows=repmat(blank_truth(),0,1);
calibrations=cell(0,1);
fitOptions=struct('MaxIter',cfg.MaxIter,'Tol',cfg.Tol,'MinRcond',cfg.MinRcond);
regime=0; condition=0; started=tic;
for ni=1:numel(cfg.Ns)
    for ji=1:numel(cfg.Jps)
        for li=1:numel(cfg.Lambdas)
            for di=1:numel(cfg.DGPs)
                regime=regime+1;
                n=cfg.Ns(ni); jp=cfg.Jps(ji); lambda=cfg.Lambdas(li);
                dgp=cfg.DGPs{di};
                calSeed=cfg.Seed+1000000+regime*20000;
                calibration=calibrate(n,jp,lambda,dgp,cfg,fitOptions,calSeed);
                calibration.N=n; calibration.Jp=jp;
                calibration.Lambda=lambda; calibration.DGP=dgp;
                calibrations{regime,1}=calibration;
                fprintf('CAL N=%d Jp=%d lambda=%.4g DGP=%s beta=%.6g valid=%d\n', ...
                    n,jp,lambda,dgp,calibration.beta,calibration.valid);
                for ei=1:numel(cfg.Effects)
                    effect=cfg.Effects{ei};
                    if strcmp(effect,'zero'), beta=0;
                    elseif strcmp(effect,'positive'), beta=cfg.PositiveBeta;
                    else, beta=calibration.beta;
                    end
                    condition=condition+1;
                    truthSeed=cfg.Seed+10000000+condition*20000;
                    if isfinite(beta)
                        truth=estimate_truth(n,jp,beta,lambda,dgp,cfg.TruthReplications, ...
                            truthSeed,fitOptions);
                    else
                        truth=invalid_truth(cfg.TruthReplications,'Calibration failed.');
                    end
                    tr=blank_truth();
                    tr.Condition=condition; tr.N=n; tr.Jp=jp; tr.Lambda=lambda;
                    tr.DGP=dgp; tr.Effect=effect; tr.Beta=beta;
                    tr.CalibrationValid=calibration.valid;
                    tr.TruthValid=truth.valid; tr.TruthMean=truth.mean;
                    tr.TruthMCSE=truth.mcse; tr.TruthFailed=truth.failed;
                    tr.TruthCompleted=truth.completed; tr.TruthRequested=truth.requested;
                    tr.CompleterMean=truth.completerMean;
                    tr.CompleterMCSE=truth.completerMCSE;
                    tr.Classification=classify_truth(truth,effect);
                    truthRows(end+1,1)=tr; %#ok<AGROW>
                    for ri=1:numel(cfg.ReferenceModes)
                        reference=cfg.ReferenceModes{ri};
                        firstRecord=numel(records)+1;
                        for rep=1:cfg.Replications
                            row=blank_record();
                            row.Condition=condition; row.N=n; row.Jp=jp;
                            row.Lambda=lambda; row.DGP=dgp; row.Effect=effect;
                            row.Beta=beta; row.Reference=reference; row.Rep=rep;
                            if isfinite(beta)
                                % Paired source datasets and bootstrap seeds across
                                % oracle/sample-reference comparisons.
                                dataSeed=cfg.Seed+20000000+condition*20000+rep;
                                bootSeed=cfg.Seed+30000000+condition*20000+rep;
                                [data,model]=gsca_correction_dgp(n,jp,beta,dgp,dataSeed);
                                opt=struct('Boot',cfg.Boot,'Outer',cfg.Outer, ...
                                    'InnerBoot',cfg.InnerBoot,'Lambda',lambda, ...
                                    'Alpha',cfg.Alpha,'Seed',bootSeed, ...
                                    'MaxIter',cfg.MaxIter,'Tol',cfg.Tol, ...
                                    'MinRcond',cfg.MinRcond);
                                if strcmp(reference,'oracle')
                                    opt.ReferenceMean=model.mean;
                                    opt.ReferenceCov=model.covariance;
                                end
                                out=gsca_predictive(data,model.W0,model.C0,model.B0,opt);
                                row.Valid=logical(out.valid);
                                if isnumeric(out.failures), row.FitFailures=sum(out.failures(:)); end
                                if row.Valid
                                    p=model.focalIndex;
                                    row.Raw=out.deltaRaw(p);
                                    row.Corrected=out.deltaCorrected(p);
                                    row.Correction=out.correction(p);
                                    row.pLegacy=out.pLegacy(p);
                                    row.pRaw=out.pRaw(p);
                                    row.pShifted=out.pShifted(p);
                                    row.LowerRaw=out.lowerRaw(p);
                                    row.LowerShifted=out.lowerShifted(p);
                                    row.pOuterRaw=out.pOuterRaw(p);
                                    row.pOuterCorrected=out.pOuterCorrected(p);
                                    row.LowerOuterRaw=out.lowerOuterRaw(p);
                                    row.LowerOuterCorrected=out.lowerOuterCorrected(p);
                                end
                            end
                            records(end+1,1)=row; %#ok<AGROW>
                        end
                        block=records(firstRecord:end);
                        summary=summarize(block,tr,cfg);
                        summaries(end+1,1)=summary; %#ok<AGROW>
                        fprintf(['SIM N=%d Jp=%d lambda=%.4g DGP=%s effect=%s ref=%s ' ...
                            'beta=%.5f truth=%.6f (MCSE %.3g) valid=%d/%d ' ...
                            'bias raw=%.5f corrected=%.5f rej legacy/raw/shift/outer=' ...
                            '%.3f/%.3f/%.3f/%.3f\n'], ...
                            n,jp,lambda,dgp,effect,reference,beta,tr.TruthMean,tr.TruthMCSE, ...
                            summary.Completed,summary.Attempted,summary.BiasRaw, ...
                            summary.BiasCorrected,summary.RejectLegacy,summary.RejectRaw, ...
                            summary.RejectShifted,summary.RejectOuterCorrected);
                        % Checkpoint after each completed condition/reference.
                        Results=assemble(cfg,records,summaries,truthRows,calibrations,toc(started));
                        save(fullfile(cfg.OutputDir,'simulation_results.mat'),'Results','-v7');
                        writetable(Results.Summary,fullfile(cfg.OutputDir,'summary.csv'));
                        writetable(Results.Truth,fullfile(cfg.OutputDir,'truth.csv'));
                        writetable(Results.Replicates,fullfile(cfg.OutputDir,'replicates.csv'));
                    end
                end
            end
        end
    end
end
Results=assemble(cfg,records,summaries,truthRows,calibrations,toc(started));
save(fullfile(cfg.OutputDir,'simulation_results.mat'),'Results','-v7');
fprintf('Completed %d condition/reference rows in %.1f seconds. Output: %s\n', ...
    height(Results.Summary),Results.ElapsedSeconds,cfg.OutputDir);
disp(Results.Summary);
end

function cfg=defaults(mode)
cfg.Mode=mode;
cfg.Ns=80; cfg.Jps=[1,3]; cfg.Lambdas=[0,.1];
cfg.DGPs={'normal'}; cfg.Effects={'zero','boundary','positive'};
cfg.ReferenceModes={'oracle'};
cfg.Replications=100; cfg.CalReplications=100; cfg.TruthReplications=300;
cfg.Boot=199; cfg.Outer=0; cfg.InnerBoot=99;
cfg.CalIterations=8; cfg.BoundaryUpper=.55; cfg.PositiveBeta=.30;
cfg.Alpha=.05; cfg.Seed=20261106; cfg.MaxIter=1000;
cfg.Tol=1e-8; cfg.MinRcond=1e-10;
switch mode
    case 'smoke'
        cfg.Replications=2; cfg.CalReplications=12; cfg.TruthReplications=12;
        cfg.Boot=49; cfg.Lambdas=0; cfg.CalIterations=4;
    case 'ci'
        % Moderate pilot; uncertainty at nominal .05 is substantial with R=100.
    case 'pilot'
        cfg.Ns=[80,200]; cfg.DGPs={'normal','heteroscedastic'};
        cfg.ReferenceModes={'oracle','sample'};
        cfg.Replications=300; cfg.CalReplications=500; cfg.TruthReplications=1500;
        cfg.Boot=499; cfg.CalIterations=10;
    case 'sensitivity'
        cfg.Ns=200; cfg.Jps=3;
        cfg.DGPs={'normal','heteroscedastic'};
    case 'nested_smoke'
        cfg.Ns=40; cfg.Jps=1; cfg.Lambdas=0;
        cfg.Effects={'boundary','positive'};
        cfg.Replications=2; cfg.CalReplications=8; cfg.TruthReplications=12;
        cfg.CalIterations=4; cfg.Boot=49; cfg.Outer=19; cfg.InnerBoot=49;
    otherwise
        error('gsca:InvalidSimulationMode','Unknown mode: %s',mode);
end
cfg.OutputDir=fullfile(fileparts(mfilename('fullpath')),'results',mode);
end

function validate_config(c)
validateattributes(c.Ns,{'numeric'},{'vector','integer','>=',8});
validateattributes(c.Jps,{'numeric'},{'vector','integer','>=',1});
validateattributes(c.Lambdas,{'numeric'},{'vector','real','finite','nonnegative'});
ints={'Replications','CalReplications','TruthReplications','CalIterations','Boot','InnerBoot','MaxIter'};
for k=1:numel(ints), validateattributes(c.(ints{k}),{'numeric'},{'scalar','integer','>=',1}); end
validateattributes(c.Outer,{'numeric'},{'scalar','integer','nonnegative'});
validateattributes(c.Alpha,{'numeric'},{'scalar','>',0,'<',1});
validateattributes(c.Seed,{'numeric'},{'scalar','integer','nonnegative'});
if c.CalReplications<2 || c.TruthReplications<2 || c.Replications<2
    error('gsca:InvalidSimulationCount','Use at least two datasets for MCSE.');
end
if ~iscell(c.DGPs) || ~all(ismember(c.DGPs,{'normal','heteroscedastic'})) || ...
        ~iscell(c.ReferenceModes) || ~all(ismember(c.ReferenceModes,{'oracle','sample'})) || ...
        ~iscell(c.Effects) || ~all(ismember(c.Effects,{'zero','boundary','positive'}))
    error('gsca:InvalidSimulationCategory','Invalid DGP, reference, or effect category.');
end
end

function cal=calibrate(n,jp,lambda,dgp,cfg,fitOptions,seed)
cal=struct('valid',false,'beta',NaN,'reason','','trace',[], ...
    'calibrationReplications',cfg.CalReplications,'iterations',0);
lo=0; hi=cfg.BoundaryUpper;
tlo=estimate_truth(n,jp,lo,lambda,dgp,cfg.CalReplications,seed,fitOptions);
thi=estimate_truth(n,jp,hi,lambda,dgp,cfg.CalReplications,seed,fitOptions);
cal.trace=[lo,tlo.mean,tlo.mcse,tlo.failed;hi,thi.mean,thi.mcse,thi.failed];
if ~tlo.valid || ~thi.valid
    cal.reason='A training fit failed; the risk boundary is undefined for this calibration.';
    return
end
if tlo.mean>0 || thi.mean<0
    cal.reason='The Monte Carlo mean risk does not bracket zero.';
    return
end
for k=1:cfg.CalIterations
    mid=(lo+hi)/2;
    tm=estimate_truth(n,jp,mid,lambda,dgp,cfg.CalReplications,seed,fitOptions);
    cal.trace(end+1,:)=[mid,tm.mean,tm.mcse,tm.failed]; %#ok<AGROW>
    cal.iterations=k;
    if ~tm.valid
        cal.reason='A training fit failed during boundary calibration.';
        return
    end
    if tm.mean<=0, lo=mid; else, hi=mid; end
end
cal.beta=(lo+hi)/2;
cal.valid=true;
cal.bracket=[lo,hi];
cal.reason='Monte Carlo mean root; independently validate its risk and MCSE.';
end

function t=estimate_truth(n,jp,beta,lambda,dgp,R,seed,fitOptions)
vals=NaN(R,1); reason='';
for rep=1:R
    [data,model]=gsca_correction_dgp(n,jp,beta,dgp,seed+rep);
    fit=gsca_fit(data,model.W0,model.C0,model.B0,lambda,fitOptions);
    if ~fit.converged
        if isempty(reason), reason=fit.reason; end
        continue
    end
    maps=gsca_risk_maps(fit,model.C0,model.B0,model.mean,model.covariance);
    if ~maps.valid
        if isempty(reason), reason='Invalid residual map.'; end
        continue
    end
    risk=gsca_exact_risk(maps,model.mean,model.covariance);
    vals(rep)=risk.deltaEPE(model.focalIndex);
end
good=isfinite(vals);
t=invalid_truth(R,reason);
t.values=vals; t.completed=sum(good); t.failed=R-t.completed;
if t.completed>0, t.completerMean=mean(vals(good)); end
if t.completed>1, t.completerMCSE=std(vals(good),0)/sqrt(t.completed); end
t.valid=all(good);
if t.valid
    t.mean=t.completerMean; t.mcse=t.completerMCSE; t.reason='';
end
end

function t=invalid_truth(R,reason)
t=struct('valid',false,'mean',NaN,'mcse',NaN,'failed',R, ...
    'completed',0,'requested',R,'completerMean',NaN,'completerMCSE',NaN, ...
    'values',NaN(R,1),'reason',reason);
end

function label=classify_truth(t,effect)
if ~t.valid, label='undefined_fit_failure'; return; end
if t.mean-1.96*t.mcse>0
    label='positive_MC_supported';
elseif t.mean+1.96*t.mcse<0
    label='null_interior_MC_supported';
elseif strcmp(effect,'boundary')
    label='boundary_MC_unresolved';
else
    label='sign_MC_unresolved';
end
end

function r=blank_record()
r=struct('Condition',0,'N',0,'Jp',0,'Lambda',0,'DGP','','Effect','', ...
    'Beta',NaN,'Reference','','Rep',0,'Valid',false,'FitFailures',NaN, ...
    'Raw',NaN,'Corrected',NaN,'Correction',NaN, ...
    'pLegacy',NaN,'pRaw',NaN,'pShifted',NaN, ...
    'LowerRaw',NaN,'LowerShifted',NaN, ...
    'pOuterRaw',NaN,'pOuterCorrected',NaN, ...
    'LowerOuterRaw',NaN,'LowerOuterCorrected',NaN);
end

function t=blank_truth()
t=struct('Condition',0,'N',0,'Jp',0,'Lambda',0,'DGP','','Effect','', ...
    'Beta',NaN,'CalibrationValid',false,'TruthValid',false,'TruthMean',NaN, ...
    'TruthMCSE',NaN,'TruthFailed',0,'TruthCompleted',0,'TruthRequested',0, ...
    'CompleterMean',NaN,'CompleterMCSE',NaN,'Classification','');
end

function s=blank_summary()
s=struct('Condition',0,'N',0,'Jp',0,'Lambda',0,'DGP','','Effect','', ...
    'Beta',NaN,'Reference','','Classification','','TruthMean',NaN, ...
    'TruthMCSE',NaN,'TruthValid',false,'Attempted',0,'Completed',0, ...
    'DatasetFailures',0,'BootstrapFitFailures',0, ...
    'MeanRaw',NaN,'MeanCorrected',NaN,'MeanCorrection',NaN, ...
    'BiasRaw',NaN,'BiasCorrected',NaN,'RMSERaw',NaN,'RMSECorrected',NaN, ...
    'RejectLegacy',NaN,'MCSELegacy',NaN,'RejectRaw',NaN,'MCSERaw',NaN, ...
    'RejectShifted',NaN,'MCSEShifted',NaN, ...
    'RejectOuterRaw',NaN,'MCSEOuterRaw',NaN,'OuterRawCompleted',0, ...
    'RejectOuterCorrected',NaN,'MCSEOuterCorrected',NaN,'OuterCorrectedCompleted',0, ...
    'LowerCoverageRaw',NaN,'LowerCoverageShifted',NaN, ...
    'LowerCoverageOuterRaw',NaN,'LowerCoverageOuterCorrected',NaN);
end

function s=summarize(rows,t,cfg)
s=blank_summary();
copy={'Condition','N','Jp','Lambda','DGP','Effect','Beta'};
for k=1:numel(copy), s.(copy{k})=t.(copy{k}); end
s.Reference=rows(1).Reference; s.Classification=t.Classification;
s.TruthMean=t.TruthMean; s.TruthMCSE=t.TruthMCSE; s.TruthValid=t.TruthValid;
s.Attempted=numel(rows); good=[rows.Valid]; s.Completed=sum(good);
s.DatasetFailures=s.Attempted-s.Completed;
fail=[rows.FitFailures]; s.BootstrapFitFailures=sum(fail(isfinite(fail)));
if s.Completed==0, return; end
r=rows(good); raw=[r.Raw]; corrected=[r.Corrected];
s.MeanRaw=mean(raw); s.MeanCorrected=mean(corrected);
s.MeanCorrection=mean([r.Correction]);
if t.TruthValid
    s.BiasRaw=mean(raw-t.TruthMean); s.BiasCorrected=mean(corrected-t.TruthMean);
    s.RMSERaw=sqrt(mean((raw-t.TruthMean).^2));
    s.RMSECorrected=sqrt(mean((corrected-t.TruthMean).^2));
    s.LowerCoverageRaw=coverage([r.LowerRaw],t.TruthMean);
    s.LowerCoverageShifted=coverage([r.LowerShifted],t.TruthMean);
    s.LowerCoverageOuterRaw=coverage([r.LowerOuterRaw],t.TruthMean);
    s.LowerCoverageOuterCorrected=coverage([r.LowerOuterCorrected],t.TruthMean);
end
[s.RejectLegacy,s.MCSELegacy]=rate([r.pLegacy],cfg.Alpha);
[s.RejectRaw,s.MCSERaw]=rate([r.pRaw],cfg.Alpha);
[s.RejectShifted,s.MCSEShifted]=rate([r.pShifted],cfg.Alpha);
[s.RejectOuterRaw,s.MCSEOuterRaw,s.OuterRawCompleted]=rate([r.pOuterRaw],cfg.Alpha);
[s.RejectOuterCorrected,s.MCSEOuterCorrected,s.OuterCorrectedCompleted]= ...
    rate([r.pOuterCorrected],cfg.Alpha);
end

function [p,se,n]=rate(values,alpha)
values=values(isfinite(values)); n=numel(values); p=NaN; se=NaN;
if n==0, return; end
p=mean(values<alpha); se=sqrt(p*(1-p)/n);
end

function p=coverage(lower,truth)
lower=lower(isfinite(lower)); p=NaN;
if ~isempty(lower), p=mean(lower<=truth); end
end

function R=assemble(cfg,records,summaries,truthRows,calibrations,elapsed)
R.Config=cfg; R.Summary=struct2table(summaries);
R.Truth=struct2table(truthRows); R.Replicates=struct2table(records);
R.Calibration=calibrations; R.ElapsedSeconds=elapsed;
R.Notes={ ...
    'Delta EPE_p: unnormalized expected squared-error difference.', ...
    'Oracle finite-N target uses learned fitted weights, fixed lambda, and population projection conditional on those weights.', ...
    'Structural loss is algorithm-dependent because fitted endogenous component weights vary.', ...
    'Sample-reference estimates approximate the oracle target; their projection and scale error remain part of their bias.', ...
    'All coverage values are one-sided lower-bound coverage, using independently estimated truth.', ...
    'A boundary label means a Monte Carlo calibrated candidate, not an exact predictive null.', ...
    'TruthMean is NaN after any failed training fit. CompleterMean is diagnostic only.', ...
    'Rejection rates condition on valid datasets; always inspect failure counts.', ...
    'Outer bootstrap entries are NaN when Outer=0; nested_smoke is an execution check, not a calibration study.', ...
    'R=100 has rejection MCSE about .022 at p=.05; zero observed rejection does not establish zero population rejection.'};
end
