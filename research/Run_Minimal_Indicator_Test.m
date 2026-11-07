function Results = Run_Minimal_Indicator_Test(N,options)
%RUN_MINIMAL_INDICATOR_TEST Original GSCA fits and one in-bag/OOB bootstrap.
% Only the evaluation changes in GSCA_Reg are used: in-bag auxiliary
% replacement and fixed observed target-indicator loss for B paths.
% No variance correction, nested bootstrap, refitting learner, or risk-truth
% calibration is added here. Original TABLE.B columns 1:4 are retained.
%
%   Run_Minimal_Indicator_Test(80)
%   Run_Minimal_Indicator_Test(200)
% Optional options: Replications, Boot, Seed, OutputDir.
% The six conditions are Lambda={0,.1} by beta_1={0,.1,.3}, Jp=3.
% Rejection means the ORIGINAL TABLE.B(:,5) lower bound exceeds zero.
% These are coefficient conditions, not calibrated predictive-null cases.
if nargin<1, N=80; end
if nargin<2, options=struct; end
validateattributes(N,{'numeric'},{'scalar','integer'});
if ~ismember(N,[80,200]), error('gsca:MinimalN','N must be 80 or 200.'); end
cfg=struct('N',N,'Jp',3,'Lambdas',[0,.1],'Betas',[0,.1,.3], ...
    'Replications',100,'Boot',199,'Seed',20261106, ...
    'OutputDir',fullfile(fileparts(mfilename('fullpath')), ...
    'results_minimal',sprintf('N%d',N)));
allowed={'Replications','Boot','Seed','OutputDir'};
names=fieldnames(options);
for k=1:numel(names)
    if ~ismember(names{k},allowed), error('gsca:MinimalOption','Unknown option %s.',names{k}); end
    cfg.(names{k})=options.(names{k});
end
validateattributes(cfg.Replications,{'numeric'},{'scalar','integer','>=',1,'<=',9999});
validateattributes(cfg.Boot,{'numeric'},{'scalar','integer','>=',100});
validateattributes(cfg.Seed,{'numeric'},{'scalar','integer','nonnegative'});
previous=rng; restore=onCleanup(@()rng(previous)); %#ok<NASGU>
if ~exist(cfg.OutputDir,'dir'), mkdir(cfg.OutputDir); end
rows=repmat(empty_row(),0,1); summary=repmat(empty_summary(),0,1);
parameterTables=cell(numel(cfg.Lambdas)*numel(cfg.Betas),cfg.Replications);
condition=0; started=tic;
for li=1:numel(cfg.Lambdas)
    for bi=1:numel(cfg.Betas)
        condition=condition+1; lambda=cfg.Lambdas(li); beta=cfg.Betas(bi);
        first=numel(rows)+1;
        for rep=1:cfg.Replications
            row=empty_row(); row.N=N; row.Jp=cfg.Jp;
            row.Lambda=lambda; row.Beta=beta; row.Rep=rep;
            row.DatasetIndex=((N==200)*6+condition-1)*10000+rep;
            row.DataSeed=cfg.Seed+row.DatasetIndex;
            row.BootstrapSeed=row.DataSeed+1000000;
            try
                [data,model]=gsca_correction_dgp(N,cfg.Jp,beta,'normal',row.DataSeed);
                % Separate fixed data-generation and bootstrap RNG seeds.
                rng(row.BootstrapSeed,'twister');
                fit=GSCA_Reg(data,model.W0,model.C0,model.B0,lambda, ...
                    [1,5],cfg.Boot,1000,1e-8,false);
                [source,target]=find(model.B0);
                focal=find(source==1 & target==3);
                assert(isscalar(focal),'Expected one focal path 1->3.');
                parameterTables{condition,rep}=fit.TABLE.B;
                tableRow=fit.TABLE.B(focal,:);
                row.PathEstimate=tableRow(1); row.PathBootstrapSE=tableRow(2);
                row.PathCI025=tableRow(3); row.PathCI975=tableRow(4);
                row.IndicatorLower95=tableRow(5); row.IndicatorP=tableRow(6);
                row.IndicatorDeltaMean=mean(fit.ETC.delOPE_b_Boot(focal,:));
                row.LegacyComponentDeltaMean=mean(fit.ETC.delOPE_b_legacy_Boot(focal,:));
                row.IndicatorEPE=fit.INI.indicatorOPE_byComponent(3);
                required=[tableRow(1:6),row.IndicatorDeltaMean, ...
                    row.LegacyComponentDeltaMean,row.IndicatorEPE];
                if any(~isfinite(required))
                    error('gsca:NonfiniteMinimalResult','A required original GSCA result is nonfinite.');
                end
                row.IndicatorReject=(row.IndicatorLower95>0);
                row.Success=true;
            catch err
                row.FailureIdentifier=err.identifier;
                row.FailureReason=err.message;
            end
            rows(end+1,1)=row; %#ok<AGROW>
            if mod(rep,25)==0 || rep==cfg.Replications
                fprintf('MINIMAL N=%d lambda=%.1f beta=%.1f rep=%d/%d success=%d\n', ...
                    N,lambda,beta,rep,cfg.Replications,row.Success);
            end
        end
        s=summarize(rows(first:end));
        summary(end+1,1)=s; %#ok<AGROW>
        Results=assemble(cfg,rows,summary,parameterTables,toc(started));
        save_outputs(Results);
        fprintf('RESULT N=%d lambda=%.1f beta=%.1f meanDelta=%.6f EPE=%.6f rejection=%.3f MCSE=%.4f completed=%d/%d\n', ...
            N,lambda,beta,s.MeanIndicatorDelta,s.MeanIndicatorEPE,s.RejectionRate, ...
            s.RejectionMCSE,s.Completed,s.Attempted);
    end
end
Results=assemble(cfg,rows,summary,parameterTables,toc(started));
save_outputs(Results);
disp(Results.Summary);
fprintf('JSON_SUMMARY_BEGIN\n');
disp(jsonencode(table2struct(Results.Summary)));
fprintf('JSON_SUMMARY_END\n');
end

function r=empty_row()
r=struct('N',0,'Jp',3,'Lambda',NaN,'Beta',NaN,'Rep',0, ...
    'DatasetIndex',0,'DataSeed',0,'BootstrapSeed',0,'Success',false, ...
    'FailureIdentifier','','FailureReason','', ...
    'PathEstimate',NaN,'PathBootstrapSE',NaN,'PathCI025',NaN,'PathCI975',NaN, ...
    'IndicatorDeltaMean',NaN,'IndicatorLower95',NaN,'IndicatorP',NaN, ...
    'IndicatorReject',false,'LegacyComponentDeltaMean',NaN,'IndicatorEPE',NaN);
end

function s=empty_summary()
s=struct('N',0,'Jp',3,'Lambda',NaN,'Beta',NaN, ...
    'Attempted',0,'Completed',0,'Failed',0,'MeanPathEstimate',NaN, ...
    'MeanPathBootstrapSE',NaN,'MeanPathCI025',NaN,'MeanPathCI975',NaN, ...
    'MeanIndicatorDelta',NaN,'SDIndicatorDelta',NaN,'MCSEIndicatorDelta',NaN, ...
    'MeanIndicatorLower95',NaN,'MeanIndicatorP',NaN, ...
    'RejectionRate',NaN,'RejectionMCSE',NaN,'RejectionWilsonLow',NaN, ...
    'RejectionWilsonHigh',NaN,'MeanLegacyComponentDelta',NaN,'MeanIndicatorEPE',NaN);
end

function s=summarize(rows)
s=empty_summary();
s.N=rows(1).N; s.Jp=rows(1).Jp; s.Lambda=rows(1).Lambda; s.Beta=rows(1).Beta;
s.Attempted=numel(rows); good=[rows.Success]; s.Completed=sum(good);
s.Failed=s.Attempted-s.Completed;
if s.Completed==0, return; end
r=rows(good);
s.MeanPathEstimate=mean([r.PathEstimate]);
s.MeanPathBootstrapSE=mean([r.PathBootstrapSE]);
s.MeanPathCI025=mean([r.PathCI025]); s.MeanPathCI975=mean([r.PathCI975]);
s.MeanIndicatorDelta=mean([r.IndicatorDeltaMean]);
if s.Completed>1
    s.SDIndicatorDelta=std([r.IndicatorDeltaMean],0);
    s.MCSEIndicatorDelta=s.SDIndicatorDelta/sqrt(s.Completed);
end
s.MeanIndicatorLower95=mean([r.IndicatorLower95]);
s.MeanIndicatorP=mean([r.IndicatorP]);
s.RejectionRate=mean([r.IndicatorReject]);
s.RejectionMCSE=sqrt(s.RejectionRate*(1-s.RejectionRate)/s.Completed);
z=1.95996398454005; n=s.Completed; p=s.RejectionRate; denominator=1+z^2/n;
center=(p+z^2/(2*n))/denominator;
half=z*sqrt(p*(1-p)/n+z^2/(4*n^2))/denominator;
s.RejectionWilsonLow=max(0,center-half);
s.RejectionWilsonHigh=min(1,center+half);
s.MeanLegacyComponentDelta=mean([r.LegacyComponentDeltaMean]);
s.MeanIndicatorEPE=mean([r.IndicatorEPE]);
end

function R=assemble(cfg,rows,summary,parameterTables,elapsed)
R.Config=cfg; R.Summary=struct2table(summary); R.Replicates=struct2table(rows);
R.ParameterTables=parameterTables; R.ElapsedSeconds=elapsed;
R.Notes={ ...
    'Original GSCA_Reg and ALS_GSCA_Reg fitting and one-level bootstrap are used.', ...
    'Only in-bag auxiliary replacement and fixed raw target-indicator evaluation are changed in GSCA_Reg.', ...
    'Each delta is the original equal-bootstrap-replicate mean of ETC.delOPE_b_Boot; no row-reweighting or correction is added.', ...
    'TABLE.B columns 1:4 retain the original path estimate, bootstrap SE, and two-sided coefficient confidence interval.', ...
    'The predictive lower bound and p-value are copied directly from original TABLE.B columns 5:6.', ...
    'Rejection uses TABLE.B column5>0, retaining the original rounded percentile order-statistic convention; it is not replaced by another ECDF cutoff.', ...
    'Beta values are coefficient conditions, not calibrated predictive-null conditions. Report rejection rates without Type I error or power labels.', ...
    'The legacy component-loss delta is a different estimand and is retained only as a labeled descriptive benchmark.', ...
    'No failed dataset or bootstrap is redrawn. Summaries condition on successful datasets and include failure counts.', ...
    'Rejection MCSE and Wilson bounds describe simulation uncertainty; they are not predictive confidence intervals.'};
end

function save_outputs(R)
folder=R.Config.OutputDir;
save(fullfile(folder,'minimal_simulation.mat'),'R','-v7');
writetable(R.Summary,fullfile(folder,'summary.csv'));
writetable(R.Replicates,fullfile(folder,'replicates.csv'));
write_json(fullfile(folder,'summary.json'),table2struct(R.Summary));
write_json(fullfile(folder,'replicates.json'),table2struct(R.Replicates));
write_json(fullfile(folder,'config.json'),R.Config);
end

function write_json(path,value)
fid=fopen(path,'w','n','UTF-8');
if fid<0, error('gsca:MinimalOutput','Cannot write %s.',path); end
closer=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'%s\n',jsonencode(value));
end

