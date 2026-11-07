function example = run_indicator_example(outputDirectory)
%RUN_INDICATOR_EXAMPLE Illustrate fixed-indicator GSCA correction on ACSI replica.
%   example = run_indicator_example()
%   example = run_indicator_example(outputDirectory)
%
%   Uses the EXISTING ACSI_Comp_Replica.csv distributed with this repository.
%   This is an illustration with a replica, not new empirical observations.
%   Model masks match Run_Example_GSCA_Reg.m. The penalty is prespecified at
%   Lambda=0; there is no CV selection and no outer bootstrap.
%
%   Outputs are point estimates and DIAGNOSTIC one-sided percentile results.
%   pRaw/pShifted do not establish calibrated finite-N predictive inference.
%   A failed fit/draw is recorded, not dropped or redrawn. Files are saved
%   even if out.valid=false. The return value exposes failure status so the
%   caller can inspect it before deciding whether to assert in validation.
researchDirectory = fileparts(mfilename('fullpath'));
repositoryDirectory = fileparts(researchDirectory);
if nargin<1 || isempty(outputDirectory)
    outputDirectory = fullfile(researchDirectory,'results_indicator','acsi_example');
end
if ~exist(outputDirectory,'dir'), mkdir(outputDirectory); end
dataPath = fullfile(repositoryDirectory,'ACSI_Comp_Replica.csv');
inputTable = readtable(dataPath);
Data = table2array(inputTable);
validateattributes(Data,{'numeric'},{'2d','real','finite','nonempty'});
if size(Data,2)~=14
    error('gsca:indicatorExampleInput','The existing ACSI replica must have 14 numeric columns.');
end

W0 = [1 1 1 0 0 0 0 0 0 0 0 0 0 0; ...
      0 0 0 1 1 1 0 0 0 0 0 0 0 0; ...
      0 0 0 0 0 0 1 1 0 0 0 0 0 0; ...
      0 0 0 0 0 0 0 0 1 1 1 0 0 0; ...
      0 0 0 0 0 0 0 0 0 0 0 1 0 0; ...
      0 0 0 0 0 0 0 0 0 0 0 0 1 1]';
C0 = W0';
B0 = [0 1 1 1 0 0; ...
      0 0 1 1 0 0; ...
      0 0 0 1 0 0; ...
      0 0 0 0 1 1; ...
      0 0 0 0 0 1; ...
      0 0 0 0 0 0];
opts = struct('Lambda',0,'Boot',199,'Outer',0,'Seed',2026110609, ...
    'Alpha',.05,'MaxIter',1000,'Tol',1e-8,'MinRcond',1e-10);
metadata = struct( ...
    'dataDescription','Existing distributed ACSI replica; not new empirical data', ...
    'dataFile','ACSI_Comp_Replica.csv', ...
    'modelSource','Run_Example_GSCA_Reg.m', ...
    'N',size(Data,1),'J',size(Data,2),'P',size(B0,1), ...
    'loss','Mean squared prediction error across fixed target indicators in raw units', ...
    'projection','Full-data plug-in moments for other-parent score projection', ...
    'regularization','Prespecified Lambda=0; no cross-validation or tuning', ...
    'inference','Point estimates and diagnostic one-sided percentiles only', ...
    'exceptionIdentifier','','exceptionMessage','');
try
    out = gsca_indicator_predictive(Data,W0,C0,B0,opts);
catch err
    % Preserve unexpected failures for inspection; do not fabricate estimates.
    out = struct('valid',false,'reason',['Exception: ',err.message], ...
        'failures',NaN,'outerValid',false);
    metadata.exceptionIdentifier = err.identifier;
    metadata.exceptionMessage = err.message;
end

[source,target] = find(B0);
E = numel(source);
IndicatorCount = sum(C0(target,:)~=0,2);
deltaRaw = field_column(out,'deltaRaw',E);
correction = field_column(out,'correction',E);
deltaVC = field_column(out,'deltaCorrected',E);
pRaw = field_column(out,'pRaw',E);
pShifted = field_column(out,'pShifted',E);
lowerRaw95 = field_column(out,'lowerRaw',E);
lowerShifted95 = field_column(out,'lowerShifted',E);

N = size(Data,1);
counts = nan(N,1);
if isfield(out,'counts'), counts = out.counts; end
minCount = NaN;
if all(isfinite(counts)), minCount = min(counts); end
drawSucceeded = false(opts.Boot,1);
drawStatus = repmat({'not attempted'},opts.Boot,1);
attempted = 0;
if isfield(out,'draws')
    drawSucceeded = all(isfinite(out.draws.raw),2);
    attempted = opts.Boot;
end
if isfield(out,'failureReasons')
    for b=1:opts.Boot
        if drawSucceeded(b)
            drawStatus{b} = 'ok';
        elseif ~isempty(out.failureReasons{b})
            drawStatus{b} = out.failureReasons{b};
        else
            drawStatus{b} = 'invalid draw with no reported reason';
        end
    end
end
if ~isempty(metadata.exceptionMessage)
    % An exception can interrupt a call before it returns progress records.
    attempted = NaN;
    drawStatus(:) = {'unknown: analysis raised an exception'};
end
successfulDraws = repmat(sum(drawSucceeded),E,1);
attemptedDraws = repmat(attempted,E,1);
failedDraws = repmat(out.failures,E,1);
minOOBCount = repmat(minCount,E,1);
valid = repmat(out.valid,E,1);
pathTable = table(source,target,IndicatorCount,deltaRaw,correction,deltaVC, ...
    pRaw,pShifted,lowerRaw95,lowerShifted95,attemptedDraws,successfulDraws, ...
    failedDraws,minOOBCount,valid);
observation = (1:N)';
countTable = table(observation,counts,'VariableNames',{'observation','OOBCount'});
bootstrap = (1:opts.Boot)';
failureTable = table(bootstrap,drawSucceeded,drawStatus, ...
    'VariableNames',{'bootstrap','succeeded','status'});

paths = struct( ...
    'pathResults',fullfile(outputDirectory,'acsi_path_results.csv'), ...
    'oobCounts',fullfile(outputDirectory,'acsi_oob_counts.csv'), ...
    'drawDiagnostics',fullfile(outputDirectory,'acsi_draw_diagnostics.csv'), ...
    'workspace',fullfile(outputDirectory,'acsi_example.mat'), ...
    'notes',fullfile(outputDirectory,'acsi_example_notes.txt'));
writetable(pathTable,paths.pathResults);
writetable(countTable,paths.oobCounts);
writetable(failureTable,paths.drawDiagnostics);
save(paths.workspace,'out','opts','metadata','pathTable','countTable','failureTable','W0','C0','B0');

fid = fopen(paths.notes,'w');
if fid<0, error('gsca:indicatorExampleWrite','Cannot write example notes.'); end
closeNotes = onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'ACSI replica fixed-indicator correction illustration\n');
fprintf(fid,'Data: %s\n',metadata.dataDescription);
fprintf(fid,'N=%d; J=%d; P=%d; structural paths=%d\n',metadata.N,metadata.J,metadata.P,E);
fprintf(fid,'Model masks: %s\n',metadata.modelSource);
fprintf(fid,'Loss: %s\n',metadata.loss);
fprintf(fid,'Projection: %s\n',metadata.projection);
fprintf(fid,'Regularization: %s\n',metadata.regularization);
fprintf(fid,'Boot=%d; Outer=%d; Seed=%d; Alpha=%.2f\n',opts.Boot,opts.Outer,opts.Seed,opts.Alpha);
fprintf(fid,'Valid=%d; status=%s\n',out.valid,out.reason);
fprintf(fid,'Attempted draws=%g; successful draws=%d; failed draws=%g; minimum OOB count=%g\n', ...
    attempted,sum(drawSucceeded),out.failures,minCount);
fprintf(fid,'The correction is a variance/bagging candidate, not a proven finite-N bias correction.\n');
fprintf(fid,'pRaw and pShifted invert diagnostic one-sided replicate percentiles.\n');
fprintf(fid,'They are not validated predictive p-values. No outer-bootstrap inference was run.\n');
fprintf(fid,'Raw indicator units are retained; block risks average squared errors, not residuals.\n');
fprintf(fid,'Failed draws are not dropped or redrawn; invalid analyses retain NaN results.\n');
if ~isempty(metadata.exceptionMessage)
    fprintf(fid,'Exception identifier: %s\nException message: %s\n', ...
        metadata.exceptionIdentifier,metadata.exceptionMessage);
end
example = struct('valid',out.valid,'reason',out.reason,'results',pathTable, ...
    'diagnostics',out,'metadata',metadata,'paths',paths);
fprintf('ACSI REPLICA illustration: valid=%d; %s\n',out.valid,out.reason);
fprintf('Fixed Lambda=0; Boot=%d; Outer=0; point estimates and diagnostic percentiles only.\n',opts.Boot);
disp(pathTable);
fprintf('JSON_ACSI_BEGIN\n');
disp(jsonencode(table2struct(pathTable)));
fprintf('JSON_ACSI_END\n');
fprintf('JSON_ACSI_METADATA_BEGIN\n');
disp(jsonencode(metadata));
fprintf('JSON_ACSI_METADATA_END\n');
if ~out.valid
    warning('gsca:indicatorExampleInvalid', ...
        'Illustration invalid: %s. Failure details and NaN results were saved.',out.reason);
end
end

function value = field_column(out,name,E)
value = nan(E,1);
if isfield(out,name) && numel(out.(name))==E
    value = reshape(out.(name),E,1);
end
end

