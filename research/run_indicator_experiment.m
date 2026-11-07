function Results = run_indicator_experiment(experiment,overrides)
%RUN_INDICATOR_EXPERIMENT Reproducible fixed-indicator comparison jobs.
% Run jobs separately to checkpoint and time bounded parts of the study.
%   run_indicator_experiment('smoke')
%   run_indicator_experiment('primary_n80_jp1')
%   run_indicator_experiment('primary_n80_jp3')
%   run_indicator_experiment('primary_n200_jp1')
%   run_indicator_experiment('primary_n200_jp3')
%   run_indicator_experiment('nested_n80_sample_boundary')
%   run_indicator_experiment('nested_n80_sample_positive')
%   run_indicator_experiment('nested_n80_sample_boundary_inner199')
%   run_indicator_experiment('nested_n80_oracle_boundary')
%   run_indicator_experiment('nested_n200_hetero_sample_boundary')
% overrides can change counts, Lambda sets, output location, etc.
% No job is automatically launched by opening this file.
if nargin<1 || isempty(experiment), experiment='smoke'; end
if nargin<2, overrides=struct; end
experiment=char(experiment);
if ~isstruct(overrides)||~isscalar(overrides)
    error('gsca:IndicatorExperiment','overrides must be a scalar structure.');
end
common=struct('OutputDir',fullfile(fileparts(mfilename('fullpath')), ...
    'results_indicator',experiment));
tokens=regexp(experiment,'^primary_n(80|200)_jp(1|3)$','tokens','once');
if ~isempty(tokens)
    common.Ns=str2double(tokens{1}); common.Jps=str2double(tokens{2});
    mode='primary';
else
    switch experiment
        case {'smoke','nested_smoke','primary','bridge','learned','pilot','sensitivity'}
            mode=experiment;
        case {'nested_n80_sample_boundary','nested_n80_sample_positive'}
            mode='nested';
            if endsWith(experiment,'boundary')
                common.Effects={'boundary'};
            else
                common.Effects={'positive'};
            end
        case 'nested_n80_sample_boundary_inner199'
            mode='nested'; common.Effects={'boundary'};
            common.Replications=100; common.InnerBoot=199;
            % Same seed/calibration/data as the first 100 primary nested
            % boundary datasets; only inner Monte Carlo precision changes.
        case 'nested_n80_oracle_boundary'
            mode='nested'; common.ReferenceModes={'oracle'};
            common.Effects={'boundary'}; common.Lambdas=0;
            common.Replications=100;
        case 'nested_n200_hetero_sample_boundary'
            mode='nested'; common.Ns=200; common.DGPs={'heteroscedastic'};
            common.Effects={'boundary'}; common.Lambdas=0;
            common.Replications=100;
        otherwise
            error('gsca:IndicatorExperiment','Unknown experiment: %s',experiment);
    end
end
fields=fieldnames(overrides);
for k=1:numel(fields), common.(fields{k})=overrides.(fields{k}); end
Results=Run_Simulation_Indicator_v4(mode,common);
write_json(fullfile(Results.Config.OutputDir,'summary.json'),table2struct(Results.Summary));
write_json(fullfile(Results.Config.OutputDir,'truth.json'),table2struct(Results.Truth));
write_json(fullfile(Results.Config.OutputDir,'config.json'),Results.Config);
fprintf('JSON_SUMMARY_BEGIN\n');
disp(jsonencode(table2struct(Results.Summary)));
fprintf('JSON_SUMMARY_END\n');
fprintf('JSON_TRUTH_BEGIN\n');
disp(jsonencode(table2struct(Results.Truth)));
fprintf('JSON_TRUTH_END\n');
end

function write_json(path,value)
fid=fopen(path,'w','n','UTF-8');
if fid<0, error('gsca:IndicatorOutput','Cannot write %s.',path); end
closer=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'%s\n',jsonencode(value));
end

