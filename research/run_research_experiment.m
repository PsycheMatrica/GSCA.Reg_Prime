function run_research_experiment(experiment)
% Reproducible bounded experiments used by the research CI workflow.
if nargin<1, experiment='primary'; end
common=struct('CalReplications',1000,'TruthReplications',2000,...
    'OutputDir',fullfile('research','results',experiment));
switch experiment
    case 'primary'
        common.CalReplications=100; common.TruthReplications=300;
        R=Run_Simulation_Correction_v3('ci',common);
    case 'sensitivity'
        common.Ns=[80,200]; common.Jps=3; common.Lambdas=0;
        common.DGPs={'normal','heteroscedastic'};
        common.ReferenceModes={'oracle','sample'};
        common.Replications=200; common.Boot=199;
        R=Run_Simulation_Correction_v3('ci',common);
    case {'nested_boundary','nested_positive'}
        common.Ns=80; common.Jps=3; common.Lambdas=0;
        common.Replications=100; common.Boot=199;
        common.Outer=99; common.InnerBoot=99;
        if strcmp(experiment,'nested_boundary')
            common.Effects={'boundary'};
        else
            common.Effects={'positive'};
        end
        R=Run_Simulation_Correction_v3('ci',common);
    otherwise
        error('gsca:experiment','Unknown experiment.');
end
fprintf('JSON_SUMMARY_BEGIN\n');
disp(jsonencode(table2struct(R.Summary)));
fprintf('JSON_SUMMARY_END\n');
fprintf('JSON_TRUTH_BEGIN\n');
disp(jsonencode(table2struct(R.Truth)));
fprintf('JSON_TRUTH_END\n');
end
