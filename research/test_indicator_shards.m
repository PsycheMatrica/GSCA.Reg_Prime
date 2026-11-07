function test_indicator_shards()
% Exact nested results must survive deterministic replication sharding.
% Persist a four-dataset full run and four one-dataset shards so the separate
% Python collector can independently test its pooled-summary recomputation.
outdir=fullfile(fileparts(mfilename('fullpath')), ...
    'results_indicator','shard_validation');
base=struct('Ns',80,'Jps',3,'Lambdas',0,'DGPs',{{'normal'}}, ...
    'ReferenceModes',{{'sample'}},'Effects',{{'boundary'}}, ...
    'Replications',4,'TotalReplications',4,'ReplicationStart',1, ...
    'Boot',49,'Outer',9,'InnerBoot',49, ...
    'CalReplications',12,'TruthReplications',12,'CalIterations',4, ...
    'ProgressEvery',0,'Seed',67019,'OutputDir',fullfile(outdir,'full'));
full=run_indicator_experiment('nested_n80_sample_boundary',base);
parts=cell(4,1);
for k=1:4
    opt=base; opt.Replications=1; opt.ReplicationStart=k;
    opt.OutputDir=fullfile(outdir,sprintf('part%d',k));
    parts{k}=run_indicator_experiment(sprintf('nested_n80_sample_boundary_part%d',k),opt);
end
joined=[parts{1}.Replicates;parts{2}.Replicates;parts{3}.Replicates;parts{4}.Replicates];
assert(isequaln(full.Replicates,joined), ...
    'Sharding changed a dataset, seed, fit result, percentile, or diagnostic.');
assert(isequal(full.Replicates.Rep,(1:4)'), ...
    'Full run did not retain global replication identifiers.');
assert(all(full.Replicates.Valid) && all(full.Replicates.OuterValid), ...
    'A deterministic shard test fit failed.');
assert(full.Config.TotalReplications==4 && full.Summary.Attempted==4);
for k=1:4
    part=parts{k};
    assert(part.Replicates.Rep==k && part.Config.TotalReplications==4 && ...
        part.Summary.Attempted==1,'Shard identifiers or denominators are incorrect.');
    assert(isequaln(full.Truth,part.Truth),'Independent truth estimates changed across shards.');
    assert(isequaln(full.TruthSamples,part.TruthSamples), ...
        'Raw independent training truth draws changed across shards.');
    assert(isequaln(full.Calibration,part.Calibration), ...
        'Boundary calibration changed across shards.');
    for name={'summary.json','truth.json','config.json','replicates.csv'}
        assert(exist(fullfile(outdir,sprintf('part%d',k),name{1}),'file')==2, ...
            'A shard fixture artifact is missing.');
    end
end
bad=struct('Replications',2,'ReplicationStart',4,'TotalReplications',4);
caught=false;
try
    Run_Simulation_Indicator_v4('smoke',bad);
catch err
    caught=strcmp(err.identifier,'gsca:InvalidReplicationRange');
end
assert(caught,'Out-of-range shard indices were not rejected.');
fprintf('test_indicator_shards: PASS (4 datasets == 1 + 1 + 1 + 1; nested draws, truth and calibration identical)\n');
fprintf('SHARD_FIXTURE_DIR=%s\n',outdir);
end

