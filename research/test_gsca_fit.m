function results=test_gsca_fit()
%TEST_GSCA_FIT Deterministic numerical checks for the research solver.
% Run from the repository after addpath('research').
previous=rng; cleaner=onCleanup(@() rng(previous)); %#ok<NASGU>
rng(24681357,'twister');
names={};
N=160;
x1=randn(N,1); x2=.4*x1+sqrt(.84)*randn(N,1);
y=.35*x1-.2*x2+.8*randn(N,1);
Data=[x1,x2,y];
W0=eye(3); C0=zeros(3); B0=zeros(3); B0(1:2,3)=1;
for lambda=[0,.2]
    f=gsca_fit(Data,W0,C0,B0,lambda);
    assert(f.converged && strcmp(f.reason,'converged'),'Singleton fit failed.');
    Z=(Data-f.mu)./f.scale;
    expected=(Z(:,1:2)'*Z(:,1:2)/N+lambda*eye(2))\(Z(:,1:2)'*Z(:,3)/N);
    close(f.W,eye(3),1e-12,'singleton weights');
    close(f.B(1:2,3),expected,1e-11,'OLS/ridge coefficients');
    direct=mean((Z(:,3)-Z(:,1:2)*expected).^2)+lambda*sum(expected.^2);
    close(f.objective,direct,1e-12,'OLS/ridge objective');
    names{end+1}=sprintf('singleton_lambda_%g',lambda); %#ok<AGROW>
end
% Reordering moves the dependent component before its predictors.
p=[3,1,2];
g=gsca_fit(Data,W0(:,p),C0(p,:),B0(p,p),.2);
assert(g.converged,'Reordered component fit failed.');
close(g.W,f.W(:,p),1e-12,'reordered weights');
close(g.B,f.B(p,p),1e-11,'reordered paths');
close(g.objective,f.objective,1e-12,'reordered objective');
names{end+1}='component_reordering';

% Measurement-only models previously used an undefined penalty variable.
u=randn(N,1);
M=[.6*u+.8*randn(N,1), .8*u+.6*randn(N,1), .9*u+sqrt(.19)*randn(N,1)];
fm=gsca_fit(M,ones(3,1),ones(1,3),0,0);
assert(fm.converged && strcmp(fm.reason,'converged'),'Measurement-only fit failed.');
assert(fm.W(1)>0,'First specified indicator must have positive weight.');
names{end+1}='measurement_only_and_sign';

% The saved objective must use the returned final W, C, B.
Zm=(M-fm.mu)./fm.scale;
Gm=Zm*fm.W;
direct=mean(sum((Zm-Gm*fm.C).^2,2));
close(fm.objective,direct,1e-11,'final reconstruction objective');
close(fm.objectiveHistory(end),fm.objective,1e-11,'final history objective');
h=fm.objectiveHistory(:);
assert(all(diff(h)<=1e-8*max(1,abs(h(1:end-1)))+1e-12), ...
    'Accepted objective history contains a substantive increase.');
assert(fm.iterations>=1 && fm.iterations<=1000 && fm.minRcond>0, ...
    'Convergence diagnostics are missing or invalid.');
names{end+1}='final_objective_and_diagnostics';

% Exhausting the requested iteration budget is a reported failure.
limited=gsca_fit(M,ones(3,1),ones(1,3),0,0,struct('MaxIter',1,'Tol',1e-16));
assert(~limited.converged && strcmp(limited.reason,'maximum_iterations'), ...
    'Iteration exhaustion must not be marked converged.');
names{end+1}='maximum_iterations';

constant=gsca_fit([M(:,1),ones(N,1)],ones(2,1),ones(1,2),0,0);
assert(~constant.converged && strcmp(constant.reason,'zero_or_nonfinite_training_scale'), ...
    'Constant training indicators must be rejected.');
names{end+1}='constant_indicator_failure';

singular=gsca_fit([M(:,1),M(:,1)],ones(2,1),ones(1,2),0,0);
assert(~singular.converged && strcmp(singular.reason,'ill_conditioned_weight_gram'), ...
    'Singular required weight Gram matrix must be rejected.');
names{end+1}='singular_gram_failure';

% Exact balanced design: the first weight update has zero anchor.
H=[1,1,1;1,-1,-1;-1,1,-1;-1,-1,1];
S=[1,-.5,-.5;-.5,1,.8;-.5,.8,1];
weak=gsca_fit(H*chol(S),ones(3,1),ones(1,3),0,0,struct('MaxIter',1));
assert(~weak.converged && strcmp(weak.reason,'weak sign anchor'), ...
    'A zero first-indicator anchor must be rejected.');
names{end+1}='weak_fixed_anchor_failure';

threw=false;
try
    gsca_fit(Data,W0,C0,B0,[0,.1]);
catch
    threw=true;
end
assert(threw,'A lambda vector must not silently select or broadcast tuning values.');
names{end+1}='scalar_lambda_validation';

results=struct('passed',numel(names),'checks',{names});
fprintf('test_gsca_fit: %d deterministic checks passed.\n',numel(names));
end

function close(actual,expected,tolerance,label)
assert(isequal(size(actual),size(expected)),['Shape mismatch: ',label]);
assert(all(isfinite(actual(:))) && ...
    norm(actual(:)-expected(:),Inf)<=tolerance,['Numerical mismatch: ',label]);
end
