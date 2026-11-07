function test_gsca_original_flow()
%TEST_GSCA_ORIGINAL_FLOW Preserve original GSCA fitting and bootstrap stream.
% The fixture is the main-branch GSCA_Reg.m at commit
% 9d586446970a47ebd420516a0451143ba608fff9 with ONLY its first function renamed.
% Both functions call the same original, unchanged ALS_GSCA_Reg implementation.
% Predictive evaluation changes are tested separately from fitted parameters.
previous=rng; cleanup=onCleanup(@()rng(previous)); %#ok<NASGU>
rng(3817,'twister');
N=120; P=3; Jp=3;
g1=randn(N,1);
g2=.4*g1+sqrt(.84)*randn(N,1);
g3=.25*g1+.6*g2+sqrt(.4575)*randn(N,1);
G=[g1,g2,g3];
loadings=[.7,.8,.9];
Data=zeros(N,P*Jp);
for p=1:P
    cols=(p-1)*Jp+(1:Jp);
    Data(:,cols)=G(:,p)*loadings+randn(N,Jp).*sqrt(1-loadings.^2);
end
% Unequal scales/nonzero means make raw-unit outcome changes observable.
Data=Data.*[1,1.2,.8,1.1,.9,1.3,2,.7,1.5]+[2,-1,3,0,4,-2,5,-3,1];
W0=kron(eye(P),ones(Jp,1)); C0=W0';
B0=zeros(P); B0(1:2,3)=1;
lambda=[0,.1]; KxK=[2,3]; Nboot=100;
seed=692031;
rng(seed,'twister');
original=GSCA_Reg_OriginalFixture(Data,W0,C0,B0,lambda,KxK,Nboot,1000,1e-8,false);
oldRng=rng;
rng(seed,'twister');
updated=GSCA_Reg(Data,W0,C0,B0,lambda,KxK,Nboot,1000,1e-8,false);
newRng=rng;
assert(isequal(oldRng,newRng),'Predictive evaluation changed the RNG stream.');

% All preexisting INI fields, including selected lambda and legacy GoF,
% must remain bit-for-bit identical. The only added INI field is separate.
f=fieldnames(original.INI);
for k=1:numel(f)
    assert(isequaln(original.INI.(f{k}),updated.INI.(f{k})), ...
        ['Original INI field changed: ',f{k}]);
end
for field={'W_Boot','C_Boot','B_Boot'}
    name=field{1};
    assert(isequaln(original.ETC.(name),updated.ETC.(name)), ...
        ['Bootstrap estimates changed: ',name]);
    assert(all(isfinite(updated.ETC.(name)(:))),['Nonfinite parameter draws: ',name]);
end
for field={'W','C','B'}
    name=field{1};
    assert(isequaln(original.TABLE.(name)(:,1:4),updated.TABLE.(name)(:,1:4)), ...
        ['Parameter estimate/SE/CI changed: ',name]);
end

% Stored legacy predictive draws reproduce the original predictive columns.
% The new method retains the ORIGINAL rank convention and one-sided bound.
rank=round(.05*(Nboot-1))+1;
assert(isequal(size(updated.ETC.delOPE_b_Boot),[nnz(B0),Nboot]));
assert(isequal(size(updated.ETC.delOPE_c_Boot),[nnz(C0),Nboot]));
check_columns(original.TABLE.B,updated.ETC.delOPE_b_legacy_Boot,rank,'legacy B');
check_columns(original.TABLE.C,updated.ETC.delOPE_c_legacy_Boot,rank,'legacy C');
check_columns(updated.TABLE.B,updated.ETC.delOPE_b_Boot,rank,'updated B');
check_columns(updated.TABLE.C,updated.ETC.delOPE_c_Boot,rank,'updated C');
assert(all(isfinite(updated.ETC.delOPE_b_Boot(:))),'Nonfinite new B losses.');
assert(all(isfinite(updated.ETC.delOPE_c_Boot(:))),'Nonfinite new C losses.');
assert(norm(updated.ETC.delOPE_b_Boot-updated.ETC.delOPE_b_legacy_Boot,'fro')>1e-6, ...
    'New indicator-block/replacement evaluation unexpectedly equals legacy score deletion.');

assert(isequal(size(updated.ETC.indicatorOPE_Boot),[P,Nboot]));
expected=mean(updated.ETC.indicatorOPE_Boot,2)';
assert(isequaln(updated.INI.indicatorOPE_byComponent,expected), ...
    'Indicator OPE summary is not the equal-repetition mean.');
assert(all(isnan(updated.INI.indicatorOPE_byComponent(1:2))), ...
    'Exogenous components should have no structural indicator prediction error.');
assert(isfinite(updated.INI.indicatorOPE_byComponent(3)), ...
    'Endogenous target indicator prediction error is missing.');

% Preserve the original no-bootstrap behavior and initialize new fields.
rng(seed,'twister');
old0=GSCA_Reg_OriginalFixture(Data,W0,C0,B0,0,KxK,0,1000,1e-8,false);
rng(seed,'twister');
new0=GSCA_Reg(Data,W0,C0,B0,0,KxK,0,1000,1e-8,false);
f=fieldnames(old0.INI);
for k=1:numel(f)
    assert(isequaln(old0.INI.(f{k}),new0.INI.(f{k})), ...
        ['No-bootstrap INI changed: ',f{k}]);
end
for field={'W','C','B'}
    name=field{1};
    assert(isequaln(old0.TABLE.(name),new0.TABLE.(name)), ...
        ['No-bootstrap parameter table changed: ',name]);
end
assert(isempty(new0.ETC.indicatorOPE_Boot) && isempty(new0.ETC.delOPE_b_Boot) && ...
    isempty(new0.ETC.delOPE_b_legacy_Boot) && isempty(new0.ETC.delOPE_c_Boot) && ...
    isempty(new0.ETC.delOPE_c_legacy_Boot),'No-bootstrap diagnostic arrays must be empty.');
assert(all(isnan(new0.INI.indicatorOPE_byComponent)), ...
    'No-bootstrap indicator OPE must be unavailable.');

fprintf(['test_gsca_original_flow: exact original parameters, CV selection, GoF, ' ...
    'bootstrap draws, RNG, parameter summaries, and percentile convention preserved.\n']);
end

function check_columns(table,draws,rank,label)
sorted=sort(draws,2);
assert(isequaln(table(:,5),sorted(:,rank)),['Lower-bound convention changed: ',label]);
assert(isequaln(table(:,6),mean(draws<=0,2)),['One-sided p-value changed: ',label]);
end

