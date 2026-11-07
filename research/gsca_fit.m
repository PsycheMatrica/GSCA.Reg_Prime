function fit = gsca_fit(Data,W0,C0,B0,lambda,opts)
%GSCA_FIT Audited fixed-lambda GSCA solver for correction research.
% fit = gsca_fit(Data,W0,C0,B0,lambda,opts)
% W0 (J x P), C0 (P x J), B0 (P x P) are nonzero/free masks.
% opts: MaxIter=1000, Tol=1e-8, MinRcond=1e-10.
%
% The criterion is mean squared reconstruction error for modeled dependent
% indicators/components plus lambda*sum(B(:).^2). Lambda is on the
% mean-loss scale, as in ALS_GSCA_Reg; do not multiply it by sample size.
% Train-only population standard deviations are returned in fit.scale.
% Numerical failures return converged=false and a nonempty reason.
% Invalid inputs throw errors. Always check fit.converged before prediction.
%
% Research implementation: no claim that a covariance correction gives
% finite-sample unbiasedness or calibrated percentile inference.
if nargin < 6 || isempty(opts), opts=struct(); end
if ~isstruct(opts) || ~isscalar(opts)
    error('gsca_fit:Options','opts must be a scalar structure.');
end
opts=default_option(opts,'MaxIter',1000);
opts=default_option(opts,'Tol',1e-8);
opts=default_option(opts,'MinRcond',1e-10);
validateattributes(opts.MaxIter,{'numeric'},{'real','finite','scalar','integer','positive'});
validateattributes(opts.Tol,{'numeric'},{'real','finite','scalar','positive'});
validateattributes(opts.MinRcond,{'numeric'},{'real','finite','scalar','>',0,'<',1});
validateattributes(lambda,{'numeric'},{'real','finite','scalar','nonnegative'});
validateattributes(Data,{'numeric'},{'real','finite','2d','nonempty'});
[N,J]=size(Data);
if N<2, error('gsca_fit:SampleSize','At least two training rows are required.'); end
check_mask(W0,'W0'); check_mask(C0,'C0'); check_mask(B0,'B0');
P=size(W0,2);
if P<1 || size(W0,1)~=J || ~isequal(size(C0),[P,J]) || ~isequal(size(B0),[P,P])
    error('gsca_fit:Dimensions','Required sizes: W0 J-by-P, C0 P-by-J, B0 P-by-P.');
end
W0=W0~=0; C0=C0~=0; B0=B0~=0;
if any(sum(W0,1)==0)
    error('gsca_fit:EmptyComponent','Each component needs at least one indicator.');
end
dep=[any(C0,1),any(B0,1)];
if ~any(dep)
    error('gsca_fit:EmptyModel','Specify at least one dependent indicator or component.');
end
mu=mean(Data,1);
scale=std(Data,1,1);
W=double(W0); C=zeros(P,J); B=zeros(P,P);
fit=struct('mu',mu,'scale',scale,'W',W,'C',C,'B',B, ...
    'lambda',lambda,'converged',false,'iterations',0,'objective',NaN, ...
    'reason','','objectiveHistory',[],'minRcond',Inf);
if any(~isfinite(scale)) || any(scale<=0)
    fit.reason='zero_or_nonfinite_training_scale'; return
end
Z=(Data-mu)./scale;
S=(Z'*Z)/N;
S=(S+S')/2;
% A singular full S is permissible: only required subproblems are checked.
for p=1:P
    idx=find(W0(:,p));
    G=S(idx,idx);
    rc=rcond(G); fit.minRcond=min(fit.minRcond,rc);
    if ~isfinite(rc) || rc<opts.MinRcond
        fit.reason='ill_conditioned_weight_gram'; return
    end
    n2=W(idx,p)'*G*W(idx,p);
    if ~isfinite(n2) || n2<=0
        fit.reason='invalid_initial_component_norm'; return
    end
    W(idx,p)=W(idx,p)/sqrt(n2);
end
[C,B,ok,reason,minrc]=coefficient_step(S,W,C0,B0,lambda,opts.MinRcond);
fit.minRcond=min(fit.minRcond,minrc);
if ~ok
    fit.W=W; fit.C=C; fit.B=B; fit.reason=reason; return
end
old=criterion(S,W,C,B,dep,lambda);
fit.objectiveHistory=old;
if ~isfinite(old)
    fit.reason='nonfinite_initial_objective'; return
end
for iter=1:opts.MaxIter
    A=[C,B];
    for p=1:P
        idx=find(W0(:,p));
        % One-indicator weights have fixed magnitude under normalization.
        % The positive sign is an identification convention, not an update.
        if numel(idx)==1, continue; end
        % Refresh V after every coordinate update (including changed targets).
        V=[eye(J),W];
        V(:,J+p)=0;
        Wminus=W; Wminus(:,p)=0;
        Delta=Wminus*A-V;
        beta=-A(p,:); beta(J+p)=beta(J+p)+1;
        G=S(idx,idx);
        rc=rcond(G); fit.minRcond=min(fit.minRcond,rc);
        if ~isfinite(rc) || rc<opts.MinRcond
            fit.reason='ill_conditioned_weight_gram'; break
        end
        % Equivalent to the original Kronecker least-squares direction;
        % its positive scalar beta*beta' cancels under unit normalization.
        rhs=S(idx,:)*(Delta*beta');
        theta=G\rhs;
        n2=theta'*G*theta;
        if any(~isfinite(theta)) || ~isfinite(n2) || n2<=0
            fit.reason='zero_or_nonfinite_weight_update'; break
        end
        W(idx,p)=theta/sqrt(n2);
    end
    fit.iterations=iter;
    if ~isempty(fit.reason), break; end
    % Refit coefficients for the current W before computing the objective.
    [C,B,ok,reason,minrc]=coefficient_step(S,W,C0,B0,lambda,opts.MinRcond);
    fit.minRcond=min(fit.minRcond,minrc);
    if ~ok, fit.reason=reason; break; end
    current=criterion(S,W,C,B,dep,lambda);
    fit.objectiveHistory(end+1,1)=current;
    if ~isfinite(current)
        fit.reason='nonfinite_objective'; break
    end
    tol=opts.Tol*max(1,abs(old));
    change=old-current;
    if change < -tol
        fit.reason='objective_increase'; break
    end
    if abs(change)<=tol
        fit.converged=true; fit.reason='converged'; break
    end
    old=current;
end
if isempty(fit.reason), fit.reason='maximum_iterations'; end
% Identification uses the FIRST specified indicator in each component.
% These anchors are the same across all fits, without a fitted-data reference.
% Orientation is unstable when an anchor approaches zero; reject such fits.
% Joint transformation preserves every squared residual and the penalty.
d=ones(P,1);
for p=1:P
    idx=find(W0(:,p));
    raw=W(idx,p)./scale(idx)';
    if any(~isfinite(raw)) || abs(raw(1))<=1e-8*norm(raw)
        fit.W=W; fit.C=C; fit.B=B;
        fit.objective=criterion(S,W,C,B,dep,lambda);
        fit.converged=false; fit.reason='weak sign anchor'; return
    end
    if raw(1)<0, d(p)=-1; end
end
D=diag(d);
W=W*D; C=D*C; B=D*B*D;
fit.W=W; fit.C=C; fit.B=B;
fit.objective=criterion(S,W,C,B,dep,lambda);
if ~isfinite(fit.objective)
    fit.converged=false; fit.reason='nonfinite_final_objective';
end
end

function opts=default_option(opts,name,value)
if ~isfield(opts,name) || isempty(opts.(name)), opts.(name)=value; end
end

function check_mask(M,name)
if ~(isnumeric(M) || islogical(M)) || ~isreal(M) || ~ismatrix(M) || any(~isfinite(M(:)))
    error('gsca_fit:Mask','%s must be a finite real numeric or logical matrix.',name);
end
end

function [C,B,ok,reason,minrc]=coefficient_step(S,W,C0,B0,lambda,threshold)
[P,J]=size(C0);
C=zeros(P,J); B=zeros(P,P);
ok=true; reason=''; minrc=Inf;
for j=find(any(C0,1))
    idx=find(C0(:,j)); U=W(:,idx);
    G=U'*S*U; G=(G+G')/2;
    rc=rcond(G); minrc=min(minrc,rc);
    if ~isfinite(rc) || rc<threshold
        ok=false; reason='ill_conditioned_loading_gram'; return
    end
    C(idx,j)=G\(U'*S(:,j));
end
for q=find(any(B0,1))
    idx=find(B0(:,q)); U=W(:,idx);
    G=U'*S*U; G=(G+G')/2+lambda*eye(numel(idx));
    rc=rcond(G); minrc=min(minrc,rc);
    if ~isfinite(rc) || rc<threshold
        ok=false; reason='ill_conditioned_path_gram'; return
    end
    B(idx,q)=G\(U'*S*W(:,q));
end
if any(~isfinite(C(:))) || any(~isfinite(B(:)))
    ok=false; reason='nonfinite_coefficients';
end
end

function value=criterion(S,W,C,B,dep,lambda)
J=size(S,1);
R=[eye(J),W]-W*[C,B];
R=R(:,dep);
value=sum(sum(R.*(S*R)))+lambda*sum(B(:).^2);
end
