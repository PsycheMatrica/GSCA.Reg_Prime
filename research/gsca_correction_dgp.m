function [Data,Model] = gsca_correction_dgp(N,Jp,b,kind,seed)
%GSCA_CORRECTION_DGP Exact-composite DGP for correction experiments.
% P=3; gamma_1 and gamma_2 correlate .30; gamma_3 has paths b and .35.
% Within each Jp-indicator block, X_p*w_p=gamma_p exactly. Indicator
% residuals lie in the orthogonal complement of w_p. This is a component
% DGP, not a reflective-factor approximation with noisy component scores.
% Normal and heteroscedastic conditions have the same population covariance.
if nargin<4 || isempty(kind), kind='normal'; end
if nargin<5, seed=20261106; end
validateattributes(N,{'numeric'},{'scalar','integer','>=',4});
validateattributes(Jp,{'numeric'},{'scalar','integer','>=',1});
validateattributes(b,{'numeric'},{'scalar','real','finite'});
kind=char(kind);
if ~any(strcmp(kind,{'normal','heteroscedastic'}))
    error('gsca:InvalidDGP','kind must be normal or heteroscedastic.');
end
old=rng; restore=onCleanup(@()rng(old)); %#ok<NASGU>
rng(seed,'twister');
P=3; J=P*Jp; rho=.30; nuisance=.35; within=.50;
sigma2=1-(b^2+nuisance^2+2*rho*b*nuisance);
if sigma2<=0
    error('gsca:InvalidDGP','The structural residual variance must be positive.');
end
u=randn(N,3);
Gamma=zeros(N,P);
Gamma(:,1)=u(:,1);
Gamma(:,2)=rho*u(:,1)+sqrt(1-rho^2)*u(:,2);
innovation=u(:,3);
if strcmp(kind,'heteroscedastic')
    innovation=innovation.*sqrt(.5+.5*Gamma(:,1).^2);
end
Gamma(:,3)=b*Gamma(:,1)+nuisance*Gamma(:,2)+sqrt(sigma2)*innovation;
covGamma=[1 rho b+rho*nuisance; ...
          rho 1 rho*b+nuisance; ...
          b+rho*nuisance rho*b+nuisance 1];
blockCov=(1-within)*eye(Jp)+within*ones(Jp);
w=ones(Jp,1)/sqrt(Jp*(1+(Jp-1)*within));
c=blockCov*w;
contrast=eye(Jp)-ones(Jp)/Jp;
Data=zeros(N,J);
for p=1:P
    ix=(p-1)*Jp+(1:Jp);
    residual=sqrt(1-within)*randn(N,Jp)*contrast;
    Data(:,ix)=Gamma(:,p)*c'+residual;
end
covData=kron(covGamma,c*c')+kron(eye(P),(1-within)*contrast);
Model.W0=kron(eye(P),ones(Jp,1));
Model.C0=Model.W0';
Model.B0=zeros(P); Model.B0(1:2,3)=1;
Model.Wtrue=kron(eye(P),w);
Model.Ctrue=kron(eye(P),c');
Model.Btrue=zeros(P); Model.Btrue(1:2,3)=[b;nuisance];
Model.mean=zeros(1,J);
Model.covariance=(covData+covData')/2;
Model.componentCovariance=covGamma;
Model.sigma2=sigma2;
Model.focalIndex=nnz(Model.C0)+find(find(Model.B0)==sub2ind([P,P],1,3));
Model.kind=kind;
Model.N=N; Model.P=P; Model.J=J; Model.Jp=Jp;
Model.beta=b;
Model.compositeIdentityMaxError=max(abs(Data*Model.Wtrue-Gamma),[],'all');
end
