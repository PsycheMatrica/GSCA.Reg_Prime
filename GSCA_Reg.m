function Results=GSCA_Reg(Data,W,C,B,vecLambda,KxK,N_Boot,Max_iter,Min_limit,Flag_Parallel)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RegGSCA() - MATLAB function to perform Regularized Generalized          %
%               Structured Component Analysis (GSCA).                     %
% Author: Gyeongcheol Cho                                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input arguments:                                                        %
%   Data = an N by J matrix of scores for N individuals on J indicators   %
%   W = a J by P matrix of weight parameters                              %
%   C = a P by J matrix of loading parameters                             %
%   B = a P by P matrix of path coefficients                              %
%   vecLambda = a vector of lambda values for regularization              %
%   KxK = a 1 by 2 vector of the number of folds and the number of        %
%        replications for cross-validation                                %
%   N_Boot = Integer representing the number of bootstrap samples for     %
%            calculating standard errors (SE) and 95% confidence          %
%            intervals (CI).                                              %
%   Max_iter = Maximum number of iterations for the Alternating Least     % 
%              Squares (ALS) algorithm                                    %
%   Min_limit = Tolerance level for ALS algorithm                         %
%   Flag_Parallel = Logical value to determine whether to use parallel    %
%                   computing for bootstrapping                           %
% Output arguments:                                                       %
%   Results: Structure array containing (1) results from the original     %
%       sample (INI); (2) summary tables with standard errors and         %
%       confidence intervals (TABLE); and (3) bootstrap estimates for     %
%       various parameter sets (ETC).                                     %    
%   .INI: Strucutre array containing goodness-of-fit values, R-squared    % 
%        values, and matrices parameter estimates                         %
%     .GoF = [FIT_D,   OPE_D;                                             %
%             FIT_M_D, OPE_M_D;                                           %
%             FIT_S_D, OPE_S_D];                                          %
%     .R2m = Vector of R-squared values for dependent variables           %
%               in the measurement model                                  %
%     .R2s = Vector of R-squared values for dependent variables           %
%               in the structural model                                   %
%     .W: a J by P matrix of weight estimates                             %
%     .C: a P by J matrix of loading estimates                            %
%     .B: a P by P matrix of path coefficient estimates                   %
%     .lambda_optimal: Optimal lambda value                               %
%     .minVE: minimum value of the validation error                       %
%     .LamVE: a matrix of lambda values and their corresponding           %
%             validation errors                                           %
%  .TABLE: Structure array containing tables of parameter estimates, their%
%         SEs, 95% CIs,and other statistics                               %
%     .W: Table for weight estimates                                      %
%     .C: Table for loading estimates                                     %
%     .B: Table for path coefficients estimates                           %
%  .ETC: Structure array including bootstrapped parameter estmates        %
%     .W_Boot: Matrix of bootstrapped weight estimates                    %
%     .C_Boot: Matrix of bootstrapped loading estimates                   %
%     .B_Boot: Matrix of bootstrapped path coefficient estimates          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
if nargin<10; Flag_Parallel=false;
    if nargin<9; Min_limit=10^(-8); 
        if nargin<8; Max_iter=1000; 
            if nargin<7; N_Boot=0; 
                if nargin<6; KxK=[1,5]; 
                    if nargin<5; vecLambda=0;
                    end
                end
            end
        end
    end
end

%% (1) Preliminary stage
    Z=Data;  
    [N,J]=size(Z);
    P=size(B,1);
    T=J+P;
                       
    % index
    W0=W~=0; Nw=sum(sum(W0,1),2);
    C0=C~=0; Nc=sum(sum(C0,1),2);
    B0=B~=0; Nb=sum(sum(B0,1),2);
    B0=B~=0;
    ind_Cdep=sum(C0,1)>0; Jy = sum(ind_Cdep,2); loc_Cdep=find(ind_Cdep);
    ind_Bdep=sum(B0,1)>0; Py = sum(ind_Bdep,2); loc_Bdep=find(ind_Bdep); Ty = Jy+Py;
    ind_Adep=[ind_Cdep, ind_Bdep];
    loc_w_t=cell(1,P);
    loc_b_t=cell(1,P);
    for p=1:P
        loc_w_t{1,p}=find(W(:,p))';
        loc_b_t{1,p}=find(B(:,p))';
    end
    loc_c_t=cell(1,J);
    for j=1:J
        loc_c_t{1,j}=find(C(:,j))';
    end
 %       ind_exo=find(structural_eq);
    % Setting the intital values for A,W    
        W(W0)=1;
%% (2) Identify the optimal lambda value
    Nfold=KxK(1,1);
    K=KxK(1,2);
    N_lambda=size(vecLambda,2);
    MatVE=zeros(N_lambda,K,Nfold);
    if N_lambda==1
        lambda_optimal=vecLambda;
        minVE=NaN;
        LamVE=[NaN,NaN];
    else
        Zf=Z;
        for nfold=1:Nfold
            if nfold>1; Zf=Z(randperm(N),:); end
            N_tr=floor(N/K);
            vecST=1:N_tr:N;
            vecST=vecST(1,1:K);
            vecED=[vecST(1,2:K)-1,N];
            for k=1:K 
                id_tr=true(1,N);
                loc_tt=vecST(1,k):vecED(1,k);
                id_tr(1,loc_tt)=false;
                N_tr=sum(id_tr,2);
                N_tt=size(loc_tt,2);

                Zf_tp=Zf;  % for parfor
                Z_tr=Zf_tp(id_tr,:);
                Z_tt=Zf_tp(loc_tt,:);
                mean_Z_tr=mean(Z_tr);
                std_Z_tr=std(Z_tr,1); 
                Z_tt=(Z_tt-ones(N_tt,1)*mean_Z_tr)./(ones(N_tt,1)*std_Z_tr);
                if Flag_Parallel
                    parfor n_l=1:N_lambda
                        vecLambda_tp=vecLambda;
                        lambda_b=vecLambda_tp(1,n_l);
                        [W_tr,C_tr,B_tr,~]=ALS_GSCA_Reg(Z_tr,W,W0,C0,B0,ind_Adep,lambda_b,Min_limit,Max_iter,N_tr,J,P,T,Jy,Py,loc_Cdep,loc_Bdep);
                        CV_tt=Z_tt*W_tr;
                        e_tt=[Z_tt,CV_tt]-CV_tt*[C_tr,B_tr];
                        vecVE=sum(e_tt(:,ind_Adep).^2,1)/N_tt;
                        MatVE(n_l,k,nfold)=mean(vecVE,2);
                    end
                else
                    for n_l=1:N_lambda
                        vecLambda_tp=vecLambda;
                        lambda_b=vecLambda_tp(1,n_l);
                        [W_tr,C_tr,B_tr,~]=ALS_GSCA_Reg(Z_tr,W,W0,C0,B0,ind_Adep,lambda_b,Min_limit,Max_iter,N_tr,J,P,T,Jy,Py,loc_Cdep,loc_Bdep);
                        CV_tt=Z_tt*W_tr;
                        e_tt=[Z_tt,CV_tt]-CV_tt*[C_tr,B_tr];
                        vecVE=sum(e_tt(:,ind_Adep).^2,1)/N_tt;
                        MatVE(n_l,k,nfold)=mean(vecVE,2);
                    end                    
                end
            end
        end
        LamVE=[vecLambda',mean(mean(MatVE,3),2)];
        %plot(LamVE)
        [minVE,locVE]=min(LamVE(:,2));
        lambda_optimal=LamVE(locVE,1);
    end
    INI.minVE=minVE;
    INI.LamVE=LamVE;
    INI.lambda_optimal=lambda_optimal;
%% (3) Estimation of paramters
    [est_W,est_C,est_B,vec_err]=ALS_GSCA_Reg(Z,W,W0,C0,B0,ind_Adep,lambda_optimal,Min_limit,Max_iter,N,J,P,T,Jy,Py,loc_Cdep,loc_Bdep);            
    R_squared_dep=ones(1,Ty)-vec_err;
    R_squared=zeros(1,T);
    R_squared(1,ind_Adep)=R_squared_dep;
    Eval=[mean(R_squared(1,[ind_Cdep,ind_Bdep])),NaN; %% FIT_D
          mean(R_squared(1,[ind_Cdep,false(1,P)])),NaN; %% FIT_M_D    
          mean(R_squared(1,[false(1,J),ind_Bdep])),NaN]; %% FIT_S_D
    INI.GoF=Eval;
    INI.R2_m = R_squared(1,[ind_Cdep,false(1,P)]);
    INI.R2_s = R_squared(1,[false(1,J),ind_Bdep]);

    INI.W=est_W;
    INI.C=est_C;
    INI.B=est_B;

%% (4) Estimation parameters for N_Boot
    if N_Boot<100
        TABLE.W=[est_W(W0),NaN(Nw,5)];
        TABLE.C=[est_C(C0),NaN(Nc,5)];
        TABLE.B=[est_B(B0),NaN(Nb,5)];
        ETC.W_Boot=[];
        ETC.C_Boot=[];
        ETC.B_Boot=[];  
    else
        W_Boot=zeros(Nw,N_Boot);
        C_Boot=zeros(Nc,N_Boot);
        B_Boot=zeros(Nb,N_Boot);
        delOPE_c_Boot=zeros(Nc,N_Boot);
        delOPE_b_Boot=zeros(Nb,N_Boot);

        OPE_Boot=zeros(3,N_Boot);    
        if Flag_Parallel
            parfor b=1:N_Boot
                [Z_ib,Z_oob]=GC_Boot(Z);
                mean_Z_ib=mean(Z_ib);
                std_Z_ib=std(Z_ib,1);
                [W_b,C_b,B_b,~]=ALS_GSCA_Reg(Z_ib,W,W0,C0,B0,ind_Adep,lambda_optimal,Min_limit,Max_iter,size(Z_ib,1),J,P,T,Jy,Py,loc_Cdep,loc_Bdep);            
                W_Boot(:,b)=W_b(W0);
                C_Boot(:,b)=C_b(C0);
                B_Boot(:,b)=B_b(B0);
        
            %   (4) Predictabiliy 
            %   (4a) Err for model
                %(4a-1) Scaling for Z_oob, est_W
                N_oob=size(Z_oob,1);
                Z_oob=(Z_oob-ones(N_oob,1)*mean_Z_ib)./(ones(N_oob,1)*std_Z_ib);
                CV_oob=Z_oob*W_b;             
                %(4a-2) to estimate Err(b)
                e_oob=[Z_oob,CV_oob]-CV_oob*[C_b,B_b];
                ope_full=sum(e_oob.^2,1)/N_oob;
        
                em_oob=ope_full(:,1:J);
                es_oob=ope_full(:,(J+1):T);
        
                em_dep_oob=em_oob(:,ind_Cdep);
                es_dep_oob=es_oob(:,ind_Bdep);
        
                sum_em_dep_oob = sum(em_dep_oob,2);
                sum_es_dep_oob = sum(es_dep_oob,2);
                OPE_Boot(:,b)=[(sum_em_dep_oob+sum_es_dep_oob)/Ty;sum_em_dep_oob/Jy;sum_es_dep_oob/Py]; % Error_oob for the measurement model

                [delOPE_c_Boot(:,b),delOPE_b_Boot(:,b)]=Gen_delOPE(Z_oob,CV_oob, ...
                                                                                    em_oob,es_oob, ...
                                                                                    C_b,B_b,C0,B0, ...
                                                                                    loc_c_t,loc_b_t,loc_Cdep,loc_Bdep, ...
                                                                                    N_oob,Nc,Nb);
            end
        else
            for b=1:N_Boot
                [Z_ib,Z_oob]=GC_Boot(Z);
                mean_Z_ib=mean(Z_ib);
                std_Z_ib=std(Z_ib,1);
                [W_b,C_b,B_b,~]=ALS_GSCA_Reg(Z_ib,W,W0,C0,B0,ind_Adep,lambda_optimal,Min_limit,Max_iter,size(Z_ib,1),J,P,T,Jy,Py,loc_Cdep,loc_Bdep);            
                W_Boot(:,b)=W_b(W0);
                C_Boot(:,b)=C_b(C0);
                B_Boot(:,b)=B_b(B0);
        
            %   (4) Predictabiliy 
            %   (4a) Err for model
                %(4a-1) Scaling for Z_oob, est_W
                N_oob=size(Z_oob,1);
                Z_oob=(Z_oob-ones(N_oob,1)*mean_Z_ib)./(ones(N_oob,1)*std_Z_ib);
                CV_oob=Z_oob*W_b;             
                %(4a-2) to estimate Err(b)
                e_oob=[Z_oob,CV_oob]-CV_oob*[C_b,B_b];
                ope_full=sum(e_oob.^2,1)/N_oob;
        
                em_oob=ope_full(:,1:J);
                es_oob=ope_full(:,(J+1):T);
        
                em_dep_oob=em_oob(:,ind_Cdep);
                es_dep_oob=es_oob(:,ind_Bdep);
        
                sum_em_dep_oob = sum(em_dep_oob,2);
                sum_es_dep_oob = sum(es_dep_oob,2);
                OPE_Boot(:,b)=[(sum_em_dep_oob+sum_es_dep_oob)/Ty;sum_em_dep_oob/Jy;sum_es_dep_oob/Py]; % Error_oob for the measurement model
    
                [delOPE_c_Boot(:,b),delOPE_b_Boot(:,b)]=Gen_delOPE(Z_oob,CV_oob, ...
                                                                                     em_oob,es_oob, ...
                                                                                     C_b,B_b,C0,B0, ...
                                                                                     loc_c_t,loc_b_t,loc_Cdep,loc_Bdep, ...
                                                                                     N_oob,Nc,Nb);
            end
        end

    %% (5) Calculation of statistics
    % Predictiablity
        Eval(:,2)=[mean(OPE_Boot(1,:),2);mean(OPE_Boot(2,:),2);mean(OPE_Boot(3,:),2)];        
        INI.GoF=Eval;
    % CI
        alpha=.05;
        CI=[alpha/2,alpha,1-alpha,1-(alpha/2)];
        loc_CI=round(CI*(N_Boot-1))+1; % .025 .05 .95 .975
          
    % basic statistics for parameter
        TABLE.W=para_stat(est_W(W0),W_Boot,loc_CI,[]);
        if Jy>0; TABLE.C=para_stat(est_C(C0),C_Boot,loc_CI,delOPE_c_Boot); end
        if Py>0; TABLE.B=para_stat(est_B(B0),B_Boot,loc_CI,delOPE_b_Boot); end
        ETC.W_Boot=W_Boot;
        ETC.C_Boot=C_Boot;
        ETC.B_Boot=B_Boot;  
    end
Results.INI=INI;
Results.TABLE=TABLE;
Results.ETC=ETC;
end
function Table=para_stat(est_mt,boot_mt,CI_mp,delOPE_p_q_Boot)
    flag_pet = true;
    if isempty(delOPE_p_q_Boot); flag_pet = false; end
    boot_mt=sort(boot_mt,2);
    SE=std(boot_mt,0,2);
    Table=[est_mt,SE,boot_mt(:,CI_mp(1,1)),boot_mt(:,CI_mp(1,4))]; 
    if flag_pet
        delOPE_p_q_Boot=sort(delOPE_p_q_Boot,2);
        Table=[Table,delOPE_p_q_Boot(:,CI_mp(1,2)),mean(delOPE_p_q_Boot<=0,2)];
    end
end
function [in_sample,out_sample,index,N_oob]=GC_Boot(Data)
    N=size(Data,1); 
    index=ceil(N*rand(N,1));
    in_sample=Data(index,:); 
    index_oob=(1:N)'; index_oob(index)=[];
    out_sample=Data(index_oob,:);
    N_oob=length(index_oob);
end
function [delOPE_c,delOPE_b]=Gen_delOPE(Z_oob,CV_oob, ...
                                         em_oob,es_oob, ...
                                         C_b,B_b,C0,B0, ...
                                         loc_c_t,loc_b_t,loc_Cdep,loc_Bdep, ...
                                         N_oob,Nc,Nb)
    %{
    i_w=0;
    delOPE_w=zeros(Nw,1);
    loc_w_t_set=loc_w_t; % for parfor
    for p=1:P
%   (4b) VIMP for each parameter
        wp = W_b(:,p);
        CVp = CV_oob(:,p);    
        for j=loc_w_t_set{1,p}
            i_w=i_w+1;                
            ind_wp_mj=W0(:,p);
            ind_wp_mj(j)=false;
            e_p_q_k=CVp-Z_oob(:,ind_wp_mj)*wp(ind_wp_mj,1);
            ope_p_q_k=(sum(e_p_q_k.^2,1)/N_oob);
            delOPE_w(i_w,1)=ope_p_q_k;  % for parfor
        end
    end
    %}
    i_c=0;
    delOPE_c=zeros(Nc,1);
    loc_c_t_set=loc_c_t; % for parfor
    for jy_dep=loc_Cdep
%   (4b) VIMP for each parameter
        c_jy = C_b(:,jy_dep);
        zy = Z_oob(:,jy_dep);    
        ope_jy_k = em_oob(1,jy_dep);
        for j=loc_c_t_set{1,jy_dep}
            i_c=i_c+1;                
            ind_cq_p0=C0(:,jy_dep);
            ind_cq_p0(j)=false;                
            e_j_q_k=zy-CV_oob(:,ind_cq_p0)*c_jy(ind_cq_p0,1);
            ope_j_q_k=(sum(e_j_q_k.^2,1)/N_oob);
            delOPE_c(i_c,1)=ope_j_q_k-ope_jy_k;  % for parfor
        end
    end

    i_b=0;
    delOPE_b=zeros(Nb,1);
    loc_b_t_set=loc_b_t; % for parfor
    for q_dep=loc_Bdep
%   (4b) VIMP for each parameter
        bq = B_b(:,q_dep);
        ry = CV_oob(:,q_dep);    
        ope_q_k = es_oob(1,q_dep);
        for p=loc_b_t_set{1,q_dep}
            i_b=i_b+1;                
            ind_bq_p0=B0(:,q_dep);
            ind_bq_p0(p)=false;                
            e_p_q_k=ry-CV_oob(:,ind_bq_p0)*bq(ind_bq_p0,1);
            ope_p_q_k=(sum(e_p_q_k.^2,1)/N_oob);
            delOPE_b(i_b,1)=ope_p_q_k-ope_q_k;  % for parfor
        end
    end
end