%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Illustration for GSCA.Reg_Prime package                                 %
%   Author: Gyeongcheol Cho                                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description:                                                            %
%   - This code aims to illustrate how to use GSCA.Reg_Prime package.     %
%   - The dataset is a replica of the ACSI data used in Cho (submitted).  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% References                                                              %
%     * Hwang, H. Regularized Generalized Structured Component Analysis.  %
%         Psychometrika 74, 517–530 (2009).                               %
%         https://doi.org/10.1007/s11336-009-9119-y                       %
%     * Cho, G. (submitted). Predictor exclusion threshold: A criterion   %
%         for determining predictor relevance in regularized generalized  % 
%         structured component analysis.                                  %   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Data=readtable('ACSI_Comp_Replica.csv');
W0=[1 1 1 0 0 0 0 0 0 0 0 0 0 0 ; ...
   0 0 0 1 1 1 0 0 0 0 0 0 0 0 ; ...
   0 0 0 0 0 0 1 1 0 0 0 0 0 0 ; ...
   0 0 0 0 0 0 0 0 1 1 1 0 0 0 ; ...
   0 0 0 0 0 0 0 0 0 0 0 1 0 0 ; ...
   0 0 0 0 0 0 0 0 0 0 0 0 1 1 ]';
C0=W0';
B0=[0 1 1 1 0 0;...
   0 0 1 1 0 0;...
   0 0 0 1 0 0;...
   0 0 0 0 1 1;...
   0 0 0 0 0 1;...
   0 0 0 0 0 0];
N_Boot=1000;
vecLambda0=[0 .1 .5 1 5 10 100];
KxK=[1,5];
Max_iter = 1000;
Min_limit = 10^(-8);
Flag_Parallel = false;
[INI,TABLE,ETC]=RegGSCA(Data{:,:},W0,C0,B0,vecLambda0,KxK,N_Boot,Max_iter,Min_limit,Flag_Parallel);
INI
INI_reg.minVE
INI_reg.LamVE
INI_reg.lambda_optimal
INI.GoF
INI.W
INI.C
INI.B
INI.R2_m
INI.R2_s
TABLE
TABLE.W
TABLE.C
TABLE.B
ETC