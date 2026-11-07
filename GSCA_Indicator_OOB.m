function stats = GSCA_Indicator_OOB(Z_ib_raw,Z_oob_raw,W_b,C_b,B_b,C0,B0)
%GSCA_INDICATOR_OOB Single-fit OOB path loss using attached observed indicators.
%   stats = GSCA_Indicator_OOB(Z_ib_raw,Z_oob_raw,W_b,C_b,B_b,C0,B0)
%
%   The supplied W_b, C_b, and B_b are kept unchanged. Standardization uses
%   the in-bag mean and population standard deviation, as in GSCA_Reg.
%   For each specified B(p,q), regress in-bag CV_p on an intercept and the
%   OTHER specified component predictors of q. Evaluate this regression on
%   the OOB parent scores and replace only CV_p's contribution to q.
%
%   Full and replaced target-component predictions are mapped through the
%   supplied target loadings to predictions of q's attached indicators.
%   Evaluate squared errors against the fixed RAW observed indicator values,
%   then average over indicators and OOB observations. The observed target
%   indicators themselves are never inputs to their prediction.
%
%   delOPE_b follows find(B0), exactly as B_b(B0~=0). This helper introduces
%   no fitting/bootstrap loops, correction, weighting change, or normalization
%   by EPE. It supports simple disjoint measurement blocks; crossloadings,
%   self paths, and score weights using another block's indicators are
%   rejected explicitly.
%   For C edges, delOPE_c retains the original STANDARDIZED-indicator
%   reconstruction loss and uses the same in-bag auxiliary replacement.
validateattributes(Z_ib_raw,{'numeric'},{'2d','real','finite','nonempty'});
validateattributes(Z_oob_raw,{'numeric'},{'2d','real','finite','nonempty'});
validateattributes(W_b,{'numeric'},{'2d','real','finite','nonempty'});
validateattributes(C_b,{'numeric'},{'2d','real','finite'});
validateattributes(B_b,{'numeric'},{'2d','real','finite'});
[J,P] = size(W_b);
N_ib = size(Z_ib_raw,1);
N_oob = size(Z_oob_raw,1);
if N_ib<2 || size(Z_ib_raw,2)~=J || size(Z_oob_raw,2)~=J || ...
        ~isequal(size(C_b),[P,J]) || ~isequal(size(B_b),[P,P]) || ...
        ~isequal(size(C0),[P,J]) || ~isequal(size(B0),[P,P])
    error('GSCA_Indicator_OOB:Dimensions','Input dimensions are incompatible.');
end
if (~isnumeric(C0) && ~islogical(C0)) || ...
        (~isnumeric(B0) && ~islogical(B0)) || ~isreal(C0) || ~isreal(B0) || ...
        any(~isfinite(C0(:))) || any(~isfinite(B0(:)))
    error('GSCA_Indicator_OOB:Masks','Model masks must be finite and real.');
end
Cmask = C0~=0;
Bmask = B0~=0;
if any(sum(Cmask,1)>1) || any(sum(Cmask,2)==0)
    error('GSCA_Indicator_OOB:UnsupportedModel', ...
        'Each component requires its own nonoverlapping indicator block.');
end
if any(diag(Bmask)) || any(W_b(~Cmask')~=0)
    error('GSCA_Indicator_OOB:UnsupportedModel', ...
        'Self paths and score weights using another block''s indicators are unsupported.');
end
if any(C_b(~Cmask)~=0) || any(B_b(~Bmask)~=0)
    error('GSCA_Indicator_OOB:UnsupportedModel','Supplied coefficients violate the model masks.');
end
mean_ib = mean(Z_ib_raw,1);
std_ib = std(Z_ib_raw,1,1);
if any(std_ib<=0)
    error('GSCA_Indicator_OOB:ConstantIndicator','An in-bag indicator has zero standard deviation.');
end
Z_ib = bsxfun(@rdivide,bsxfun(@minus,Z_ib_raw,mean_ib),std_ib);
Z_oob = bsxfun(@rdivide,bsxfun(@minus,Z_oob_raw,mean_ib),std_ib);
CV_ib = Z_ib*W_b;
CV_oob = Z_oob*W_b;
[source,target] = find(Bmask);
Nb = numel(source);
stats.delOPE_b = nan(Nb,1);
stats.epeFull_by_target = nan(P,1);
stats.epeRemoved_by_edge = nan(Nb,1);
stats.source = source;
stats.target = target;
stats.edgeIndex = find(Bmask);
stats.indicatorCount = sum(Cmask(target,:),2);
stats.indicatorIndices = cell(P,1);
stats.fullPrediction_by_target = cell(P,1);
stats.removedPrediction_by_edge = cell(Nb,1);
stats.auxiliaryCoefficients = cell(Nb,1);
stats.auxiliaryParents = cell(Nb,1);
stats.replacementScore = cell(Nb,1);
stats.trainingMean = mean_ib;
stats.trainingSD = std_ib;
stats.N_oob = N_oob;

% Measurement C edges retain the original standardized outcome/loss scale.
[Csource,Ctarget] = find(Cmask);
Nc = numel(Csource);
stats.delOPE_c = nan(Nc,1);
stats.epeFull_C_by_indicator = nan(J,1);
stats.epeRemoved_C_by_edge = nan(Nc,1);
stats.Csource = Csource;
stats.Ctarget = Ctarget;
stats.CedgeIndex = find(Cmask);
stats.auxiliaryCoefficients_C = cell(Nc,1);
stats.auxiliaryParents_C = cell(Nc,1);
for j=find(any(Cmask,1))
    parents = find(Cmask(:,j));
    fullPrediction = CV_oob(:,parents)*C_b(parents,j);
    fullError = Z_oob(:,j)-fullPrediction;
    stats.epeFull_C_by_indicator(j) = mean(fullError.^2);
    edges = find(Ctarget==j);
    for e=reshape(edges,1,[])
        p = Csource(e);
        other = parents(parents~=p);
        coefficients = [ones(N_ib,1),CV_ib(:,other)]\CV_ib(:,p);
        replacement = [ones(N_oob,1),CV_oob(:,other)]*coefficients;
        replacedPrediction = fullPrediction+(replacement-CV_oob(:,p))*C_b(p,j);
        removedEPE = mean((Z_oob(:,j)-replacedPrediction).^2);
        stats.epeRemoved_C_by_edge(e) = removedEPE;
        stats.delOPE_c(e) = removedEPE-stats.epeFull_C_by_indicator(j);
        stats.auxiliaryCoefficients_C{e} = coefficients;
        stats.auxiliaryParents_C{e} = other;
    end
end

targets = find(any(Bmask,1));
for q=targets
    indicators = find(Cmask(q,:));
    parents = find(Bmask(:,q));
    if any(any(W_b(indicators,parents)~=0))
        error('GSCA_Indicator_OOB:TargetLeakage', ...
            'Target indicators enter a parent score used for prediction.');
    end
    fullScore = CV_oob(:,parents)*B_b(parents,q);
    fullPrediction = raw_prediction(fullScore,C_b(q,indicators), ...
        mean_ib(indicators),std_ib(indicators));
    observed = Z_oob_raw(:,indicators);
    fullError = observed-fullPrediction;
    stats.epeFull_by_target(q) = mean(fullError(:).^2);
    stats.indicatorIndices{q} = indicators;
    stats.fullPrediction_by_target{q} = fullPrediction;
    edges = find(target==q);
    for e=reshape(edges,1,[])
        p = source(e);
        other = parents(parents~=p);
        % Ordinary least squares with an intercept, fitted on in-bag scores.
        aux_ib = [ones(N_ib,1),CV_ib(:,other)];
        coefficients = aux_ib\CV_ib(:,p);
        replacement = [ones(N_oob,1),CV_oob(:,other)]*coefficients;
        replacedScore = fullScore+(replacement-CV_oob(:,p))*B_b(p,q);
        replacedPrediction = raw_prediction(replacedScore,C_b(q,indicators), ...
            mean_ib(indicators),std_ib(indicators));
        removedError = observed-replacedPrediction;
        removedEPE = mean(removedError(:).^2);
        stats.epeRemoved_by_edge(e) = removedEPE;
        stats.delOPE_b(e) = removedEPE-stats.epeFull_by_target(q);
        stats.removedPrediction_by_edge{e} = replacedPrediction;
        stats.auxiliaryCoefficients{e} = coefficients;
        stats.auxiliaryParents{e} = other;
        stats.replacementScore{e} = replacement;
    end
end
end

function prediction = raw_prediction(score,loadings,location,scale)
prediction = bsxfun(@plus,bsxfun(@times,score*loadings,scale),location);
end

