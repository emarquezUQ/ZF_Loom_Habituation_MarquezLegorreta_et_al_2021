function [ModelResults,SelectedSamples] = Test_Regress(DataSet,Regressor,idxKmeans,varargin)
%function [ModelResults] = Test_Regress(DataSet,Regressor,varargin)
% Test DataSet against Regressor with multiple linear regression
% DataSet and Regressor of the form (Samples x Time)
numvarargs = length(varargin);
if numvarargs > 1
    error('Test_Regress:TooManyInputs', ...
        'requires at most 1 optional inputs');
end

optargs = {0.5};
optargs(1:numvarargs) = varargin;
[rsq_limit] = optargs{:};

parfor i=1:size(DataSet,1)
    %mdl=stepwiselm(Regressor',DataSet(i,:),'linear','Criterion','adjrsquared','Intercept',false,'Upper','interactions','Verbose',0);
    mdl=stepwiselm(Regressor',DataSet(i,:),'linear','Criterion','bic','Upper','linear','Verbose',0);
    ModelResults(i).coef=mdl.Coefficients;
    ModelResults(i).MSE=mdl.MSE;
    ModelResults(i).Fitted=mdl.Fitted;
    ModelResults(i).rsquared=mdl.Rsquared.Adjusted;
end

[SelectedSamples]=Draw_Regress(ModelResults,DataSet,idxKmeans,rsq_limit);
