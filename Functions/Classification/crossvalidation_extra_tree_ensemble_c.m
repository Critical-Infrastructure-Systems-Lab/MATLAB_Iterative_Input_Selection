function [model] = crossvalidation_extra_tree_ensemble_c(subset,M,k,nmin,ns,inputType,sampleWeights,flag)


% This function cross-validates an ensemble of Exra-Trees 
%
% Inputs:
% subset    = observations
% M         = number of trees in the ensemble
% k         = number of random cut-directions 
% nmin      = minimum number of points per leaf 
% ns        = number of folds in the k-fold cross-validation process 
% inputType = binary vector indicating feature type (0:categorical, 1:numerical)
% flag       = if flag == 1, the model is then evaluated (and saved) on the
% full dataset
%
% Output: 
% model  = structure containing models and performance 
%
%
%
% Copyright 2015 Ahmad Alsahaf
% Research fellow, Politecnico di Milano
% ahmadalsahaf@gmail.com
% Copyright 2014 Stefano Galelli
% Assistant Professor, Singapore University of Technology and Design
% stefano_galelli@sutd.edu.sg
% http://people.sutd.edu.sg/~stefano_galelli/index.html
% Please refer to README.txt for further information.
%
%
% This file is part of MATLAB_IterativeInputSelection.
% 
%     MATLAB_IterativeInputSelection is free software: you can redistribute 
%     it and/or modify it under the terms of the GNU General Public License 
%     as published by the Free Software Foundation, either version 3 of the 
%     License, or (at your option) any later version.     
% 
%     This code is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with MATLAB_IterativeInputSelection.  
%     If not, see <http://www.gnu.org/licenses/>.ction




% 0) SET THE PROBLEM PARAMETERS FOR THE ENSEMBLE CROSS-VALIDATION

% Number of lines characterizing an alternative (a single fold)
l = floor(length(subset)/ns);

% Re-define the subset matrix and the sampleWeights
subset = subset(1:l*ns,:);
sampleWeights = sampleWeights(1:l*ns);

% 1) INITIALIZATION OF THE OUTPUT VARIABLES 

% Initialize Performance VECTORS
classPerf_cal_pred = zeros(ns,1); 
classPerf_val_pred = zeros(ns,1); 

% Initialize the function output
model.cross_validation.performance = [];


% 2) MODEL CONSTRUCTION AND EVALUATION OF THE PERFORMANCES (k-fold cross-validation)
% Counter
% disp('Start cross-validation:')
for i = 1:ns

    % Counter
    % disp('Start cross-validation:'); disp(i);

    % Define the calibration and validation data-set
    % Calibration
    if (i > 1) && (i < ns)
        subset_tar = [subset(i*l+1:end,:) ; subset(1:(i-1)*l,:)];
        sampleWeights_tar = [sampleWeights(i*l+1:end,:) ; sampleWeights(1:(i-1)*l,:)];
    else if i == 1
            subset_tar = subset(i*l+1:end,:);
            sampleWeights_tar = sampleWeights(i*l+1:end,:);
        else
            subset_tar = subset(1:(i-1)*l,:);
            sampleWeights_tar = sampleWeights(1:(i-1)*l,:);
        end
    end
    
    % Validation
    subset_val = subset((i-1)*l+1:i*l,:);
    sampleWeights_val = sampleWeights((i-1)*l+1:i*l,:);

    % Ensemble building + test the ensemble on the calibration dataset
    [ensemble,finalResult_cal_pred] = buildAnEnsemble_c(M,k,nmin,subset_tar,inputType,sampleWeights_tar);
    classPerf_cal_pred(i)           = class_perf(subset_tar(:,end),finalResult_cal_pred,sampleWeights_tar);    
    
    % Test the ensemble on the validation data-set     
    [finalResult_val_pred] = predictWithAnEnsemble_c(ensemble,subset_val(:,1:end-1));        
    classPerf_val_pred(i)        = class_perf(subset_val(:,end),finalResult_val_pred,sampleWeights_val);

end

% Average performance
model.cross_validation.performance.classPerf_cal_pred = classPerf_cal_pred;
model.cross_validation.performance.classPerf_val_pred = classPerf_val_pred;
model.cross_validation.performance.classPerf_cal_pred_mean = mean(classPerf_cal_pred);
model.cross_validation.performance.classPerf_val_pred_mean = mean(classPerf_val_pred);


% 3) MODEL CONSTRUCTION ON THE WHOLE DATA-SET

% Check if is necessary to test the model on the whole data-set
if flag == 1

    % Add new fields to the function output
    model.complete_model.ensemble      = [];
    model.complete_model.trajectories  = [];
    model.complete_model.performance   = [];

    % Counter
    % disp('Building and testing the final model');

    % Model construction
    [ensemble,finalResult_pred] = buildAnEnsemble_c(M,k,nmin,subset,inputType,sampleWeights);
    model.complete_model.ensemble = ensemble;
    model.complete_model.trajectories = finalResult_pred;

    % Evaluate performance             
    [model.complete_model.performance.classPerf, ...
       model.complete_model.performance.misClassified, ...
       model.complete_model.performance.Classified] ...
       = class_perf(subset(:,end),finalResult_pred);
    
else
    
    return
    
end







