function [model] = crossvalidation_extra_tree_ensemble(subset,M,k,nmin,ns,flag)


% This function cross-validate an ensemble of Exra-Trees 
%
% Inputs: 
% subset = observations
% M      = number of trees in the ensemble
% k      = number of random cut-directions 
% nmin   = minimum number of points per leaf 
% ns     = number of folds in the k-fold cross-validation process 
% flag   = if flag == 1, the model is then evaluated (and saved) on the
% full dataset
%
% Output: 
% model  = structure containing models and performance 
%
%
% Copyright 2014 Stefano Galelli
% Assistant Professor, Singapore University of Technology and Design
% stefano_galelli@sutd.edu.sg
% http://people.sutd.edu.sg/~stefano_galelli/index.html
%
%
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
%     If not, see <http://www.gnu.org/licenses/>.
% 



% 0) SET THE PROBLEM PARAMETERS FOR THE ENSEMBLE CROSS-VALIDATION

% Number of lines characterizing an alternative (a single fold)
l = floor(length(subset)/ns);

% Re-define the subset matrix
subset = subset(1:l*ns,:);


% 1) INITIALIZATION OF THE OUTPUT VARIABLES 

% Initialize R2 AND RMSE VECTORS
Rt2_cal_pred = zeros(ns,1); 
Rt2_val_pred = zeros(ns,1); 

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
    else if i == 1
            subset_tar = subset(i*l+1:end,:);
        else
            subset_tar = subset(1:(i-1)*l,:);
        end
    end
    % Validation
    subset_val = subset((i-1)*l+1:i*l,:);

    % Ensemble building + test the ensemble on the calibration dataset
    [ensemble,finalResult_cal_pred] = buildAnEnsemble(M,k,nmin,subset_tar);
     Rt2_cal_pred(i)                = Rt2_fit(subset_tar(:,end),finalResult_cal_pred);    
    
    % Test the ensemble on the validation data-set     
    [finalResult_val_pred] = predictWithAnEnsemble(ensemble,subset_val);        
    Rt2_val_pred(i)        = Rt2_fit(subset_val(:,end),finalResult_val_pred);

end

% Average R2
model.cross_validation.performance.Rt2_cal_pred = Rt2_cal_pred;
model.cross_validation.performance.Rt2_val_pred = Rt2_val_pred;
model.cross_validation.performance.Rt2_cal_pred_mean = mean(Rt2_cal_pred);
model.cross_validation.performance.Rt2_val_pred_mean = mean(Rt2_val_pred);


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
    [ensemble,finalResult_pred] = buildAnEnsemble(M,k,nmin,subset);
    model.complete_model.ensemble = ensemble;
    model.complete_model.trajectories = finalResult_pred;

    % Evaluate R2              
    model.complete_model.performance.Rt2   = Rt2_fit(subset(:,end),finalResult_pred);

else
    
    return
    
end


% This code has been written by Stefano Galelli.
 









