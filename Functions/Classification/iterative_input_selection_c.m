
function [result] = iterative_input_selection_c(subset,M,nmin,ns,p,inputType,epsilon,max_iter)

% This function implements the IIS algorithm for classification problems
%
% subset   = observations
% M        = number of trees in the ensemble
% nmin     = minimum number of points per leaf 
% ns       = number of folds in the k-fold cross-validation process 
% p        = number of SISO models (it must be smaller than the number of
%            candidate inputs).
% inputType = binary vector indicating feature type (0:categorical, 1:numerical)
% epsilon  = tolerance
% max_iter = maximum number of iterations
%
%
% Output
% result   = structure containing the result for each iteration
%
% Copyright 2015 Ahmad Alsahaf
% Research fellow, Politecnico di Milano
% ahmadalsahaf@gmail.com
%
% Copyright 2014 Stefano Galelli and Matteo Giuliani
% Assistant Professor, Singapore University of Technology and Design
% stefano_galelli@sutd.edu.sg
% http://people.sutd.edu.sg/~stefano_galelli/index.html
% Research Fellow, Politecnico di Milano
% matteo.giuliani@polimi.it
% http://giuliani.faculty.polimi.it
%
% Copyright 2015 Ahmad Alsahaf
% Research fellow, Politecnico di Milano
% ahmadalsahaf@gmail.com
%
% Please refer to README.txt for bibliographical references on Extra-Trees!
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

% 0) SET THE PARAMETERS

% Initialize the counter and the exit condition flag
iter     = 1;    % iterations counter
diff     = 1;    % exit flag 

% Re-define the subset matrix
l = floor(length(subset)/ns);
subset = subset(1:l*ns,:);


% Define the MISO model output
miso_output = subset(:,end);

% Define the set of candidate input variables
input  = subset(:,1:end-1);

% Other variables to be initialized      
miso_input = [];  % initialize an empty set to store the input variables to be selected%


% 1) IIS ALGORITHM
while (diff > epsilon) && (iter <= max_iter)
    
    % Visualize the iteration
    disp('ITERATION:'); disp(iter);
    
    % Define the output variable to be used during the ranking
    rank_output = miso_output; 
    rank_input  = input;
    
    % Define sample weights
    if iter == 1       
        sampleWeights = (1/(l*ns))*ones(l*ns,1);              % first iteration, weight samples evenly
    else         
        sampleWeights = new_weights*(1/sum(new_weights));    % weight according to previous classifier performance    
    end
   
    
    % Define the ranking matrix
    matrix_ranking = [rank_input rank_output];
    
    % Run the feature ranking
    disp('Ranking:');
    k = size(input,2);
    [ranking] = input_ranking_c(matrix_ranking,M,k,nmin,inputType,sampleWeights);     
    eval(['result.iter_' num2str(iter) '.ranking' '=' 'ranking;']);
    disp(ranking);
    
    % Select and cross-validate p SISO models (the first p-ranked models)
    disp('Evaluating SISO models:');
    features = ranking(1:p,2);                             % p features to be considered           
    performance = zeros(p,1);	                           % initialize a vector for the performance of the p SISO models%	
    
    for i = 1:p
        [siso_model] = crossvalidation_extra_tree_ensemble_c([subset(:,features(i)) rank_output(1:l*ns,:)],M,1,nmin,ns,inputType,sampleWeights,0);        
		performance(i) = siso_model.cross_validation.performance.classPerf_val_pred_mean;
    end
    eval(['result.iter_' num2str(iter) '.SISO' '=' '[features performance];']);
    disp([features performance]);
    
    % Choose the SISO model with the best performance
	[val,idx_siso] = max(performance);
	best_siso_input = features(idx_siso);
    eval(['result.iter_' num2str(iter) '.best_SISO' '=' '[best_siso_input val];']);
    disp('Select variable:'); disp(best_siso_input);
    
    % Check the exit condition 
    if (all(miso_input - best_siso_input) == 0) 
        result.exit_condition = 'An input variable was selected twice';
        result.iters_done = iter;
        return
    end
        
%     % force algorithms to keep going
%     if (all(miso_input - best_siso_input) == 0) 
%         miso_input = unique(miso_input);
%     end



	% Build a MISO model with the selected inputs	
    disp('Evaluating MISO model:');
	miso_input = [miso_input best_siso_input];				 
	k = length(miso_input);	
    misoWeight = (1/(ns*l))*ones(ns*l,1);
	[miso_model] = crossvalidation_extra_tree_ensemble_c([subset(:,miso_input) miso_output],M,k,nmin,ns,inputType,misoWeight,1);				 
    eval(['miso_model_' num2str(iter) '= miso_model;']);
    eval(['result.iter_' num2str(iter) '.MISO' '=' 'miso_model;']);
    disp(miso_model.cross_validation.performance.classPerf_val_pred_mean);
	
    % Evaluate the performance of the MISO model and calculate the difference with respect to the previous MISO model
    if iter == 1   % at the first iteration, use a default value
        diff = 1;
    else
        diff = miso_model.cross_validation.performance.classPerf_val_pred_mean ...
            - eval(['miso_model_' num2str(iter-1) '.cross_validation.performance.classPerf_val_pred_mean'])
    end	
        
	% Compute the MISO model residual by weighting classified samples
    residual_idx = miso_model.complete_model.performance.misClassified;
	correct_idx = miso_model.complete_model.performance.Classified;

    new_weights = zeros(ns*l,1);
    eps_r = 1-miso_model.complete_model.performance.classPerf;  % error = ratio of misclassification weights
    beta_r = eps_r/(1-eps_r);                                   % odds of a fake classification
    new_weights(residual_idx)=sampleWeights(residual_idx);      % weight of misclassified samples fixed
    new_weights(correct_idx)=sampleWeights(correct_idx)*beta_r; % multiply classified sample weights by beta_r
    
    
	% Update the counter				 
	iter = iter + 1;
    
    % Check the exit condition
    if iter > max_iter 
        result.exit_condition = 'The maximum number of iterations was reached'; 
        result.iters_done = iter-1;
    end
    if diff <= epsilon  
        result.exit_condition = 'The tolerance epsilon was reached';    
        result.iters_done = iter-1;
    end
    
end

    
% This code has been written by Stefano Galelli, Matteo Giuliani, Ahmad Alsahaf















