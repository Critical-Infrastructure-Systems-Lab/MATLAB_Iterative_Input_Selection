function [model] = crossvalidation_extra_tree_ensemble(subset,M,k,nmin,ns,flag,problemType,inputType)
% This function cross-validates an ensemble of Exra-Trees 
%
% Inputs: 
% subset = observations
% M      = number of trees in the ensemble
% k      = number of random cut-directions 
% nmin   = minimum number of points per leaf 
% ns     = number of folds in the k-fold cross-validation process 
% flag   = if flag == 1, the model is then evaluated (and saved) on the full dataset
% problemType = specify problem type (1 for regression, zero for classification)
% inputType   = binary vector indicating feature type(0:categorical,1:numerical)
% only include input type for classification problems
%
% Output: 
% model  = structure containing models and performance 
%
% Copyright 2015 Ahmad Alsahaf
% Research fellow, Politecnico di Milano
% ahmadalsahaf@gmail.com
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

if problemType == 0
    [model] = crossvalidation_extra_tree_ensemble_r(subset,M,k,nmin,ns,flag);
    
else
    [model] = crossvalidation_extra_tree_ensemble_c(subset,M,k,nmin,ns,inputType,flag);

end

