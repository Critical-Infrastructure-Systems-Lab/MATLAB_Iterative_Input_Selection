function [result] = iterative_input_selection(subset,M,nmin,ns,p,epsilon,max_iter,problemType,Vflag,inputType)
% This function implements the IIS algorithm for regression or classification problems
%
% subset      = observations
% M           = number of trees in the ensemble
% nmin        = minimum number of points per leaf 
% ns          = number of folds in the k-fold cross-validation process 
% p           = number of SISO models (it must be smaller than the number of
%            candidate inputs).
% epsilon     = tolerance
% max_iter    = maximum number of iterations
% problemType = specify problem type (0 for regression, 1 for classification)
% Vflag       = selection of the type of validation, 
%               1 = k-fold(default)
%               2= repeated random sub-sampling
% inputType   = binary vector indicating feature type(0:categorical,1:numerical)
% only include input type for classification problems
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
    % regression problems
    [result] = iterative_input_selection_r(subset,M,nmin,ns,p,epsilon,max_iter,Vflag);
    
else
    % classification problems
    [result] = iterative_input_selection_c(subset,M,nmin,ns,p,inputType,epsilon,max_iter);

end
