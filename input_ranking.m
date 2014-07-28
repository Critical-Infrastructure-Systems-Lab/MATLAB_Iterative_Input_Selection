function [result] = input_ranking(subset,M,k,nmin)

% This function builds an ensemble of Exra-Trees and then ranks 
% the input variables according to their importance
%
% Inputs:
% subset = observations 
% M      = number of trees 
% k      = number of random cut-directions 
% nmin   = minimum number of points per leaf
%
% Output: 
% result = ranked score of each attribute
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



% Build and ensemble of Extra Trees and get the score of each variable (for
% each tree)
rank_results = zeros(M,size(subset,2)-1);
for i = 1 : M
    [tree,output,rank_results(i,:)] = buildAnExtraTree(k,nmin,subset);   
end

% Compute, for each Extra-Tree, the score of each variable
cum_score = sum(rank_results,2);
for i = 1:(size(subset,2)-1)
    rank_results_(:,i) = (rank_results(:,i)./cum_score)*100;
end

% Compute the average (over M) score of each variable
contr = mean(rank_results_,1);

% Rank the variable in decreasing order
[X,I] = sort(contr',1,'descend');
result = [X,I];


% This code has been written by Stefano Galelli.





