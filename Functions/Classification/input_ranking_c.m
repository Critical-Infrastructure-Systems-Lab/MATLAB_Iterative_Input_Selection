function [result] = input_ranking_c(subset,M,k,nmin,inputType)

% This function builds an ensemble of Exra-Trees and then ranks 
% the input variables according to their importance
%
% Inputs:
% subset = observations 
% M      = number of trees 
% k      = number of random cut-directions 
% nmin   = minimum number of points per leaf
% inputType = binary vector indicating feature type (0:categorical, 1:numerical)
%
% Output: 
% result = ranked score of each attribute
%
%
% Copyright 2015 Ahmad Alsahaf
% Research fellow, Politecnico di Milano
% ahmadalsahaf@gmail.com
%
%
% Copyright 2014 Stefano Galelli
% Assistant Professor, Singapore University of Technology and Design
% stefano_galelli@sutd.edu.sg
% http://people.sutd.edu.sg/~stefano_galelli/index.html
%
% Copyright 2015 Ahmad Alsahaf
% Research fellow, Politecnico di Milano
% ahmadalsahaf@gmail.com
%
% Please refer to README.txt for more information
%
% This file is part of MATLAB_IterativeInputSelection_classification, an input 
% selection in classification problems. For the same toolbox in regression
% problems, go to: https://github.com/stefano-galelli/MATLAB_IterativeInputSelection





% Build and ensemble of Extra Trees and get the score of each variable (for
% each tree)

rank_results = zeros(M,size(subset,2)-1);
for i = 1 : M
    [tree,output,rank_results(i,:)] = buildAnExtraTree_c(k,nmin,subset,inputType);   
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


% This code has been written by Stefano Galelli





