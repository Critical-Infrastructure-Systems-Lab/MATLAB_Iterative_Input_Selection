
% This script shows how to use the functions available in the
% MATLAB_IterativeInputSelection toolbox on a classification sample dataset. These
% functions require the MATLAB_ExtraTrees_classification toolbox, which can be found at
% https://github.com/rtaormina/MATLAB_ExtraTrees.
%
%
% Copyright 2015 Ahmad Alsahaf
% Research fellow, Politecnico di Milano
% ahmadalsahaf@gmail.com
%
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


% Test with census dataset from UCI machine learning repository
% https://archive.ics.uci.edu/ml/datasets/Adult

clear all
clc

% import data
data = importdata('census.mat');
[n m] = size(data);
data = data(randperm(n),:);

% split dataset into training and testing
trainSet = data(1:1000,:);
testSet = data(2001:end,:);

%  %Add aritifial features to the data set before IIS
% nff = 1;     %number of artificial features to add
% fake_features = zeros(size(trainSet,1),nff);
% for i=1:nff
%     x = randi([1 4]);
%     fake_features(:,i) = randi([0 x],size(trainSet,1),1);
% end
% 
% 
% % append artificial features after real features
% trainSet_and_fake = [trainSet(:,1:end-1) fake_features trainSet(:,end)];


% IIS algorithm training parameters
subset = trainSet;
M = 50;
nmin = 3;
ns = 10;
p = 4;
epsilon = eps;
max_iter = 10;
problemType = 1;
Vflag = 1;
inputType = logical(zeros(size(subset,2)-1,1));


% Run IIS algorithm 
result = iterative_input_selection(subset,M,nmin,ns,p,epsilon,max_iter,problemType,Vflag,inputType)



%------------------------------------------------------------------------------------------------------------------
%%  Second example:  Heart disease data from UCI repository

clear all
clc

% import and define dataset
data = csvread('Heart.csv',1,0);
subset = data(2:end,:);
inputType = logical(data(1,1:end-1));
[n m] = size(subset);
subset = subset(randperm(n),:);


% IIS algorithm training parameters
M = 100;
nmin = 3;
ns = 10;
p = 5;
epsilon = -1;
max_iter = 10;
problemType = 1;
Vflag = 1;



% Run IIS algorithm 
result = iterative_input_selection(subset,M,nmin,ns,p,epsilon,max_iter,problemType,Vflag,inputType)

