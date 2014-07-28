function [subset_out,idx] = shuffle_data(subset_in)

% This function builds an ensemble of Exra-Trees and then ranks 
% the input variables according to their importance
%
% Input: 
% subset = observations
%
% Output: 
% subset_out = shuffled observations
% idx = indexes used in the permutation
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



% create a random permutation
[N,M] = size(subset_in);
idx   = randperm(N);
idx   = idx';

% initialize the output vector
subset_out = nan(size(subset_in));

% shuffle
for j = 1:N
    subset_out(j,:) = subset_in(idx(j),:);
end


% This code has been written by Stefano Galelli.


