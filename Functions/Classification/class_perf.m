function    [R, idx] = class_perf(yo,ys)

% This function calculates measures the performance of classification
% as the percentace of correct classifications (to be modified)
%
% Input:
% yo = observed data
% ys = simulated (or predicted) data
%
% Output
% R   = percentage of correct classifications
% idx = indices of incorrect classifications
%
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



if(nargin ~= 2)
  disp(  'error: wrong number of inputs'  )
  return;
end;


temp = ys - yo;
R = numel(find(~temp))/numel(yo);


% find indices of mis-classified entries.
idx = find(temp);           


% This code has been written by Ahmad Alsahaf