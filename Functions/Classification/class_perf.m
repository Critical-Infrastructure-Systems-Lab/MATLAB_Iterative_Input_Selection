function    [R, idx1, idx2] = class_perf(yo,ys,sampleWeights)

% This function calculates measures the performance of classification
% as the percentace of correct classifications (to be modified)
%
% Input:
% yo = observed data
% ys = simulated (or predicted) data
%
% Output
% R   = percentage of correct classifications
% idx1 = indices of incorrect classifications
% idx2 = indices of correct classifications
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



if(nargin >3)
  disp(  'error: wrong number of inputs'  )
  return;
end;


temp = ys - yo;

if nargin <3
    R = numel(find(~temp))/numel(yo);
else
    R = sum(sampleWeights(find(~temp)))/sum(sampleWeights);
end


% find indices of mis-classified and classified entries.
idx1 = find(temp);           
idx2 = find(~temp); 





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% classes  = unique(yo);
% nClasses = numel(unique(yo));
% if nClasses == 2
%     nClasses = 1; % handle binary classification
% end
%     
% 
% Acc = zeros(1,nClasses);
% for j = 1 : nClasses
% %     if nClasses == 1
% %         thisClass = 2; % the H1
% %     else
%         thisClass = classes(j);
% %     end
%     % ixes 
%     ixes1 = (yo == thisClass);        
%     ixes2 = (ys == thisClass); 
%     % compute confusion matrix
%     tp = sum((ixes1==ixes2)&(ixes1==1));
%     tn = sum((ixes1==ixes2)&(ixes1==0));
%     fn = sum((ixes1-ixes2)==1);
%     fp = sum((ixes1-ixes2)==-1);    
%     % compute accuracy
%     Acc(j) = (tp+tn)/(tp+fn+fp+tn);
% end
% % get average accuracy
% R = mean(Acc);
% 
% 




