function    R = Rt2_fit(yo,ys)

% This function calculates the value of the coefficient of determinarion R2
%
% Input:
% yo = observed data
% ys = simulated (or predicted) data
%
% Output
% R = coefficient of determination
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



if(nargin ~= 2)
  disp(  'error: wrong number of inputs'  )
  return;
end;


R = 1-[(cov(yo-ys)/cov(yo))];


% This code has been written by Stefano Galelli