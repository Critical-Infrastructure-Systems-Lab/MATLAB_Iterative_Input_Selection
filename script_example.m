
% This script shows how to use the functions available in the
% MATLAB_IterativeInputSelection toolbox on a sample dataset. These
% functions require the MATLAB_ExtraTrees toolbox, which can be found at
% https://github.com/rtaormina/MATLAB_ExtraTrees.
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


%% Set workspace
clear all
clc

% Addpath for the MATLAB_ExtraTrees toolbox
addpath('..\MATLAB_ExtraTrees\');



%% Load and prepare data

% load data
load -ascii Friedman_dataset.txt

% rename data
data = Friedman_dataset;
clear Friedman_dataset

% definition of the calibration and validation data-set
subset_cal = data(1:180,:);      
subset_val = data(181:end,:);



%% Set the parameters for the Extra-Trees

M    =  25;  % Number of trees in the ensemble
nmin =   5;  % Number of points per leaf
k    =  10;  % Number of random cuts (it should be equal to the number of inputs)%



%% Calibrate and validate an ensemble of Extra-Trees

% Build an ensemble of Extra-Trees and return the predictions on the
% calibration dataset
[ensemble,finalResult_cal] = buildAnEnsemble(M,k,nmin,subset_cal);

% Run the ensemble on a validation dataset
[finalResult_val] = predictWithAnEnsemble(ensemble,subset_val);

% Calculate the model performance in calibration and validation
Rt2_fit(subset_cal(:,end),finalResult_cal)  % 
Rt2_fit(subset_val(:,end),finalResult_val)  % 

% Graphical analysis
figure;
subplot(221)
plot(subset_cal(:,end),'.-'); hold on; plot(finalResult_cal,'.-r'); grid on;
axis([1 length(subset_cal) min(subset_cal(:,end)) max(subset_cal(:,end))]);
xlabel('time'); ylabel('output');
legend('measured','predicted');
title('calibration - trajectory');
subplot(222)
plot(finalResult_cal,subset_cal(:,end),'.'); grid on
axis([min(subset_cal(:,end)) max(subset_cal(:,end)) min(subset_cal(:,end)) max(subset_cal(:,end))]); 
xlabel('measured'); ylabel('predicted');
title('calibration - scatter plot');
subplot(223)
plot(subset_val(:,end),'.-'); hold on; plot(finalResult_val,'.-r'); grid on;
axis([1 length(subset_val) min(subset_val(:,end)) max(subset_val(:,end))]);
xlabel('time'); ylabel('output');
legend('measured','predicted');
title('validation - trajectory');
subplot(224)
plot(finalResult_val,subset_val(:,end),'.'); grid on
axis([min(subset_val(:,end)) max(subset_val(:,end)) min(subset_val(:,end)) max(subset_val(:,end))]);
xlabel('measured'); ylabel('predicted');
title('validation - scatter plot');



%% k-fold cross-validation

% Define the parameters for the cross-validation
ns   = 5; % number of folds
flag = 1; % if flag == 1, an ensemble is built on the whole dataset at the end of the cross-validation. 
          % Otherwise (flag == 0), such model is not built.

% Shuffle the data
data_sh = shuffle_data(data);

% Run the cross-validation
[model] = crossvalidation_extra_tree_ensemble(data_sh,M,k,nmin,ns,flag);

% Model performance in calibration and validation
model.cross_validation.performance.Rt2_cal_pred_mean  % 
model.cross_validation.performance.Rt2_val_pred_mean  % 

% Graphical analysis
figure;
subplot(211)
plot(data_sh(:,end),'.-'); hold on; plot(model.complete_model.trajectories,'.-r'); grid on;
axis([1 length(data_sh) min(data_sh(:,end)) max(data_sh(:,end))]);
xlabel('time'); ylabel('output');
legend('measured','predicted');
title('calibration - trajectory');
subplot(212)
plot(model.complete_model.trajectories,data_sh(:,end),'.'); grid on
axis([min(data_sh(:,end)) max(data_sh(:,end)) min(data_sh(:,end)) max(data_sh(:,end))]); 
xlabel('measured'); ylabel('predicted');
title('calibration - scatter plot');



%% Input ranking

% Shuffle the data
data_sh = shuffle_data(data);

% Run the ranking algorithm
[result_rank] = input_ranking(data_sh,M,k,nmin);

% Results: Contribution (first column, %), Variable num. (secon column, [])
%
%    11.4821    2.0000
%    11.2524    1.0000
%    11.1244    4.0000
%    10.2580    5.0000
%     9.9643    6.0000
%     9.6793    3.0000
%     9.6423    7.0000
%     9.5031    9.0000
%     9.3559   10.0000
%     7.7382    8.0000



%% IIS algorithm

% Define the parameters
ns = 5;         % number of folds
p  = 5;         % number of SISO models evaluated at each iteration (this number must be smaller than the 
                % number of candidate inputs.
epsilon  = 0;   % tolerance
max_iter = 6;   % maximum number of iterations

% Shuffle the data
data_sh = shuffle_data(data);

% Run the IIS algorithm
[result_iis] = iterative_input_selection(data_sh,M,nmin,ns,p,epsilon,max_iter);

% Exit condition
result_iis.exit_condition
% --> An input variable was selected twice
result_iis.iters_done
% --> 6

% Selected variables (by iteration):
sel_variables    = nan(max_iter-1,1);
sel_variables(1) = result_iis.iter_1.best_SISO(1);
sel_variables(2) = result_iis.iter_2.best_SISO(1);
sel_variables(3) = result_iis.iter_3.best_SISO(1);
sel_variables(4) = result_iis.iter_4.best_SISO(1);
sel_variables(5) = result_iis.iter_5.best_SISO(1);

% Cumulated R2 of the MISO model
R2    = nan(max_iter-1,1);
R2(1) = result_iis.iter_1.MISO.cross_validation.performance.Rt2_val_pred_mean;
R2(2) = result_iis.iter_2.MISO.cross_validation.performance.Rt2_val_pred_mean;
R2(3) = result_iis.iter_3.MISO.cross_validation.performance.Rt2_val_pred_mean;
R2(4) = result_iis.iter_4.MISO.cross_validation.performance.Rt2_val_pred_mean;
R2(5) = result_iis.iter_5.MISO.cross_validation.performance.Rt2_val_pred_mean;
deltaR2 = [R2(1) ; diff(R2)];

% Graphical analysis
figure; 

subplot(321)
bar([1,2,3,4,5],deltaR2,'FaceColor','b'); hold on;
plot(R2,'o-','Color','k','LineWidth', 1, 'MarkerSize', 8, 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'k'); grid on;
axis([0.5 5.5 -0.2 1.0]);
set(gca,'XTick',[1,2,3,4,5]); 
set(gca,'XTickLabel', {num2str(sel_variables(1)), num2str(sel_variables(2)), num2str(sel_variables(3)), num2str(sel_variables(4)), num2str(sel_variables(5))},'Ylim',[0.00 1.0]);
xlabel('selected variables'); ylabel('R^2');
title('IIS');

subplot(322)
plot(data_sh(:,end),result_iis.iter_1.MISO.complete_model.trajectories,'.y','MarkerSize',12); grid on;
axis([min(data_sh(:,end)) max(data_sh(:,end)) min(data_sh(:,end)) max(data_sh(:,end))]);
xlabel('measured');  ylabel('predicted'); title('First iter');
legend('1 input');

subplot(323)
plot(data_sh(:,end),result_iis.iter_1.MISO.complete_model.trajectories,'.y','MarkerSize',12); hold on; grid on;
plot(data_sh(:,end),result_iis.iter_2.MISO.complete_model.trajectories,'.r','MarkerSize',12);
axis([min(data_sh(:,end)) max(data_sh(:,end)) min(data_sh(:,end)) max(data_sh(:,end))]);
xlabel('measured');  ylabel('predicted'); title('Second iter');
legend('1 input','2 inputs');

subplot(324)
plot(data_sh(:,end),result_iis.iter_1.MISO.complete_model.trajectories,'.y','MarkerSize',12); hold on; grid on;
plot(data_sh(:,end),result_iis.iter_2.MISO.complete_model.trajectories,'.r','MarkerSize',12);
plot(data_sh(:,end),result_iis.iter_3.MISO.complete_model.trajectories,'.b','MarkerSize',12);
axis([min(data_sh(:,end)) max(data_sh(:,end)) min(data_sh(:,end)) max(data_sh(:,end))]);
xlabel('measured');  ylabel('predicted'); title('Third iter');
legend('1 input','2 inputs','3 inputs');

subplot(325)
plot(data_sh(:,end),result_iis.iter_1.MISO.complete_model.trajectories,'.y','MarkerSize',12); hold on; grid on;
plot(data_sh(:,end),result_iis.iter_2.MISO.complete_model.trajectories,'.r','MarkerSize',12);
plot(data_sh(:,end),result_iis.iter_3.MISO.complete_model.trajectories,'.b','MarkerSize',12);
plot(data_sh(:,end),result_iis.iter_4.MISO.complete_model.trajectories,'.g','MarkerSize',12);
axis([min(data_sh(:,end)) max(data_sh(:,end)) min(data_sh(:,end)) max(data_sh(:,end))]);
xlabel('measured');  ylabel('predicted'); title('Fourth iter');
legend('1 input','2 inputs','3 inputs','4 inputs');

subplot(326)
plot(data_sh(:,end),result_iis.iter_1.MISO.complete_model.trajectories,'.y','MarkerSize',12); hold on; grid on;
plot(data_sh(:,end),result_iis.iter_2.MISO.complete_model.trajectories,'.r','MarkerSize',12);
plot(data_sh(:,end),result_iis.iter_3.MISO.complete_model.trajectories,'.b','MarkerSize',12);
plot(data_sh(:,end),result_iis.iter_4.MISO.complete_model.trajectories,'.g','MarkerSize',12);
plot(data_sh(:,end),result_iis.iter_5.MISO.complete_model.trajectories,'.c','MarkerSize',12);
axis([min(data_sh(:,end)) max(data_sh(:,end)) min(data_sh(:,end)) max(data_sh(:,end))]);
xlabel('measured');  ylabel('predicted'); title('Fifth iter');
legend('1 input','2 inputs','3 inputs','4 inputs','5 inputs');


% Results (by iteration):
%
% ITERATION:
%      1
% 
% Ranking:
%    11.4322    2.0000
%    11.2239    4.0000
%    11.0137    1.0000
%    10.9990    9.0000
%    10.0515    3.0000
%     9.9002    5.0000
%     9.2541    6.0000
%     9.1899    8.0000
%     9.0488   10.0000
%     7.8868    7.0000
% 
% Evaluating SISO models:
%     2.0000   -0.0133
%     4.0000    0.0982
%     1.0000    0.1128
%     9.0000   -0.2625
%     3.0000   -0.1791
% 
% Select variable:
%      1
% 
% Evaluating MISO model:
%     0.0967
% 
% ITERATION:
%      2
% 
% Ranking:
%    11.6567    2.0000
%    11.2476    4.0000
%    10.6002    1.0000
%    10.3026   10.0000
%     9.6088    6.0000
%     9.5620    5.0000
%     9.4690    3.0000
%     9.4298    7.0000
%     9.2929    9.0000
%     8.8304    8.0000
% 
% Evaluating SISO models:
%     2.0000   -0.0568
%     4.0000   -0.0833
%     1.0000   -0.7199
%    10.0000   -0.1401
%     6.0000   -0.3703
% 
% Select variable:
%      2
% 
% Evaluating MISO model:
%     0.4380
% 
% ITERATION:
%      3
% 
% Ranking:
%    11.2693   10.0000
%    10.6105    6.0000
%    10.4458    3.0000
%    10.2647    1.0000
%    10.1912    4.0000
%     9.9456    2.0000
%     9.8561    5.0000
%     9.6732    9.0000
%     9.0529    8.0000
%     8.6906    7.0000
% 
% Evaluating SISO models:
%    10.0000   -0.2185
%     6.0000   -0.2360
%     3.0000   -0.0595
%     1.0000   -0.4419
%     4.0000    0.2212
% 
% Select variable:
%      4
% 
% Evaluating MISO model:
%     0.6860
% 
% ITERATION:
%      4
% 
% Ranking:
%    11.8074    3.0000
%    10.8342    4.0000
%    10.7690    1.0000
%    10.3771    5.0000
%    10.2577    9.0000
%     9.7050   10.0000
%     9.3843    6.0000
%     9.3242    8.0000
%     8.8326    7.0000
%     8.7086    2.0000
% 
% Evaluating SISO models:
%     3.0000    0.1143
%     4.0000   -0.4346
%     1.0000   -0.2407
%     5.0000    0.1157
%     9.0000   -0.1934
% 
% Select variable:
%      5
% 
% Evaluating MISO model:
%     0.7432
% 
% ITERATION:
%      5
% 
% Ranking:
%    10.9356    3.0000
%    10.9120    1.0000
%    10.7994   10.0000
%    10.2602    4.0000
%     9.8262    6.0000
%     9.8162    9.0000
%     9.5558    2.0000
%     9.5396    5.0000
%     9.4433    7.0000
%     8.9117    8.0000
% 
% Evaluating SISO models:
%     3.0000    0.2874
%     1.0000   -0.2605
%    10.0000   -0.1437
%     4.0000   -0.4145
%     6.0000   -0.2599
% 
% Select variable:
%      3
% 
% Evaluating MISO model:
%     0.7774
% 
% ITERATION:
%      6
% 
% Ranking:
%    11.0108   10.0000
%    10.8553    3.0000
%    10.7338    5.0000
%    10.5315    1.0000
%    10.4461    4.0000
%     9.9265    9.0000
%     9.7200    8.0000
%     9.6861    6.0000
%     8.5453    2.0000
%     8.5446    7.0000
% 
% Evaluating SISO models:
%    10.0000   -0.2566
%     3.0000   -0.1651
%     5.0000    0.0771
%     1.0000   -0.3002
%     4.0000   -0.3576
% 
% Select variable:
%      5


%% Multiple runs of the IIS algorithm (with different shuffled datasets)

% Define the parameters
ns = 5;         % number of folds
p  = 5;         % number of SISO models evaluated at each iteration (this number must be smaller than the 
                % number of candidate inputs.
epsilon  = 0;   % tolerance
max_iter = 6;   % maximum number of iterations
                %
mult_runs = 10; % number of runs for the IIS algorithm               

% Shuffle the data
for i = 1:mult_runs
    eval(['data_sh_' num2str(i) '=' 'shuffle_data(data);']);
end

% Run the IIS algorithm
for i = 1:mult_runs    
    eval(['data_sh' '=' 'data_sh_' num2str(i) ';']);
    eval(['result_iis_' num2str(i) '=' 'iterative_input_selection(data_sh,M,nmin,ns,p,epsilon,max_iter);']);
    eval(['results_iis_n{i} = result_iis_',num2str(i),';']);
    clear data_sh    
end

% Plot the results
[X, R2] = visualize_inputSel(results_iis_n, size(data,2), mult_runs, max_iter );

% plot only the first max_iter variables
[X, R2] = visualize_inputSel(results_iis_n, max_iter, mult_runs, max_iter );
% change colormap
[X, R2] = visualize_inputSel(results_iis_n, max_iter, mult_runs, max_iter, 'Jet' );


% This code has been written by Stefano Galelli.







        
        

