function [sqrN_error]=function_error_estimation_fv0_no_output(Reserved_length,experiment_interval,interceptive_experiment_MSD,psi,zeta,H,m,delta)
%% FUNCTION_ERROR_ESTIMATION_FV0_NO_OUTPUT Intro
% This scrip is used to turn part of main scrip (the part calculating relationship 
% between error and N) into a function.
%% Defining function
%% 
% 
error_test_N=100;
error_residual_history=[];
N_simulation=50;
for error_step_number=1:N_simulation
    if error_step_number/10==floor(error_step_number/10)
        disp(['error_N relationship progress rate',num2str(floor(error_step_number)*2),'%'])
    end
    residual_0=function_residual_calculating_fv3_less_display(Reserved_length,experiment_interval,interceptive_experiment_MSD,psi,zeta,H,m,delta,error_test_N);
    error_residual_history=[error_residual_history,residual_0];
end
sumR=0;
sumR_square=0;
for i=1:length(error_residual_history)
    sumR=sumR+error_residual_history(i);
    sumR_square=sumR_square+error_residual_history(i)^2;
end
error=((sumR_square-sumR^2/N_simulation)/(N_simulation-1))^0.5;
sqrN_error=(error_test_N)^0.5*error*(N_simulation)*0.5;