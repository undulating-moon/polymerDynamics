%% Intro
% In this scrip, we try to estimate error of a specific MSD.
%% Simulation
% Import data and previewing

figure
data_name='p47_5s'
load(['C:\Users\zhouquan\OneDrive\research\DNA随机游走模拟\数据模拟\MSDData4Quan_20210905\MSD_',data_name,'_All.mat'])
%% 
% cut down the declining part, we delete all data after the maximum one.


experiment_interval=5;
%% 
selecting_xyz=3;
Reserved_length=find(MSD_p47_5s_All.MSD_PP(:,selecting_xyz)==max(MSD_p47_5s_All.MSD_PP(:,selecting_xyz)));
interceptive_experiment_MSD=1e-18*MSD_p47_5s_All.MSD_PP(1:Reserved_length,selecting_xyz);

%% 
% Initialization

input=[3.0170e-09,0.2032,1.9977e-08,0.7321];;%[2.9304e-09,0.1929,2.0114e-08,0.6798];
m=input(1);
psi=input(2);
zeta=input(3);
H=input(4);

error=[];

for Total_experiment_number=[50,100,200]
    delta=0.1;


%% Main loop

    residual_history=[];
    N_simulation=100;
    
    for step_number=1:N_simulation
        disp(['------step',num2str(step_number),'--------'])
        residual_0=function_residual_calculating_fv2(Reserved_length,experiment_interval,interceptive_experiment_MSD,psi,zeta,H,m,delta,Total_experiment_number);
        residual_history=[residual_history,residual_0];
    end
%% Residual error calculation

    sumR=0
    sumR_square=0
    for i=1:length(residual_history)
        sumR=sumR+residual_history(i);
        sumR_square=sumR_square+residual_history(i)^2;
    end

    error=[error;((sumR_square-sumR^2/N_simulation)/(N_simulation*(N_simulation-1)))^0.5];
end