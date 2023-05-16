%% Introduction
% In this scrip, we just try to test the calcualtion of residual£¨²Ð²î£©. This version  
% uses 3 parameters: m, zeta and psi. The ranges of parameters come from fv3.1_3p_grid_search.
% 
% With the efficient of grid searching, its hard to get to precis parameter. 
% We use *Gradient descent *to get the precise parameter values.
%% Import data and previewing

figure
data_name='p47_5s'
load(['C:\Users\dell\Desktop\zq\MSDData4Quan_20210905\MSD_',data_name,'_All.mat'])
plot(MSD_p47_5s_All.MSD_PP,'LineWidth',2)
legend('x','y','z','Location',"best");
%% 
% cut down the declining part, we delete all data after the maximum one.

Reserved_lengths=zeros(3,1);
experiment_interval=5
%% 
% X part

%experiment_MSDX=eval(['MSD_',data_name,'_All']).MSD_PP(:,1);
%max_experiment_MSDX=max(experiment_MSDX);
%Reserved_lengths(1,1)=find(MSD_p47_5s_All.MSD_PP(:,1)==max(MSD_p47_5s_All.MSD_PP(:,1)));
%interceptive_experiment_MSDX=1e-18*MSD_p47_5s_All.MSD_PP(1:Reserved_lengths(1,1),1);
%% 
% Y part

%Reserved_lengths(2,1)=find(MSD_p47_5s_All.MSD_PP(:,2)==max(MSD_p47_5s_All.MSD_PP(:,2)));
%interceptive_experiment_MSDY=1e-18*MSD_p47_5s_All.MSD_PP(1:Reserved_lengths(2,1),2);
%% 
% Z part

%Reserved_lengths(3,1)=find(MSD_p47_5s_All.MSD_PP(:,3)==max(MSD_p47_5s_All.MSD_PP(:,3)))
%interceptive_experiment_MSDZ=1e-18*MSD_p47_5s_All.MSD_PP(1:Reserved_lengths,3);
%% Data Importing
% import Time interval

delta=0.1;
%% 
% import Total_experiment_number

Total_experiment_number=100;
%% 
% data import

Reserved_length=Reserved_lengths(3,1);
interceptive_experiment_MSD=interceptive_experiment_MSDZ;
%% Initilizing parameters
% input the initial parameter

input=[6E-9	,0.5,	2.0E-8	,0.75];
m_0=input(1);
psi_0=input(2);
zeta_0=input(3);
H_0=input(4);


%% 
% put it in use

m=m_0
psi=psi_0
zeta=zeta_0
H=H_0
%% 
% decide step length and sample length (proportional)

step_length=1/50;
sample_length=1/100;
%% 
% add parameter history and gradient history

parmeter_history=[m,psi,zeta,H];
residual_delta=[];
%% Calculating decline in three direction
% 

for step_number=1:5
    disp(['------step',num2str(step_number),'--------'])
    residual_0=function_residual_calculating_fv2(Reserved_length,experiment_interval,interceptive_experiment_MSD,psi,zeta,H,m,delta,Total_experiment_number);
    residual_m=function_residual_calculating_fv2(Reserved_length,experiment_interval,interceptive_experiment_MSD,psi,zeta,H,m*(1+sample_length),delta,Total_experiment_number);
    residual_psi=function_residual_calculating_fv2(Reserved_length,experiment_interval,interceptive_experiment_MSD,psi*(1+sample_length),zeta,H,m,delta,Total_experiment_number);
    residual_zeta=function_residual_calculating_fv2(Reserved_length,experiment_interval,interceptive_experiment_MSD,psi,zeta*(1+sample_length),H,m,delta,Total_experiment_number);
    residual_H=function_residual_calculating_fv2(Reserved_length,experiment_interval,interceptive_experiment_MSD,psi,zeta,H*(1+sample_length),m,delta,Total_experiment_number);
%% 
% saving residual
    residual_delta=[residual_delta;[residual_0,residual_m,residual_psi,residual_zeta,residual_H]];
end
%% saving
%%
saving_location='C:\Users\dell\Desktop\zq\Êý¾ÝÄ£Äâ\P.S.G.D'
save([saving_location,'\data\history data N=',num2str(Total_experiment_number),' step=',num2str(step_length),datestr(now,' yyyy-mm-dd HH'),'.mat'],'residual_change')

%% Plot and plot saving: residual_change
%%figure
for i=1:5
    plot((residual_delta(:,i))*1e36,'k-o','MarkerFaceColor','k')
    hold on
end