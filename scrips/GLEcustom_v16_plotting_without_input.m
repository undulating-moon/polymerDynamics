%% Intro
% In this scrip i try to plotting the figure with data directly get  from
% funtion fv1, instead of  importing saving data. 
% And this scrip changing all four parameters at a time.
%% Input parameters
% constants

kB=1.3806505e-23;
T=293;%confirmed
delta=0.1;
                                                                                                                                                                                                                                                                                               
Total_experiment_number=200;
storage_location='C:\Users\zhouquan\OneDrive\research\DNA随机游走模拟\阶段报告2022.04图像';
%% 

data_name='p77_5s';
experiment_interval=5;
multiple=experiment_interval/delta;
%load(['C:\Users\zhouquan\OneDrive\research\DNA随机游走模拟\数据模拟\MSDData4Quan_20210905\MSD_',data_name,'_All.mat'])

%plot(experiment_t,eval(['MSD_',data_name,'_All']).MSD_PP(:,1),'LineWidth',2,'displayname',' experiment data x')
%plot(experiment_t,eval(['MSD_',data_name,'_All']).MSD_PP(:,2),'LineWidth',2,'displayname',' experiment data y')
%plot(experiment_t,eval(['MSD_',data_name,'_All']).MSD_PP(:,3),'LineWidth',2,'displayname',' experiment data z')

load(['C:\Users\zhouquan\OneDrive\research\scrips\processed_data\svd smooth data ',data_name,'.mat']);
%load(['C:\Users\zhouquan\OneDrive\research\scrips\processed_data\expomential fitness data ',data_name,'.mat']);
experiment_t=multiple:multiple:multiple*length(expomential_fitness_data);
figure
%plot(experiment_t,expomential_fitness_data,'LineWidth',2,'displayname','expomential fitness data')
hold on
plot(experiment_t,smoothed_one_curve,'LineWidth',2,'displayname','svd smoothed data')
Reserved_length=length(expomential_fitness_data);
interceptive_experiment_MSD=1e-18*expomential_fitness_data;
Tmax=Reserved_length*experiment_interval;   

%% 
% variable parameters
n_parameter_sets=1;
input=zeros(n_parameter_sets,4);
%legend_name=zeros(n_parameter_sets,1);
input(1,:)=[6.0727E-9,0.5127,2.0212E-8,0.69714];
%input(2,:)=[2.38233593716721e-09,0.0536106688898815,1.99569264498908e-08,0.753456070825410];

sample_m=[input(:,1)];
sample_psi=[input(:,2)];
sample_zeta=[input(:,3)];
sample_H=[input(:,4)];

%% Main loop

for i=1:length(sample_m)
    m=sample_m(i);
    psi=sample_psi(i);
    zeta=sample_zeta(i);
    H=sample_H(i);
    %simulation_MSD=function_GLEsubdiffusion_fv1(m,kB,T,psi,zeta,H,delta,Tmax,Total_experiment_number,storage_location,false,false);
    simulation_MSD=function_GLEsubdiffusion_fv4_less_memory(m,kB,T,psi,zeta,H,delta,Tmax,Total_experiment_number,1);
    legend_name=strcat('simulation MSD',' m= ',num2str(m),' \psi=',num2str(psi),' \zeta= ',num2str(zeta),' H= ',num2str(H));
    plot(simulation_MSD*1E18,"LineWidth",2,'displayname',legend_name);
    hold on
end

%% 
%x:1, y:2, z:3
selecting_xyz=1;
if selecting_xyz==1
    selecting_curve='x';
elseif selecting_xyz==2
    selecting_curve='y';
elseif selecting_xyz==3
    selecting_curve='z';
end
%Reserved_length=find(eval(['MSD_',data_name,'_All']).MSD_PP(:,selecting_xyz)==max(eval(['MSD_',data_name,'_All']).MSD_PP(:,selecting_xyz)));
%interceptive_experiment_MSD=1e-18*eval(['MSD_',data_name,'_All']).MSD_PP(1:Reserved_length,selecting_xyz);


residual=sum((simulation_MSD(experiment_interval/delta*(1:Reserved_length),1)-interceptive_experiment_MSD').^2);
%% 
xlabel(['t / ',num2str(delta),' s'],'FontSize',14);
ylabel('MSD / nm^2','FontSize',14);
title({['processed experiment MSD & simulation MSD -t'],['processed experiment data: ',strrep(data_name,'_','\_'),'\_',selecting_curve],[' \Deltat= ',num2str(delta),' N= ',num2str(Total_experiment_number),' residual = ',num2str(residual),'m^4']},'Color','r','FontSize',10);
legend('Location','southoutside');
%% saving

saveas(gcf,strrep([storage_location,'\image\test.MSD','.m=',num2str(m),'.psi=',num2str(psi),'.zeta=',num2str(zeta),'.H=',num2str(H),'.delta=',num2str(delta),'.N=',num2str(Total_experiment_number),'.jpg'],'       ',' '))

