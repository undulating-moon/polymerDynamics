%% Intro
% In this scrip i try to plotting the figure with data directly get  from
% funtion fv1, instead of  importing saving data.
% And this scrip respectively changing 4 parameter and plot figure.
%% Input parameters
% constants

kB=1.3806505e-23;
T=293;%confirmed
delta=0.1;
Tmax=500;
Total_experiment_number=2000;
storage_location='C:\Users\zhouquan\OneDrive\research\DNA随机游走模拟\阶段报告2022.04图像';
%% 
% variable parameters

sample_m=5e-9;%4e-9:1e-9:9e-9
sample_psi=[0.1];%center 0.5
sample_zeta=2e-8:5e-9:4e-8;%center 3e-8
sample_H=[0.75];%center 0.75
changing_parameter=sample_zeta;
parameter_letter='\zeta=';
%% Main loop

figure
for i_m=1:length(sample_m)
    m=sample_m(i_m);
    for i_psi=1:length(sample_psi)
        psi=sample_psi(i_psi);
        for i_zeta=1:length(sample_zeta)
            zeta=sample_zeta(i_zeta);
            for i_H=1:length(sample_H)
                H=sample_H(i_H);
                MSDX=function_GLEsubdiffusion_fv1(m,kB,T,psi,zeta,H,delta,Tmax,Total_experiment_number,storage_location,false,false);
                plot(MSDX*1E18,"LineWidth",2,'displayname',[parameter_letter,num2str(sample_m(i_m))]);
                hold on
            end
        end
    end
end
%% 

xlabel(['t / ',num2str(delta),' s'],'FontSize',14);
ylabel('MSD / nm^2','FontSize',14);
title({'MSD-t',strrep([' m= ',num2str(sample_m),' \psi=',num2str(sample_psi),' \zeta= ',num2str(sample_zeta),' H= ',num2str(sample_H)],'       ',' '),[' \Deltat= ',num2str(delta),' N= ',num2str(Total_experiment_number)]},'Color','r','FontSize',10);
legend;
%% saving

saveas(gcf,strrep([storage_location,'\image\MSD','.m=',num2str(sample_m),'.psi=',num2str(sample_psi),'.zeta=',num2str(sample_zeta),'.H=',num2str(sample_H),'.delta=',num2str(delta),'.N=',num2str(Total_experiment_number),'.jpg'],'       ',' '))

