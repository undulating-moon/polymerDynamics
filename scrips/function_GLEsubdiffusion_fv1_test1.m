function [MSDX]=function_GLEsubdiffusion_fv1_test1(m,kB,T,psi,zeta,H,delta,Tmax,Total_experiment_number,storage_location,plotting_control,save_image_control)
%% FUNCTION_GLESUBDIFFUSION_FV1 Introduction
% This scrip is the one serves for formal version of scrip. We delete part of 
% unnecessary fuctions.
% 
% fv1 add the figure saving function. You can turn plotting and saving function 
% on by putting true in last two parameters of function. 

%fv1_test1 test the reliability of scrip
%% Defining fuction
%% 
% control image saving
%plotting_control=true
%save_image_control=true;
%% Main loop
% initialization
imax=Tmax/delta;
v=zeros(Total_experiment_number,imax+1);
x=zeros(Total_experiment_number,imax+1);
a=zeros(Total_experiment_number,imax+1);
%% 
% loop
for n=1:Total_experiment_number
    v(n,1)=normrnd(0,(kB*T/m)*0.5);%should be random;
    x(n,1)=0;%no need to be random
    a(n,1)=0;%no need to be random
    BH=(delta^H)*wfbm(H,imax+1);
    for i=1:imax
        KH=i:-1:1;
        x(n,i+1)=x(n,i)+v(n,i)*delta;
        v(n,i+1)=v(n,i)+a(n,i)*delta;
        a(n,i+1)=-(zeta/m)*v(n,i)-psi*x(n,i)+(2*zeta*kB*T)^0.5*(BH(i+1)-BH(i))/(m*delta);%ODE
    end
%% 
% monitoring the progress rate
    progress_rate=n/Total_experiment_number*100;
    if progress_rate/10==floor(progress_rate/10)
        disp(['progress rate  ',num2str(progress_rate),' %'])
    end
end
%% Generating
averageX=zeros(imax+1,1);
MSDX=zeros(imax+1,1);
for i=1:imax+1
    averageX(i,1)=sum(x(:,i))/Total_experiment_number;
    MSDX(i,1)=sum(x(:,i).^2)/Total_experiment_number-averageX(i,1)^2;
end
%% Plotting
if plotting_control==true
    figure
    plot(MSDX*1E18,"LineWidth",2);
    xlabel(['t / ',num2str(delta),' s'],'FontSize',14)
    ylabel('MSD / nm^2','FontSize',14);
    title({'MSD-t',[' H= ',num2str(H),' m= ',num2str(m),' \zeta= ',num2str(zeta),' \psi=',num2str(psi)],[' \Deltat= ',num2str(delta),' N= ',num2str(Total_experiment_number)]},'Color','r','FontSize',10);
else
    save_image_control=false;
end
if save_image_control==true
    saveas(gcf,[storage_location,'\image\MSD_N.P.','H=',num2str(H),'.m=',num2str(m),'.zeta=',num2str(zeta), ...
        '.psi=',num2str(psi),'.delta=',num2str(delta),'.N=',num2str(Total_experiment_number),'.jpg'])
end