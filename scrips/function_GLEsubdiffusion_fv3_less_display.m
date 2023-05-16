function [RMSX]=function_GLEsubdiffusion_fv3_less_display(m,kB,T,psi,zeta,H,delta,Tmax,Total_experiment_number,explosion_limit)
%%FUNCTION_GLESUBDIFFUSION_FV2_EXPLOSION_ALARMING Introduction
% This scrip is the one serves for formal version of scrip. We delete part of 
% unnecessary fuctions, like saving fuction and plotting.
% 
% But adding new function like explosion judgement.
%% Defining fuction
%% 
% control image saving
explosion_judgement=false;
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
        a(n,i+1)=-(zeta/m)*2*H*(2*H-1)*delta^(2*H-2)*((KH.^(2*H-2))*v(n,1:i)')*delta-psi*x(n,i)+(2*zeta*kB*T)^0.5*(BH(i+1)-BH(i))/(m*delta);%ODE
        if 10*i/imax==floor(10*i/imax)
            if max(x(n,i-imax/10+1:i))>=explosion_limit
                explosion_judgement=true;
                disp('explosion!')
                break
            end
        end
    end
%% 
% monitoring the progress rate
    if explosion_judgement==true
        break
    end
    if max(x(n,:))>=explosion_limit
        explosion_judgement=true;
        disp('explosion!!!!!!')
        break
    end
end
%% Generating
if explosion_judgement==true
    RMSX=true;
else
    averageX=zeros(imax+1,1);
    RMSX=zeros(imax+1,1);
    for i=1:imax+1
        averageX(i,1)=sum(x(:,i))/Total_experiment_number;
        RMSX(i,1)=sum(x(:,i).^2)/Total_experiment_number-averageX(i,1)^2;
    end
end