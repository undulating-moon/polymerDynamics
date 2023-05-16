%% GLEcustom_v13 Intro
% This scrip is about to research about the meating time of two particular genetic 
% locus. And in this scrip I will use a new kind of random force. In case of the 
% unreliable result of v9 (1 dimension scrip), in this scrip we use 3 dimension. 
% And i change the intensity of fbm.
% 
% And in order to expedite the scrip, i ignore part of the early history in 
% GLE equation.
%% v0 Intro
% I try to adapt the v13 scrip, for there are several problems, for example, 
% the distance when t=0, and delete the "ignored" function.
%% v1 Intro
% For there is hard to meet, i try to simulate 1 dimension FPT, get one dimensional 
% FPT distribution. And then calculate the 3 dimension FPT.
%% v2 Intro
% I try to discard the inetia item, and simulate the FPT distribution. Parameter 
% using the fitness result of GLENII.
%% Data Importing
% import parameters

a1=[4.6565e-08,2.0355e-07,0.78971];%GLENII p47 parameter
mpsi=a1(1);%psi=0---> no potential
zeta=a1(2);
H=a1(3);
kB=1.3806505e-23;
T=293;
meeting_distance=5e-8;
%% 
% import Time interval

delta=0.1;
%% 
% import Total Time

Tmax=500;
imax=Tmax/delta;
%% 
% import Total_experiment_number

Total_experiment_number=10000;
%% 
% import distribution intervals

distribution_interval=10;
%% Random Force
% We create the random force in this part and smooth it.
% 
% The interval of each point is delta

%BH=wfbm(H,imax+1);
%X=1:imax;
%x_length = linspace(1,Tmax);
%BHt=interp1(X,BH(2:imax+1),x_length,'pchip');%smooth the random
%plot(x_length,BHt)%draw the figure to verify
%% Initialization
% initialization v, x, and a data type

v=zeros(1,imax+1);
x=zeros(1,imax+1);
BH=zeros(1,imax+1);%may be imax+2 
first_passage_step=zeros(Total_experiment_number,1);
%% Main loop
% initialization v, x, and a when t =0

for n=1:Total_experiment_number
    v(1,1)=0; %normrnd(0,(kB*T/m)^0.5);
    % 这里理论上是要根据能均分定理给出一个随机的初始速度，但是由于这样会引入额外参数，这里就暂时定为0
    x(1,1)=normrnd(0,(kB*T/(mpsi))^0.5);%should be random
    BH(1,:)=(delta^H)*wfbm(H,imax+1);
    for i=1:imax
        KH=i:-1:1;
        x(1,i+1)=(-zeta*2*H*(2*H-1)*delta^(2*H-2)*((KH.^(2*H-2))*v(1:i)')*delta+(2*zeta*kB*T)^0.5*(BH(i+1)-BH(i))/delta)/mpsi;%ODE GLENII
        v(1,i+1)=(x(1,i+1)-x(1,i))/delta;
        if abs(x(1,i+1))<=meeting_distance
            break
        end
    end
%% 
% monitoring the progress rate

    progress_rate=n/Total_experiment_number*100;
    first_passage_step(n,1)=i;%means betwean t=(i-1)*delta adn t=i*delta
    if progress_rate==floor(progress_rate)
        disp(['progress rate  ',num2str(progress_rate),' %'])
    end
end
%% Data output

FPT_distripution=zeros(floor(max(first_passage_step)*delta/distribution_interval)+1,1);
for i=1:Total_experiment_number
    FPT_distripution(floor(first_passage_step(i)*delta/distribution_interval)+1,1)=FPT_distripution(floor(first_passage_step(i)*delta/distribution_interval)+1,1)+1;
end
%% Plotting and saving

storage_location='C:\Users\zhouquan\OneDrive\research\DNA随机游走模拟\阶段报告2022.04图像';
save([storage_location,'\data\FPT of GLENII','.m=',num2str(m),'.psi=',num2str(psi),'.zeta=',num2str(zeta),'.H=',num2str(H),'.d=',num2str(meeting_distance),'.N=',num2str(Total_experiment_number),'.mat'],"first_passage_step","FPT_distripution");
%%
figure
bar(1:floor(max(first_passage_step)*delta/distribution_interval)+1,FPT_distripution(:,1));
set(gca,'yscale','log');
xlabel(['passage time / ',num2str(distribution_interval),'s'],'FontSize',14);
ylabel('number of particle meet','FontSize',14);
title({'passage time distribution of GLENII',[' \Deltat=',num2str(delta),' N=',num2str(Total_experiment_number),' d=',num2str(meeting_distance)],['[m*\psi,\zeta,H] = [',num2str(mpsi),', ',num2str(zeta),', ',num2str(H),']']},'Color','r','FontSize',16)
saveas(gcf,[storage_location,'\image\FPT_distribution of GLENII','.m=',num2str(m),'.psi=',num2str(psi),'.zeta=',num2str(zeta),'.H=',num2str(H),'.d=',num2str(meeting_distance),'.N=',num2str(Total_experiment_number),'.jpg']);