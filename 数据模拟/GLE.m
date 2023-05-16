%%
%参数设定
m=1;%这组数据仅仅用作测试
psi=1;
zeta=1;
H=0.75;
kB=1;
T=1;
%%
%初始值设定
tspan = [0 5];
y0 = [0 0.01];
%%
%方案1：使用拆分法
[t,y] = ode45(@(t,y) odefcn(t,y,m,psi,zeta,kB,T,KH,H), tspan, y0);

function dydt = odefcn(t,y,m,psi,zeta,kB,T,H)
dydt = zeros(4,3,2,1);
dydt(1) = y(2);
dydt(4)=(m*psi*y(1)+zeta*y(3))/((2*zeta*kB*T)^0.5);
dydt(3) =2*H*(2*H-1)*(t-u)^(2H-1)*y(2);
end

%y(1)=x,y(2)=v,y(3)=P,y(4)=BH

