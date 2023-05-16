%% 
syms y(x);%请在一开始确定x为方程变量而非普通字符或自变量
dsolve(y+(x^2-4*x)*diff(y,x)==0);
%%
FBM = wfbm(0.3,101,'plot');
%%
A = 1;
B = 2;
tspan = [0 5];
y0 = [0 0.01];
[t,y] = ode45(@(t,y) odefcn(t,y,A,B), tspan, y0);

plot(t,y(:,1),'-o',t,y(:,2),'-.')

%%
a=[0:7]
b=a(1)
%%
i=1
loopnumber=1
v=zeros(2,2)
v(loopnumber,i)
%%
function dydt = odefcn(t,y,A,B)
dydt = zeros(2,1);
dydt(1) = y(2);
dydt(2) = (A/B)*t.*y(1);
end
