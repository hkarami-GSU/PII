% <============================================================================>
% < Author: Gerardo Chowell  ==================================================>
% <============================================================================>

function dx=SEIRasymp1(t,x,params0,extra0)

beta0=params0(1);
beta1=params0(2);
k=params0(3);
rho=params0(4);
gamma1=params0(5);
N=params0(6);

dx=zeros(6,1);  % define the vector of the state derivatives: S, E, I, A, R, C

dx(1,1)= -(beta0*x(3,1)+beta1*x(4,1)).*x(1,1)./N; %S

dx(2,1)= (beta0*x(3,1)+beta1*x(4,1)).*x(1,1)./N - k*x(2,1); %E

dx(3,1)= k*rho*x(2,1) - gamma1*x(3,1); %I

dx(4,1)= k*(1-rho)*x(2,1) - gamma1*x(4,1); %A

dx(5,1)= gamma1*x(3,1) + gamma1*x(4,1); %R

dx(6,1)= k*rho*x(2,1); %C
