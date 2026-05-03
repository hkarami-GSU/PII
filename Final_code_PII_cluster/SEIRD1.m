% <============================================================================>
% < Author: Gerardo Chowell  ==================================================>
% <============================================================================>

function dx=SEIR1(t,x,params0,extra0)

beta0=params0(1);
k=params0(2);
rho=params0(3);
gamma1=params0(4);
N=params0(5);

dx=zeros(5,1);  % define the vector of the state derivatives: S, E, I, R, D

dx(1,1)= -beta0*x(1,1).*x(3,1)./N; %S

dx(2,1)= beta0*x(1,1).*x(3,1)./N - k*x(2,1); %E

dx(3,1)= k*x(2,1) - gamma1*x(3,1); %I

dx(4,1)= (1-rho)*gamma1*x(3,1); %R

dx(5,1)= gamma1*rho*x(3,1); %D
