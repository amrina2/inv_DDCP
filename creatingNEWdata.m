%%
% Creating our own forward problem i.e. u(x). Also z(w). 
% x are in n data point & it should be in a colums
% w is discreatize with m point & should be in a row
H=10;
N=100;
m=100;
W=100;
x_N = linspace(0,W,N)'; % xdata has N data points
w = linspace(0,W,m);
%%
% Creating our own z(w)
maxz=2.5; % According to the Tarantola paper, maxz=2.5 km. 
a=int8(2*m/5);
b=int8(m/2);
c=int8(3*m/5);
ztrue=zeros(1,m);
% Now considering z(w) as continuous.Specifically,z(w) as exponential
ztrue(a:c)=maxz*exp(-5*(w(a:c)-w(b)).^2 /(m) );
u = g_small_fun(x_N,w,ztrue);
%%
sigma_d=sqrt(0.001); %%% we took 0.01 for the other works.This data considers 0.001
n1=14; %(say)
n2=int8(N/n1);
d1=u(1:n2:N,1);
n3=length(d1);
d1= d1+sigma_d*rand(n3,1);
% Making all positive d1 
d1(d1<0)=-d1(d1<0);
x1=x_N(1:n2:N,1);
plot(x1,d1,'*'); 
%%
% Saving our new data set
dlmwrite('dataNEW.txt','d1')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%