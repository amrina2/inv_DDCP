%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Computing Artifact
%%%%% Author: Amrina Ferdous
%%%%% File name: frontier.m 
%%%%% Purpose: This is the file that researchers need to run only. All
%%%%% other files in this folder are the supporting files due to frontier.m
%%%%% file.
%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% An example: This is the file where researchers need to formulate the
%%%%% problem which is in our case: calculating the frontier of the
%%%%% \cite{latexTarantola} paper's geophysics ill-posed nonlinear inverse 
%%%%% problem. This is the main background of the coding. 
%%%%%%%%%%%%
%% Setup parameters (for example)
% Defining the integration grid w
m= 100; % Number of points of integration grid
wmin=0; % [km] Lower limit of the integration domain
wmax=100; % [km] Upper limit of the integration domain
w = linspace(wmin,wmax,m); % w is discreatized with m points and should be in a row. 
% Loading data 
d= load("data.txt"); %[km] data vector (column vector)
n = length(d); %Total number of data 
x = linspace(wmin,wmax,n)' ; % Defining the data domain. x are in n number of data points and it should be in a column.
sigma_d= sqrt(0.001) ; % [km] Constant standard deviation of the data error.
sigma = 5 ; % [km] Prior uncertainty of $z_{0}(w)$ in units.
theta = 1; % [km] Spread in the prior covariance. 
K=10; % Number of iteration 
%%%%
%% Main routine
% Calling zhat function i.e., recovering the frontier using the inv_DDCP.m file. For more
% details about zhat, please see the inv_DDCP.m file and read the Tarantola paper. 
zhat=inv_DDCP(w,x,d,sigma_d,sigma,theta,K,@Gfun,@ffun)
%%%%
%%%%%%%%%%%%%%%%%%%
%% Display results
%%%%%%%%%%%%%%%%%%%
%% Plot of the frontier i.e., zhat
figure(8); clf; hold on;
for i=1:K
    plot(w,zhat(:,i),'LineWidth', 5) %setting the line width of a plot
    set(gca, 'Fontsize', 14) %setting a title of the plot
    set(gca,'FontWeight','bold') %%setting a font weight of the plot
    title('Line plot of frontier') %setting a title of the plot
    xlabel('w') %labeling x-axis
    ylabel('z(w)') %labeling x-axis
end
hold off
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Validation of our frontier results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Creating true z(w)from which the data in data.txt were calculated
maxz=2.5; % [km] 
a=int8(2*m/5);
b=int8(m/2);
c=int8(3*m/5);
ztrue=zeros(1,m);
% Now considering z(w) as continuous. Specifically,
ztrue(a:c)=maxz*exp(-5*(w(a:c)-w(b)).^2 /(m));
%%
%Creating true u (i.e., utrue) and predicted u (i.e., uhat)
N=100; %(say)
x_N = linspace(wmin,wmax,N)'; % xdata has N data points
utrue = g_small_fun(x_N,w,ztrue);
% predicted u = u2 = considering g(zhat)
uhat = g_small_fun(x_N,w,zhat(:,end)');
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%
%% Display Validation
%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Checking accuracy
figure(9)
subplot(2,1,1)
%%%% Plot to validate our estimated zhat for frontier
%plot of true z(w) vs estimated z(w) i.e., zhat
plot(w,ztrue',w,zhat(:,end),'LineWidth', 5)
set(gca, 'Fontsize', 14) %setting a title of the plot
set(gca,'FontWeight','bold') %%setting a font weight of the plot
title('Plot of ztrue vs zhat') %setting a title of the plot
xlabel('w')
legend('ztrue','zhat')
subplot(2,1,2)
%plot of utrue vs uhat 
plot(x_N,utrue,x_N,uhat,'LineWidth', 5)
title('Plot of utrue vs uhat')
xlabel('x')
set(gca, 'Fontsize', 14)
set(gca,'FontWeight','bold') 
legend('utrue','uhat')
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Display Relative Errors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(10)
subplot(2,1,1)
%%%% Relative error in parameter estimate, z
plot(w,abs(ztrue' - zhat(:,end))/norm(ztrue'),'LineWidth', 5)
set(gca, 'Fontsize', 14) 
set(gca,'FontWeight','bold') 
title('Relative error in parameter estimate, z') 
xlabel('w')
set(gca, 'Fontsize', 14)
%%%%%%%%
subplot(2,1,2)
%%%%%% Relative error in prediction
plot(x_N,abs(uhat - utrue)/norm(utrue),'LineWidth', 5)
title('Relative error in prediction') 
xlabel('x')
set(gca, 'Fontsize', 14)
set(gca,'FontWeight','bold') 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% THE END 