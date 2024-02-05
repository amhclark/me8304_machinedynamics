clear
% this code is for to linear csee Fig4.3 (a) on p176
% of A Kahraman
% Non-linear Linear Dynamic Analysis of Geared System by Ahmet Kahraman
% Figure 4.3 (a) n p176, Equation 4.6 on p173
z0=zeros(2,1); % initial conditipon
z0(1)=0 ;z0(2)=0;
omgu=2.5;omgl=0.2;
N1=200; % number of frequency points
step=(omgu-omgl)/N1;
N=20; % segment numbers in one cycle
NN=30; % cycle number
NN1=N*NN+1; %total number of points
%x=zeros(NN1,2);
ebs=0.25; % epsala small parameter

miudamp=0.05; %
f0=0.1; % Fm in the original equation
fah=0.05; % fah in the original equation
T=2*pi/1;
%T=(2*pi)/omg;
dt=T/N;
ts=0:dt:NN*T;
a=zeros(N1,1);w=zeros(N1,1);
for ii=1:N1
    omg=omgl+ii*step; % excitation frequency
    %f1=fah*omg*omg;
    f1=0.15;
    w(ii)=omg;
    [t,z]=ode45(@fun3,ts,z0,[],omg,ebs,miudamp,f0,f1);
    x=z(NN1-2*NN:NN1,1);
    xu=max(x);xl=min(x);
    a(ii)=(xu-xl)/2;
end

plot(w,a);grid;
%plot(t(a0:a1)/T,x(a0:a1,1));grid on;
%subplot(2,1,1);plot(t/T,x(:,1));grid on;
%subplot(2,1,2);plot(t/T,x(:,2));grid on;
%plot(x(1000:10:50001,1),x(1000:10:50001,2),'.');
% two different solution are ford with two initial conditions
%X0=[4.2010 -0.9602]; [-3.3068 -4.9625];