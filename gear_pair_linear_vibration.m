clear
% this code is for to linear vibration response of gear system

z0=zeros(2,1); % initial conditipon
z0(1)=0 ;z0(2)=0;
omg=2.5;
N1=200; % number of frequency points
N=20; % segment numbers in one cycle
NN=30; % cycle number
NN1=N*NN+1; %total number of points
%x=zeros(NN1,2);
ebs=0.25; % epsala small parameter
miudamp=0.05; %
f0=0.1; % the nondimensional static force
f1=0.15; % the nondimensional dynamic excitation
T=2*pi/1; % nondimension time period
dt=T/N; % time step
t= [0 128.5]; % discretized time
[t,z]=ode45(@fun3,t,z0,[],omg,ebs,miudamp,f0,f1);
figure (1), plot(t,z(:,1));grid; % plot the displacement
figure (2),plot(t,z(:,2));grid; % plot the velocity