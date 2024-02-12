% Adapted from provided code written by Dr. J. Yang 
% Machine Dynamics - Memorial University

% This code is to analyze forward dynamics of a four bar linkage:

clear
close all
clc

load fourbardate.mat                       % load kinematic dataset    

%% Input Data

% CGs
rG2 = 1.0; delt2 = 0;                           %bar 2
rG3 = 2.5; delt3 = 30*pi/180;                   %bar 3
rG4 = 1.5; delt4 = 0;                           %bar 4

% external load for bar 3
RP3 = 5*exp((30*pi/180)*(1j));                  % Force point bar 3
FP3 = 12*exp(270*pi/180*1j);                    % force (lbs)
T3 = -20;                                       % torque (lb*in)
% FP3=0;T3=0;

% external load for bar 4
RP4 = 5*exp((90*pi/180)*(1j));                  % Force point bar 4
FP4 = 60*exp(-45*pi/180*1j);                    % force (lbs)
T4 = 25;                                        % torque (lb*in)
%P4=0;T4=25;

m2 = 0.002*r2; m3 = 0.030*r3; m4 = 0.010*r4;    % mass
I2 = 0.004*r2; I3 = 0.060*r3; I4 = 0.020*r4;    % inertia


%% Matrix Init

A = zeros(9,9); b = zeros(9,1);  
X = zeros(9,n);

A(1,1) = 1; A(1,3) = 1;
A(2,2) = 1; A(2,4) = 1;
A(3,9) = 1;
A(4,3) = -1; A(4,5) = 1;
A(5,4) = -1; A(5,6) = 1;
A(7,5) = -1; A(7,7) = 1;
A(8,6) = -1; A(8,8) = 1;


%% Calculations

for i1 = 1:n
    % acceleration of CG
    aG2 = rG2*alph2*(1j)*exp((theta2(i1)+delt2)*(1j))-rG2*omega2^2*exp((theta2(i1)+delt2)*(1j));
    aG4 = rG4*alph4(i1)*(1j)*exp((theta4(i1)+delt4)*(1j))-rG4*omega4(i1)^2*exp((theta4(i1)+delt4)*(1j));
    aA = r2*alph2*(1j)*exp((theta2(i1)*(1j)))-r2*omega2^2*exp((theta2(i1)*(1j)));
    aPA = rG3*alph3(i1)*(1j)*exp((theta3(i1)+delt3)*(1j))-rG3*omega3(i1)^2*exp((theta3(i1)+delt3)*(1j));
    aG3 = aA+aPA;

    % vectors for each CG
    R12 = -rG2*exp((theta2(i1)+delt2)*(1j));
    R32 = R12+r2*exp(theta2(i1)*(1j));
    R23 = -rG3*exp((theta3(i1)+delt3)*(1j));
    R43 = R23+r3*exp(theta3(i1)*(1j));
    R14 = -rG4*exp((theta4(i1)+delt4)*(1j));
    R34 = R14+r4*exp(theta4(i1)*(1j));

    A(3,1) = -imag(R12);A(3,2)=real(R12);A(3,3)=-imag(R32);A(3,4)=real(R32);
    A(6,3) = imag(R23);A(6,4)=-real(R23);A(6,5)=-imag(R43);A(6,6)=real(R43);
    A(9,5) = imag(R34);A(9,6)=-real(R34);A(9,7)=-imag(R14);A(9,8)=real(R14);
    
    b(1) = m2*real(aG2); b(2)=m2*imag(aG2);b(3)=I2*alph2;
    b(4) = m3*real(aG3)-real(FP3); b(5)=m3*imag(aG3)-imag(FP3);
    b(6) = I3*alph3(i1)-real(RP3)*imag(FP3)+imag(RP3)*real(FP3)-T3;
    b(7) = m4*real(aG4)-real(FP4); b(8)=m4*imag(aG4)-imag(FP4);
    b(9) = I4*alph4(i1)-real(RP4)*imag(FP4)+imag(RP4)*real(FP4)-T4;
    X(:,i1) = A\b;
end


Ms = r1*X(8,:);
%Ms= r1*X(8,:)- I3*transpose(alph3) - I4*transpose(alph4); % shaking moment
Fs = sqrt((X(1,:)+X(7,:)).^2 + (X(2,:)+X(8,:)).^2); % magnitude of shaking force
Ys = -X(2,:)-X(8,:); Xs = -X(1,:)-X(7,:);
thetas = atan2(Ys,Xs);

%% Graphing

figure (1) % JOINT FORCE
subplot(2,2,1);plot(theta2, X(1,:), theta2, X(2,:));grid on; xlabel('crank angle'); subplot(2,2,2);plot(theta2, X(3,:),theta2, X(4,:) );grid on;xlabel('crank angle'); subplot(2,2,3);plot(theta2, X(5,:),theta2, X(6,:) );grid on;xlabel('crank angle'); subplot(2,2,4);plot(theta2, X(7,:),theta2, X(8,:) );grid on;xlabel('crank angle'); 

figure (2)
plot(theta2, -X(9,:)); grid on; xlabel('crank angle'); ylabel('T12');

figure(3) % shaking moment around O1, make external force and moment zero.
plot(theta2, Ms);grid on;
xlabel('crank angle');
ylabel('Shaking momemt')
hold on

figure(4) % shaking force, need make the external force and torque zero.
subplot (2,1,1);polarplot(thetas,Fs)
hold on
subplot (2,1,2);plot(theta2,Fs)
xlabel('crank angle')
ylabel('Shaking Force Magnitude')
