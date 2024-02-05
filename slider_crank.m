%DYNAMIC FORCE ANALYSIS OF SLIDER CRANK MECHANISM
% based on the drawing on p539
clear;
load scdate
% CGs in local frame
rG2=r2/2; %bar 2
rG3=r3/2; %bar 3
miu=0.2; %friction coefficient
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% external load for bar 3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% external load for bar 4
FP4=0*exp(-45*pi/180*1j); %force on slider;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% masses and moment inertia of all bars
m2=2; m3=5;m4=10; %masses
I2=0.167; I3=0.817; %moment inertia of bar 2 and 3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% X=[F12x F12y F32x F32y F43x F43y F14y t12]
A=zeros(8,8);b=zeros(8,1);
X=zeros(8,n);
F12 = zeros(1,n);
F32=zeros(1,n);
F43=zeros(1,n);
F14=zeros(1,n);

A(1,1)=1; A(1,3)=1;
A(2,2)=1;A(2,4)=1;
A(3,8)=1;
A(4,3)=-1;A(4,5)=1;
A(5,4)=-1;A(5,6)=1;
A(7,5)=-1;
A(8,6)=-1;A(8,7)=1;
JJ=2; % the jjth configuration
for i1=1:n
    % acceleration of CG points
    
    aG2=rG2*alph2*(1j)*exp(theta2(i1)*1j)-rG2*omega2^2*exp(theta2(i1)*1j); % acceleration o
    aG4=r1ddot(JJ,i1); % acceleration of cG4
    aB=r2*alph2*(1j)*exp(theta2(i1)*1j)...
    -r2*omega2^2*exp(theta2(i1)*1j); % acceleration of point B
    aG3B=rG3*alph3(JJ,i1)*1j*exp((theta3(JJ,i1)+pi)*1j)...
    +rG3*omega3(JJ,i1)^2*exp(theta3(JJ,i1)*1j); % acceleration of point Cg3 relative to
    aG3=aB+aG3B;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5%%%
    %% vectors for each CG see FBD on p590
    R12=-rG2*exp(theta2(i1)*1j);
    R32=R12+r2*exp(theta2(i1)*1j);
    R43=-rG3*exp(theta3(JJ,i1)*1j);
    R23=R43+r3*exp(theta3(JJ,i1)*1j);
    %%%%%
    A(3,1)=-imag(R12);A(3,2)=real(R12);A(3,3)=-imag(R32);A(3,4)=real(R32);
    A(6,3)=imag(R23);A(6,4)=-real(R23);A(6,5)=-imag(R43);A(6,6)=real(R43);
    A(7,7)= -sign(r1dot(JJ,i1))*miu;
    b(1)=m2*real(aG2); b(2)=m2*imag(aG2);b(3)=I2*alph2;
    b(4)=m3*real(aG3); b(5)=m3*imag(aG3);
    b(6)=I3*alph3(JJ,i1);
    b(7)=m4*real(aG4)-real(FP4); b(8)=-imag(FP4);
    X(:,i1)=A\b;
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% computer and draw pin forces
F12 =sqrt(X(1,:).^2+X(2,:).^2);
F32=sqrt(X(3,:).^2+X(4,:).^2);
F43=sqrt(X(5,:).^2+X(6,:).^2);
F14=sqrt(X(7,:).^2+(miu*X(7,:)).^2);
%%% figures of pin forces
figure (1)
subplot(2,2,1)
plot(theta2,F12);grid;
xlim([0 2*pi]);
xlabel('theta 2 [rad]')
ylabel('F12 [N]')
subplot(2,2,2)
plot(theta2,F32);grid;
xlim([0 2*pi]);
xlabel('theta 2 [rad]')
ylabel('F32 [N]')
subplot(2,2,3);
plot(theta2,F43);grid;
xlim([0 2*pi]);
xlabel('theta 2 [rad]')
ylabel('F43 [N]')
%%
subplot(2,2,4)
plot(theta2,F14);grid;
xlim([0 2*pi]);
xlabel('theta 2 [rad]')
ylabel('F14 [N]')
% figure of input torque
figure (2)
plot(theta2, X(8,:)); grid;
xlabel('theta 2 [rad]')
ylabel('T12 [Nm]')
%END OF CODE
