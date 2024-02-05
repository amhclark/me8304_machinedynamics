% This code is to analyze kinematics of four bar linkage:
% Method: Newton-Raphson
global r1 r2 r3 r4
global theta2inst
% Constants
r1 = 5.5;r2 =2.0; r3 = 6.0; r4 = 3.0; % length of the four bars
n=361; % number of increments
%knowns for links
omega2=50; % rad/s
alph2=0; % angular acceleration
T2= 2*pi/omega2; % period of bar 2
t = (linspace(0,T2,n))'; % discritize time
theta2=omega2*t;
%initialize omega3 and 4
omega3=zeros(n,1);omega4=zeros(n,1); % anular speed of bar 3 and 4
alph3=zeros(n,1);alph4=zeros(n,1); % anular acceleration of bar 3 and 4
% theta30=acos((-r4^2+r3^2+(r1-r2)^2)/(2*r3*(r1-r2))); % guess degree
theta30= 23*pi/180;
theta40= 108*pi/180;
%theta40=pi-acos(((r1-r2)^2+r4^2-r3^2)/(2*r4*(r1-r2))) + 10*pi/180; % guess degree
theta3=zeros(n,1);theta4=zeros(n,1);
theta3(1)=theta30;
theta4(1)=theta40;

% position analysis using Newton-Raphson
for i1=1:n
    theta2inst=theta2(i1);
    % Move unknown coordinates to array x
    if i1==1;
        x=[theta30; theta40];
        else
        x = [theta3(i1-1); theta4(i1-1)];
        end

    for n1 = 1:80
        % Compute 2 x 1 array of functions
        f = constraints(x);
        % condition for termination:
        normf = norm(f);
        if ( normf <= 1e-7 ) break; end;
            % Construct 2 x 2 Jacobian
            D = jacobian(x);
            % Compute corrections
            delta_x = D\f;
            % Correct the coordinates
            x = x - delta_x;
            end

    % Report solution in degrees
    theta3(i1) = x(1); % in rad
    theta4(i1) = x(2); % in rad
end

% velocity analysis: P273
for i1=1:n
    A1=[-r3*sin(theta3(i1)), r4*sin(theta4(i1));r3*cos(theta3(i1)), -r4*cos(theta4(i1))]
    b1=[r2*omega2*sin(theta2(i1)); -r2*omega2*cos(theta2(i1))];
    xx=A1\b1;
    omega3(i1)=xx(1);omega4(i1)=xx(2);
end

% acceleration analysis :P310-311
for i1=1:n
    A1=[-r3*sin(theta3(i1)), r4*sin(theta4(i1));r3*cos(theta3(i1)), -r4*cos(theta4(i1))]
    b1(1,1)=r2*alph2*sin(theta2(i1))+r2*omega2^2*cos(theta2(i1))+r3*omega3(i1)^2*cos(theta3)
    -r4*omega4(i1)^2*cos(theta4(i1));
    b1(2,1)=-r2*alph2*cos(theta2(i1))+r2*omega2^2*sin(theta2(i1))...
    + r3*omega3(i1)^2*sin(theta3(i1))-r4*omega4(i1)^2*sin(theta4(i1));
    xx=A1\b1;
    alph3(i1)=xx(1);alph4(i1)=xx(2);
end

% Plot for position
figure(1)
thetaplot=plot(theta2,theta3,theta2,theta4 );
xlabel('theta 2 [rad]');
ylabel('theta 3 (theta4) [rad]');
hkids=get(gca,'child');
legend('theta3','theta4');
grid on;

% Plot for velocity
figure(2)
plot(theta2,omega3/omega2,theta2,omega4/omega2 );
xlabel('theta2 [rad]');
ylabel('angular speed ratio');
legend('omega3/omega2','omega4/omega2');
grid on;
% Plot for acceleration
figure(3)
plot(theta2,alph3,theta2,alph4);
xlabel('theta2 [rad]');
ylabel('alph3 (4) [rad/s^2]');
legend('alph3','alph4');
grid on;
save fourbardate n t r1 r2 r3 r4 omega2 omega3 omega4 alph2 alph3 alph4 theta2 theta3 theta4

function f = constraints(x)
% Evaluate the constraints
global r1 r2 r3 r4
global theta2inst
theta3 = x(1);
theta4 = x(2);
f1 = r2*cos(theta2inst) + r3*cos(theta3) - r4*cos(theta4) - r1;
f2 = r2*sin(theta2inst) + r3*sin(theta3) - r4*sin(theta4);
f = [f1 f2]';
end

function D = jacobian(x)
% Evaluate the Jacobian
global r1 r2 r3 r4
theta3 = x(1);
theta4 = x(2);
d11 = -r3*sin(theta3);
d12 = r4*sin(theta4);
d21 = r3*cos(theta3);
d22 = -r4*cos(theta4);
D = [d11 d12; d21 d22];
end
