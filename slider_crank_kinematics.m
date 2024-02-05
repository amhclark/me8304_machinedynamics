% Adapted from Dr. J. Yang - Memorial University
% Machine Dynamics

close all
clear

%% Input Data

r2 = 1;
r3 = 1.41;
r4 = 0;                                     % offset if it is non-offset r4 = 0
                                            % R1 is unknown

omega2 = 360 * 2 * pi / 60;
alph2 = 0;                                  % angular acceleration
T = 2 * pi / omega2;                        % period of the crank
n = 361;                                    % number of time values
t = linspace(0,T,n);                        % time vector 0 to 2 seconds, 400 points
theta2 = omega2 * t;                        % array of crank angles in rad

theta3 = zeros(2,n); r1 = zeros(2,n);       % th(i,:), r1(i,:), the ith solution
omega3 = zeros(2,n); r1dot = zeros(2,n);    % initialize the velocity output
alph3 = zeros(2,n); r1ddot = zeros(2,n);    % initialize the acceleration output


%% Position Analysis

for i = 1:n

    % Position Analysis (p161)
    theta3(1,i) = asin( (r2 * sin(theta2(i)) - r4) /r3 );
    theta3(2,i) = asin( -(r2 * sin(theta2(i))-r4) /r3 ) + pi;
    r1(:,i) = r2 * cos(theta2(i)) - r3 * cos(theta3(:,i));

    % Velocity Analysis (p275)
    omega3(1:2,i) = r2 .*cos(theta2(i)) * omega2 / r3 ./cos(theta3(1:2,i));
    r1dot(1:2,i) = -r2 .*omega2 .*sin(theta2(i)) + r3 .*omega3(1:2,i) .*sin(theta3(1:2,i));

end


%% Acceleration Analysis

alph3(1,:) = (r2 * alph2 .*cos(theta2) - r2 * omega2^2 .*sin(theta2) + r3 .*omega3(1,:) .*omega3(1,:) .*sin(theta3(1,:))) ./r3 ./cos(theta3(1,:));
alph3(2,:) = (r2 * alph2 .*cos(theta2) - r2 * omega2^2 .*sin(theta2) + r3 .*omega3(2,:) .*omega3(2,:) .*sin(theta3(2,:))) ./r3 ./cos(theta3(2,:));

r1ddot(1,:) = -r2 * alph2 .*sin(theta2) - r2 * omega2^2 .*cos(theta2) +r3 .*alph3(1,:) .*sin(theta3(1,:)) + r3.*omega3(1,:) .*omega3(1,:) .*cos(theta3(1,:));
r1ddot(2,:) = -r2 * alph2 .*sin(theta2) - r2 * omega2^2 .*cos(theta2) +r3 .*alph3(2,:) .*sin(theta3(2,:)) + r3.*omega3(1,:) .*omega3(2,:) .*cos(theta3(2,:));


%% Position Plot
% Plot positions vs theta

figure(1)
subplot(2,2,1)
plot(theta2,r1(1,:));
xlabel('theta 2 [rad]')
ylabel('r1 [1,:]')

subplot(2,2,2)
plot(theta2,theta3(1,:));
xlabel('theta 2 [rad]')
ylabel('theta 3 (1,:) [rad]')

subplot(2,2,3)
plot(theta2,r1(2,:));
xlabel('theta 2 [rad]')
ylabel('theta 3(2,:) [rad]')

subplot(2,2,4)
plot(theta2,theta3(2,:));
xlabel('theta 2 [rad]')
ylabel('theta 3(2,:) [rad]')


%% Velocity Plot
% Plot veloccities vs theta

figure(2)
subplot(2,2,1)
plot(theta2,omega3(1,:));
xlabel('theta 2 [rad]')
ylabel('r1dot (1,:) [m/s]')

subplot(2,2,2)
plot(theta2,r1dot(1,:));
xlabel('theta 2 [rad]')
ylabel('r1dot (1,:) [m/s]')

subplot(2,2,3)
plot(theta2,omega3(2,:));
xlabel('theta 2 [rad]')
ylabel('omega3(2,:) [rad/s]')

subplot(2,2,4)
plot(theta2,r1dot(2,:));
xlabel('theta2 [rad]')
ylabel('r1dot(2,:) [m/s]')


%% Acceleration Plot
% Plot accelerations vs theta

figure(3)
subplot(2,2,1)
plot(theta2,alph3(1,:));
xlabel('theta 2 [rad]')
ylabel('alph3(1,:) [rad/s^2')

subplot(2,2,2)
plot(theta2,r1ddot(1,:));
xlabel('theta 2 [rad]')
ylabel('r1ddot (1,:) [m/s^2]')

subplot(2,2,3)
plot(theta2,alph3(2,:));
xlabel('theta 2 [rad]')
ylabel('alph3(2,:) [rad/s^2]')

subplot(2,2,4)
plot(theta2,r1ddot(2,:));
xlabel('theta 2 [rad]')
ylabel('r1ddot(2,:) [m/s^2]')

%% Save Data
save scdate n r1 r1dot r1ddot r2 r3 r4 omega2 omega3 alph2 alph3 theta2 theta3 t;
