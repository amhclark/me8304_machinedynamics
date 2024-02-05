close all
clear

r2 = 1;
r3 = 1.41;
r4 = 0;

omega2 = 360 * 2 * pi / 60;
alph2 = 0;
T = 2 * pi / omega2;
n = 361;
t = linspace(0,T,n);
theta2 = omega2 * t;

theta3 = zeros(2,n); r1 = zeros(2,n);
omega3 = zeros(2,n); r1dot = zeros(2,n);
alph3 = zeros(2,n); r1ddot = zeros(2,n);


for i = 1:n

    theta3(1,i) = asin( (r2 * sin(theta2(i)) - r4) /r3 );
    theta3(2,i) = asin( -(r2 * sin(theta2(i))-r4) /r3 ) + pi;
    r1(:,i) = r2 * cos(theta2(i)) - r3 * cos(theta3(:,i));

    omega3(1:2,i) = r2 .*cos(theta2(i)) * omega2 / r3 ./cos(theta3(1:2,i));
    r1dot(1:2,i) = -r2 .*omega2 .*sin(theta2(i)) + r3 .*omega3(1:2,i) .*sin(theta3(1:2,i));

end

alph3(1,:) = (r2 * alph2 .*cos(theta2) - r2 * omega2^2 .*sin(theta2) + r3 .*omega3(1,:) .*omega3(1,:) .*sin(theta3(1,:))) ./r3 ./cos(theta3(1,:));
alph3(2,:) = (r2 * alph2 .*cos(theta2) - r2 * omega2^2 .*sin(theta2) + r3 .*omega3(2,:) .*omega3(2,:) .*sin(theta3(2,:))) ./r3 ./cos(theta3(2,:));

r1ddot(1,:) = -r2 * alph2 .*sin(theta2) - r2 * omega2^2 .*cos(theta2) +r3 .*alph3(1,:) .*sin(theta3(1,:)) + r3.*omega3(1,:) .*omega3(1,:) .*cos(theta3(1,:));
r1ddot(2,:) = -r2 * alph2 .*sin(theta2) - r2 * omega2^2 .*cos(theta2) +r3 .*alph3(2,:) .*sin(theta3(2,:)) + r3.*omega3(1,:) .*omega3(2,:) .*cos(theta3(2,:));
