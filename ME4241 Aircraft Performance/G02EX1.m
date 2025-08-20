% ME4241 Exercise 1
% Group 02

% Dai Baizhou Patrick A0266515B
% Tio Shao Wei A0255841B
% Justin Kyle Leonard A0266108E

% RUN THIS FILE TO OBTAIN ALL ANSWERS.
% Please run this single file.

% Answers can also be found at the last few rows.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[t, out] = ode45(@(t,h) sglide(t,h), [0 30*60], 1000);
vq = interp1(t,out,0:3000);
[CL,range] = Qa(1000);

[t, out] = ode45(@(t,x) unsglide(t,x), [0 30*60], [33; -pi/180; 0; 1000]);

[a, b] = ode45(@(t,x) unsglide_equa(t,x), [0 30*60], [33; -pi/180; 0; 1000]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Question (a)')
figure
plot(vq);
xlabel('Time(s)'); ylabel('Glider Height (m)'); title('Steady Flight');

disp('CL for maximum aerodynamic efficiency')
disp(CL);
disp('Glide range')
disp(range);
disp('Estimated time is 1520.8s (25 min 20.8s)')
fprintf('\n');

disp('Question (b)')

figure; 
subplot(4,1,1);
plot(t, out(:,1));
xlabel('Time (s)'); ylabel('Airspeed (m/s)'); title('Unsteady Flight for Standard Temperate Atmosphere');

subplot(4,1,2);
plot(t, out(:,2));
xlabel('Time (s)'); ylabel('Glide Angle (rad)');

subplot(4,1,3);
plot(t, out(:,3));
xlabel('Time (s)'); ylabel('Glide Distance (m)');

subplot(4,1,4);
plot(t, out(:,4));
xlabel('Time (s)'); ylabel('Glider Height (m)');

disp('Maximum variation in glide angle:')
disp('Approximately 0.056301 radian: from -0.0048864 radian to 0.007437 radian')
disp('Maximum variation in airspeed:')
disp('Approximately 1.3798 m/s: from 32.9989 m/s to 34.3787 m/s');
fprintf('\n');

disp('Question (c)')

disp('Glide range:');
disp('Approximately 50.149 km')
disp('Time:');
disp('Approximately 1524 seconds (25 minutes 24 seconds)')
fprintf('\n');

disp('Question (d)')
disp('Yes, performance is affected. Variations in glide angle and airspeed are more significant.')
disp('Maximum variation in glide angle:')
disp('0.114256 radian: from -0.079084 radian to 0.035172 radian')
disp('Maximum variation in airspeed:')
disp('2.8719 m/s: from 33 m/s to 35.8719 m/s')
disp('Glide range:');
disp('Apprxoimately 50.007 km')
disp('Time:');
disp('Approximately 1484 seconds (24 minutes 44 seconds)')

figure
subplot(4,1,1);
plot(a, b(:,1));
xlabel('Time (s)'); ylabel('Airspeed (m/s)'); title('Unsteady Flight for Standard Equatorial Atmospheric Conditions');

subplot(4,1,2);
plot(a, b(:,2));
xlabel('Time (s)'); ylabel('Glide Angle (rad)');

subplot(4,1,3);
plot(a, b(:,3));
xlabel('Time (s)'); ylabel('Glide Distance (m)');

subplot(4,1,4);
plot(a, b(:,4));
xlabel('Time (s)'); ylabel('Glider Height (m)');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function dhdt = sglide(t,h)

S = 11.81;
mass = 525.0;


CL = sqrt(0.00689/0.0145);
CD = 0.00689 + 0.0145*CL^2;


E = CL/CD;
angle = -atan(1/E);
range = h*E;

dhdt = sqrt((mass*9.81)/(0.5*ISA(h/1000)*S*(CL*cos(angle)+CD*sin(angle))))*sin(angle);


end



function nonsteadyglide = unsglide(t,x)

S = 11.81;
mass = 525;
CL = sqrt(0.00689/0.0145);
CD = 0.00689 + 0.0145 * CL^2;

nonsteadyglide = [((-1/2*ISA(x(4)/1000)*(x(1))^2*S*CD) - mass*9.81*sin(x(2)))/mass; 
                  ((1/2*ISA(x(4)/1000)*(x(1))^2*S*CL) - mass*9.81*cos(x(2)))/(mass*x(1));
                  x(1)*cos(x(2)); 
                  x(1)*sin(x(2))];

end



function nonsteadyglide = unsglide_equa(t,x)

S = 11.81;
mass = 525;
CL = sqrt(0.00689/0.0145);
CD = 0.00689 + 0.0145 * CL^2;

nonsteadyglide = [((-1/2*equa(x(4)/1000)*(x(1))^2*S*CD) - mass*9.81*sin(x(2)))/mass; 
                  ((1/2*equa(x(4)/1000)*(x(1))^2*S*CL) - mass*9.81*cos(x(2)))/(mass*x(1));
                  x(1)*cos(x(2)); 
                  x(1)*sin(x(2))];

end


function [den] = equa(h) % h in km
c = 6.51;
TSL = 303;
PSL = 101325;
R = 287;

T = TSL - c*h; %TSL: temp at stratosphere (constant)
P = PSL.*(T/TSL)^5.2506;
den = P/(R*T);

end


function density = ISA(h)

Psl = 101325;
Tsl = 288;
T = Tsl - 6.51 * h;
P = (T/Tsl)^5.2506 * Psl;

density = P/(T*287);

end


function [CL,range] = Qa(h)
CL = sqrt(0.00689/0.0145);
CD = 0.00689 + 0.0145*CL^2;


E = CL/CD;
angle = -atan(1/E);
range = h*E;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Exercise 1 Answers

% (a) CL for maximum aerodynamic efficiency is 0.6893. Glide range is
%     50.024km. Estimated time for gliding is 25 min and 20.8s (1520.8s).

% (b) The variation of the glide angle is approximately 0.056301 radian from -0.0048864
%     radian to 0.007437 radian
%     The variation of the airspeed is approximately 1.3798 m/s from 32.9989
%     m/s to 34.3787 m/s

% (c) The glide time for the unsteady glide is approximately 1524 seconds or
%     around 25 minutes 24 seconds and the glide range is approximately 50.149
%     km

% (d) Yes, if the glide glide occurs under standard equatorial equatorial
%     conditions, it will affect the glide performance. The airspeed variation
%     is 2.8719 m/s from 33 m/s to 35.8719 m/s. The glide angle angle variation
%     is 0.114256 radian from -0.079084 radian to 0.035172 radian. The glide
%     time decreases to approximately around 24 minutes 44 seconds and the glide
%     range decreases to apprxoimately 50.007 km
