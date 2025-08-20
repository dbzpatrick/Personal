% Exercise 2
% Group 02
% Dai Baizhou Patrick A0266515B
% Justin Leonard A0266108E
% Tio Shao Wei A0255841B

% Run this single file to obtain figures
% Answers to exercise in last few rows

clear all

S = 17.1;
mass = 1248.5;
CLmax = 3;
W = mass * 9.81;
pset = 0.5; %0.3 0.7

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

guess1 = [0.1, 50];
sol1 = fsolve(@(var) findclsea(var, 0), guess1);
CLguess1 = sol1(1);
Vguess1 = sol1(2); 

disp(['CL1 at sea level = ', num2str(CLguess1)]);
disp(['V1 at sea level = ', num2str(Vguess1)]);

guess2 = [2, 10];
sol2 = fsolve(@(var) findclsea(var, 0), guess2);
CLguess2 = sol2(1);
Vguess2 = sol2(2); 

disp(['CL2 at sea level = ', num2str(CLguess2)]);
disp(['V2 at sea level = ', num2str(Vguess2)]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CL = [];
v = [];
h = [];
guess = [50, 30]; %initial guess
for i = CLguess1:0.05:min(CLmax,CLguess2)
    CL = [CL;i];

    solution = fsolve(@(var) findenvelope(var, i), guess);
    guess = [solution(1), solution(2)];
    v = [v;guess(1)];
    h = [h;guess(2)];

end
plot(v, h);
xlim([0 100]);
ylim([0 7000]);

hold on

%Add labels
xlabel('Velocity (m/s)'); title('Question 2');
ylabel('Height (m)');


redline = zeros(300,2);%30
for i = 1:300 %i = 1:30 
    vel = sqrt((mass*9.81)/(0.5*P2((i*10)/1000)*S*CLmax));%100
    redline(i,2) = i*10;
    redline(i,1) = vel;
end
hold on;
plot(redline(:,1),redline(:,2));

spline_interp = spline(v,h);                       % connect the points
f_to_minimize = @(x) -ppval(spline_interp, x);     % function of negative spline (guess @ variable x)
max_x = fminsearch(f_to_minimize, 40);             % find min of function, aka max of real spline
max_y = -f_to_minimize(max_x);                     % find y at x
fprintf("ABS Ceiling is "+max_y/1000 +"km \n")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function F2 =  findenvelope(var, CL)

v = var(1);
h = var(2);

S = 17.1;
mass = 1248.5;
pset = 0.5; %0.3 0.7
W = mass * 9.81;
Tsl = 288;
T = Tsl - 6.51 * (h/1000);
%Speed of sound
a = sqrt(1.4 * 287 * T);

%Calculating the Thrust on sea level
c0 = 7106.47 - 0.781351 * h;
c1 = -2094.99 + 0.216186 * h;
c2 = -3638.48 + 0.45727 * h;
T = pset * (c0 + (v/a) * (c1 + c2 * (v/a)));

%Equations to solve
F2(1) = 0.5*CL*P2(h/1000)*v^2*S - W;
F2(2) = 0.5*(0.036 + 0.06 * CL^2)*P2(h/1000)*v^2*S - T;

end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function F =  findclsea(var, h)

clg = var(1);
vg = var(2);

S = 17.1;
mass = 1248.5;
pset = 0.5; %0.3 0.7
W = mass * 9.81;
Tsl = 288;
T = Tsl - 6.51 * h;
%Speed of sound
a = sqrt(1.4 * 287 * T);

%Calculating the Thrust on sea level
c0 = 7106.47 - 0.781351 * h;
c1 = -2094.99 + 0.216186 * h;
c2 = -3638.48 + 0.45727 * h;
T = pset * (c0 + (vg/a) * (c1 + c2 * (vg/a)));

%Equations to solve
F(1) = 0.5*clg*P2(h/1000)*vg^2*S - W;
F(2) = 0.5*(0.036 + 0.06 * clg^2)*P2(h/1000)*vg^2*S - T;

end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function density = P2(h)

Tsl = 288;
Psl = 101325;
T = Tsl - 6.51 * h;

P = (T/Tsl)^5.2506 * Psl;

density = P/(T*287);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Q1
% CL Trim 1 at sea level = 0.14488
% The velocity trim 1 at sea level = 89.8101 m/s
% CL trim 2 at sea level = 4.6338
% The velocity trim 2 at sea level = 15.8802 m/s

% Q2
% The flight envelope can be observed after running the code

% Q3
%The absolute ceiling of the UAV is 5.9794 km or 5979.4 m

% Q4
% When the power setting is increased to 0.7
% CL Trim 1 at sea level = 0.10504
% The velocity trim 1 at sea level = 105.4761 m/s
% CL trim 2 at sea level = 6.6949
% The velocity trim 2 at sea level = 13.3113 m/s
% The absolute ceiling is increased to 6.8559 km or 6855.9 m
% From the flight envelope graph, it can also be seen that the higher power setting will have a bigger boundary for the flight envelope, which means the flight envelope enlarges

% When the power setting is increased to 0.3
% CL Trim 1 at sea level = 0.24792
% The velocity trim 1 at sea level = 68.6539 m/s
% CL trim 2 at sea level = 2.6126
% The velocity trim 2 at sea level = 21.14191 m/s
% The absolute ceiling is increased to 3.9596 km or 3959.6 m
% From the flight envelope graph, it can also be seen that the lower power setting will have a smaller boundary for the flight envelope, which means the flight envelope becomes smaller
