
%Group 02
%Exercise 3
% Tio Shao Wei A0255841B
% Dai Baizhou Patrick A0266515B
% Justin Kyle Leonard A0266108E
clear all


CLm = 3;

%Guessing Clb
%Change the altitude according to the set altitude
guess1 = [0.1, 50];
sol1 = fsolve(@(var) levelflight(var, 0), guess1);
Clb = sol1(1);
v1 = sol1(2);

%Guessing Cla
%Change the altitude according to the set altitude
guess2 = [2, 10];
sol2 = fsolve(@(var) levelflight(var, 0), guess2);
Cla = sol2(1);
v2 = sol2(2);

disp(['CLb = ', num2str(Clb)]);
disp(['V1 = ', num2str(v1)]);

disp(['CLa = ', num2str(Cla)]);
disp(['V2 = ', num2str(v2)]);

RC = [];
lv = [];
lgamma =[];

guess3 = [pi/180, 30];
%Generating the climb envelope
for CL = Clb:0.01:min(Cla,CLm)
    sol3 = fsolve(@(var) climbflight(var, CL), guess3);
    guess3 = [sol3(1), sol3(2)];
    disp(guess3);
    RC = [RC; sol3(2) * sin(sol3(1))];
    lv = [lv; sol3(2)];
    lgamma = [lgamma; sol3(1)];

end

%Calculating the maximum climb rate
spline_interp = spline(lv, RC);
f_to_minimize = @(x) -ppval(spline_interp, x);
x_max = fminsearch(f_to_minimize, 50);
y_max = -f_to_minimize(x_max);

%Calculating steepest climb angle
g1_max = max(lgamma); %Radian 
g1_max_deg = g1_max*(180/pi);

figure(1)
plot(lv,RC);
xlim([0, 100]);
ylim([0, 10]);
hold on
plot([0,40],[0,40*tan(g1_max)]); % plot tangent line

%Add labels
xlabel('Velocity (m/s)');
ylabel('Rate of Climb');
hold off

%Generating Height vs ROC
lh = [];
RCm = [];
lg = [];

guess4 = [pi/180, 30];
for h = 0:100:5979
    lv1 = []; %New list for velocity at particular height
    RC1 = []; %New list for RC at particular height
    lg1 = []; %New list for gamma at particular height

    for CL = Clb:0.1:min(Cla, CLm)
        sol4 = fsolve(@(var) climbflight2(var, CL, h), guess4);
        guess4 = [sol4(1), sol4(2)];
        lv1 = [lv1, sol4(2)];
        RC1 = [RC1; sol4(2) * sol4(1)];
        lg1 = [lg1; sol4(1)];
       
    end

    %Calculating the maximum climb rate
    spline_interp1 = spline(lv1, RC1);
    f_to_minimize = @(x) -ppval(spline_interp1, x);
    x_max = fminsearch(f_to_minimize, 50);
    RC_max = -f_to_minimize(x_max);
    lh = [lh; h];
    RCm = [RCm; RC_max];

    %Calculating steepest climb angle
    g_max = max(lg1); %Radian 
    g_max_deg = g_max*(180/pi); %Degree
    lg = [lg; g_max_deg];

end

figure(2)
plot(RCm, lh);
hold on;
xq = [0.7; 0.6 ; 0.5; 0.4; 0.3; 0.2; 0.1]; 
yq = interp1(RCm, lh, xq, 'spline');
plot(xq,yq);
y_service = interp1(RCm, lh, 0.503); %Height of service ceiling
xline(0.503);
hold off;
text(0.5,5661.3,'\leftarrow Service Ceiling');
%Add labels
xlabel('Rate of Climb (m/s)');
ylabel('Altitude (h)');

%Calculating the flight time
RC1 = [RCm ; xq];
HSC = [lh ; yq]; 
pp = spline(RC1, HSC);
xq = 0:0.5:9;
figure
plot(xq,ppval(pp,xq));

%RC at service ceiling is 100 ft/min (30.5m/min), (0.5083m/s)
RC_ceil = 0.5083;
yatintersect = ppval(pp,RC_ceil);
disp(yatintersect);


%yatintersect= 5661;

pp1 = spline(HSC,RC1);
%plot(spline1);
%plot(spline, linspace())
yinterpolated = integral(@(z) (1./ppval(pp1,z)), 0, yatintersect);
ytime = yinterpolated/60;


disp(['Max ROC = ', num2str(y_max), ' m/s'])
disp(['Steepest Climb Angle = ', num2str(g1_max_deg), ' deg'])
disp(['The Surface Ceiling = ' , num2str(y_service), ' m'])
disp(['Estimated Time to Reach Service Ceiling = ', num2str(ytime), ' mins'])
 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function F = levelflight(var, h)

clg = var(1);
vg = var(2);

S = 17.1;
mass = 1248.5;
pset = 0.5;
Tsl = 288;
T = Tsl - 6.51 * (h/1000);
W = 9.81 * mass;


%Calculating Speed of Sound
a = sqrt(1.4 * 287 * T);

%Calculating Thrust
c0 = 7106.47 - 0.781351 * h;
c1 = -2094.99 + 0.216186 * h;
c2 = -3638.48 + 0.45727 * h;

T = pset * (c0 + (vg/a) * (c1 + c2 * (vg/a)));

%Equations to solve
F(1) = 0.5*clg*P2(h/1000)*vg^2*S - W;
F(2) = 0.5*(0.036 + 0.06 * clg^2)*P2(h/1000)*vg^2*S - T;

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function F2 = climbflight(var, CL)
gamma = var(1);
velocity = var(2);

h = 0; %change according to the set altitude
S = 17.1;
mass = 1248.5;
pset = 0.5;
Tsl = 288;
T = Tsl - 6.51 * (h/1000);
W = 9.81 * mass;

%Calculating Speed of Sound
a = sqrt(1.4 * 287 * T);

%Calculating Thrust
c0 = 7106.47 - 0.781351 * h;
c1 = -2094.99 + 0.216186 * h;
c2 = -3638.48 + 0.45727 * h;

T = pset * (c0 + (velocity/a) * (c1 + c2 * (velocity/a)));

%Equations to solve
F2(1) = T - 0.5*(0.036 + 0.06 * CL^2)*P2(h/1000)*(velocity^2)*S - W * sin(gamma);
F2(2) = 0.5*CL*P2(h/1000)*(velocity^2)*S - W * cos(gamma);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function F3 = climbflight2(var, CL, h)
gamma = var(1);
velocity = var(2);

S = 17.1;
mass = 1248.5;
pset = 0.5;
Tsl = 288;
T = Tsl - 6.51 * (h/1000);
W = 9.81 * mass;

%Calculating Speed of Sound
a = sqrt(1.4 * 287 * T);

%Calculating Thrust
c0 = 7106.47 - 0.781351 * h;
c1 = -2094.99 + 0.216186 * h;
c2 = -3638.48 + 0.45727 * h;

T = pset * (c0 + (velocity/a) * (c1 + c2 * (velocity/a)));

%Equations to solve
F3(1) = T - 0.5*(0.036 + 0.06 * CL^2)*P2(h/1000)*(velocity^2)*S - W * sin(gamma);
F3(2) = 0.5*CL*P2(h/1000)*(velocity^2)*S - W * cos(gamma);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function density = P2(h)

Tsl = 288;
Psl = 101325;
T = Tsl - 6.51 * h;

P = (T/Tsl)^5.2506 * Psl;

density = P/(T*287);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%The Maximum Rate of Climb = 8.6153 m/s
%The Steepest Climb Angle = 10.8049 deg
%The Surface Ceiling = 5659.3472 m
%Estimated Time to Reach Service Ceiling = 31.8352 mins