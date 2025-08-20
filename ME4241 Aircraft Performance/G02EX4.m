% Group 02
% Exercise 4
% Tio Shao Wei A0255841B
% Dai Baizhou Patrick A0266515B
% Justin Kyle Leonard A0266108E

% Please run this file to see all answers.
% Full answers attached at the end of this file.

clear all;

S = 17.1;
m = 1248.5;
Tsl = 288;
Psl = 101325;
clm = 3; %Maximum CL
Nm = 3; %Maximum Load factor

options = optimoptions('fsolve','Display','off'); % turn off fsolve status

%Guessing clb
%Change the altitude according to the set altitude
guess1 = [0.1, 50];
sol1 = fsolve(@(var) levelflight(var, 0), guess1, options);
clb = sol1(1);
vb = sol1(2);

%Guessing Cla
%Change the altitude according to the set altitude
guess2 = [2, 10];
sol2 = fsolve(@(var) levelflight(var, 0), guess2, options);
cla = sol2(1);
va = sol2(2);

%Generating turn envelope
pt = 285; %Number of iteration points
cln = linspace(clb, min(cla, clm), pt);
A = zeros(pt,5); %Generating nx5 array

guess3 = [pi/180, vb];
for i = 1:pt
  sol3 = fsolve(@(var) turningflight(var, cln(i)), guess3, options);
  [rho, temp] = P2(0/1000);
  w = (0.5*rho*sol3(2)^2*S*cln(i)*sin(sol3(1)))/ (m * sol3(2));
  L = 0.5*rho*sol3(2)^2*S*cln(i);
  R = sol3(2)/w;
  A(i, 1) = sol3(1);
  A(i, 2) = sol3(2);
  A(i, 3) = w * (180/pi); %Change the turn rate to deg/s
  A(i, 4) = L;
  A(i, 5) = R;
  guess3 = [sol3(1), sol3(2)]; 
end

figure(1)
plot(A(:,2) , A(:,3),'DisplayName', sprintf('Turn Envelope'));
title('Sea Level Turn Envelope');
xlabel('Velocity (m/s)');
ylabel('Omega (deg/s)');
hold on

%Maximum turn rate
pp = spline(A(:,2) , A(:,3));
func = @(x) -ppval(pp, x);
xmax = fminsearch(func, 50); 
tmax = ppval(pp, xmax);

%Smallest turn radius
rmin = min(A(:,5));
index = find(A(:,5) == rmin); 
vrmin = A(index, 2);
wrmin = A(index, 3);

yline(tmax , 'DisplayName', 'Max turn rate');
plot([0,40],[0, 40/(rmin) * 180/pi], 'DisplayName', 'Tangent Line');

%Largest load factor
lmax = max(A(:,4));
nmax = lmax/(m*9.81);
index2 = find(A(:,4) == lmax); 
vlmax = A(index2, 2); 
wlmax = A(index2, 3); 

%Plotting the load factor line
nline = zeros(pt,2);
Vn = linspace(10, 100, pt); %List of velocities
n = 3;
L = 3*9.81*m;
for n = 1.5:0.5:Nm
    for i = 1:pt
        wn = (9.81/Vn(i)) * sqrt(n^2 -1);
        nline(i, 1) = wn * (180/pi); %Change to deg/s for turn rate
        nline(i, 2) = Vn(i);
    end
    
    plot(nline(:, 2), nline(:, 1),'.-', 'DisplayName', sprintf('Load factor =  %.1f', n));

end

pt2 = 500; %Number of iteration points
Vn2 = linspace(10, 100, pt2);
Lm = zeros(pt2, 1);
[rho, temp] = P2(0/1000);
nm = zeros(pt2, 1);
for j = 1:pt2
    Lm(j,1) = 0.5*rho*Vn2(j)^2*S*clm;
    nm(j,1) = Lm(j,1)/(m*9.81);
end

n3line = zeros(pt2,2);
for k = 1:pt2
    if (nm(k,1))^2 < 1
        continue
    else
        w3 = (9.81/Vn2(k)) * sqrt((nm(k,1))^2 -1);
    n3line(k, 1) = w3 * (180/pi); %Change to deg/s for turn rate
    n3line(k, 2) = Vn2(k);
    end

end

plot(n3line(:,2) , n3line(:,1), '--r', 'DisplayName', sprintf('CL = %.1f', clm));
xlim([0, 100]);
ylim([0, 60]);

legend show;
hold off;

%Find the corner speed and turn rate from the intersection
[vcorner, tcorner] = polyxpoly(n3line(:,2), n3line(:,1), nline(:,2), nline(:,1));
text(vcorner, tcorner, '\leftarrow Turn Rate at Corner Speed');

%Calculate specific excess power under the intersection
pset = 0.5;

%Calculating speed of sound
a = sqrt(1.4 * 287 * temp);

%Calculating Thrust
c0 = 7106.47 - 0.781351 * 0; 
c1 = -2094.00 + 0.216186 * 0;
c2 = -3638.48 + 0.45727 * 0;
T1 = pset * (c0 + vcorner/a * (c1 + c2 * (vcorner/a)));
D1 = 0.5*rho*(vcorner)^2*S*(0.036 + 0.06 * clm^2);
Px = vcorner*(T1 - D1)/(m*9.81);


fprintf('The max turn rate is = %f deg/s at V = %f m/s \n',tmax, xmax);
fprintf('The smallest turn radius is = %f m at V = %f m/s, w = %f deg/s \n',rmin, vrmin, wrmin);
fprintf('The max load factor is = %f at V = %f m/s, w = %f deg/s \n',nmax, vlmax, wlmax);
fprintf('The corner speed is = %f m/s and the corner turn rate is %f deg/s \n', vcorner, tcorner);    
fprintf('The specific excess power is = %f m/s \n', Px);

%Plot line of different turn radius
for R = 40:50:300
    wlist = [];
    Vlist = [];
    for V = 10:100
        w = (180/pi)*V/R;
        wlist = [wlist,w];
        Vlist = [Vlist,V];
    end
    hold on
    plot(Vlist,wlist,":",'DisplayName',sprintf('Turn Radius = %.1f m',R));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function F1 = levelflight(var, h)

cl1 = var(1);
v1 = var(2); 

S = 17.1;
m = 1248.5;
pset = 0.5;
[rho, temp] = P2(h/1000);
  

%Calculating speed of sound
a = sqrt(1.4 * 287 * temp);

%Calculating Thrust
c0 = 7106.47 - 0.781351 * h;
c1 = -2094.00 + 0.216186 * h;
c2 = -3638.48 + 0.45727 * h;
T = pset * (c0 + v1/a * (c1 + c2 * (v1/a)));

%Equations to solve
F1(1) = T - 0.5*rho*v1^2*S*(0.036 + 0.06 * cl1^2);
F1(2) = 0.5*rho*v1^2*S*cl1 - m*9.81;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function F2 = turningflight(var, cl)

miu1 = var(1); %Bank angle
v2 = var(2);   %Velocity

h = 0;
S = 17.1;
m = 1248.5;
pset = 0.5;
[rho, temp] = P2(h/1000);
  

%Calculating speed of sound
a = sqrt(1.4 * 287 * temp);

%Calculating Thrust
c0 = 7106.47 - 0.781351 * h;
c1 = -2094.00 + 0.216186 * h;
c2 = -3638.48 + 0.45727 * h;
T = pset * (c0 + v2/a * (c1 + c2 * (v2/a)));

%Equations to solve
F2(1) = T - 0.5*rho*v2^2*S*(0.036 + 0.06 * cl^2);
F2(2) = 0.5*rho*v2^2*S*cl*cos(miu1) - m*9.81;


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function F3 = loadfactor(var, L, miu3)

w3 = var(1); %Turn rate
v3 = var(2); %Velocity

h = 0;
S = 17.1;
m = 1248.5;
pset = 0.5;
[rho, temp] = P2(h/1000);

%Calculating speed of sound
a = sqrt(1.4 * 287 * temp);

%Calculating Thrust
c0 = 7106.47 - 0.781351 * h;
c1 = -2094.00 + 0.216186 * h;
c2 = -3638.48 + 0.45727 * h;
T = pset * (c0 + v3/a * (c1 + c2 * (v3/a)));

%Equations to solve
F3(1) = L*sin(miu3) - m*v3*w3;
F3(2) = L*cos(miu3) - m*9.81;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [density, T] = P2(h)

Tsl = 288;
Psl = 101325;
T = Tsl - 6.51 * h;

P = (T/Tsl)^5.2506 * Psl;

density = P/(T*287);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Question a 
% Envelope can be obtained by running this file.

% Question b
% The max turn rate is = 28.926805 deg/s at 37.503586 m/s 
% The smallest turn radius is = 53.991050 m at V = 23.975882 m/s, w = 25.443418 deg/s 
% The max load factor is = 2.886277 at V = 64.194840 m/s, w = 23.706153 deg/s 

% Question c
% It is within the boundaries
% The corner speed is = 34.184431 m/s and the corner turn rate is 46.506784 deg/s 
% The specific excess power is = -10.118178 m/s 

% Question d
% When UAV operates at corner speed, it will still able to turn but its altitude will decrease.
% This is because its specific excess power is negative at this condition.