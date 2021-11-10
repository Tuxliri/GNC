% Spacecraft Guidance and Navigation (2021/2022)
% Assignment # 1
% Author: Davide Iafrate

%% Ex 1 - OKAY
clearvars; close all; clc

% Load kernels
cspice_furnsh('assignment01.tm');

%%%%%%%%% #1 PROPAGATOR VALIDATION %%%%%%%%%%
primary.label = 'Earth';
primary.GM = cspice_bodvrd(primary.label,'GM',1);   % Gravitational param
                                                    %   [km^3/s^2]
r0 = [ -7128.137, 0, 0 ]';   % Initial radius vector             [ km ]
v0 = [ 0, -9.781, 0 ]';      % Initial velocity vector           [ km/s ]

% Specific energy (constant)                                    [ km^2/s^2 ]
epsilon0 = 0.5*norm(v0)^2 - primary.GM/norm(r0);

a = - primary.GM/(2*epsilon0);        % Semi-major axis                   [ km ]
T = 2*pi*sqrt(a^3/primary.GM);       % Period                            [ s ]

% Initial conditions
x0 = [r0;v0];
t0 = 0;

xend = flow2BP(t0,x0,T,primary.GM);

% If error = 0 after one period our propagator is working correctly
error = norm(xend-x0);

%% 2 -  LAMBERT PROBLEM SHOOTING SOLVER %%%%%%%%%
%(Earth centered solution) - Validated
% Primary parameters 
primary.label = 'Earth';
primary.GM = cspice_bodvrd(primary.label,'GM',1);   % Gravitational param

% Retrieve initial and final positions of the planets
% Initial position vector [km]
r1 = [ 7000; 0; 0];
v1_0 = [0; 8.5; 0];
t1 = 0;
ToF = 2000;

% Final position vector [km]
r2 = [200; 7000; -100];

% Initial guess as Hohmann transfer velocity
Earth.objective = @(x) objectiveFun(x,0,ToF,primary.GM,r2,r1);
Earth.x0 = v1_0;
Earth.solver = 'fsolve';

% Set options for using finite difference method
optsFinDiff = optimoptions('fsolve');
optsFinDiff.Display='iter';
optsFinDiff.FiniteDifferenceType = 'central';

Earth.options = optsFinDiff;

solFD=fsolve(Earth);

% Set options for using the Jacobian computed through the STM
optsJacobian = optimoptions('fsolve');
optsJacobian.Display='iter';
optsJacobian.SpecifyObjectiveGradient = true;

Earth.options = optsJacobian;

solSTM = fsolve(Earth);

%%%%%%%%%%%% Earth-to-Mars Lambert validator - Validated %%%%%%%%%%%
primary.label = 'Sun';
primary.GM = cspice_bodvrd(primary.label,'GM',1);   % Gravitational param
            
% select frame string (SPICE naming convention)
frame = 'ECLIPJ2000';
% select departure and arrival strings
departure.Label = 'Earth';
arrival.Label = 'Mars';

departure.EpochStr = '1996-Jul-7 00:00:00.0000 TDB';
arrival.EpochStr = '1996-Dec-12 00:00:00.0000 TDB';

% Times
departure.time = cspice_str2et(departure.EpochStr);
arrival.time = cspice_str2et(arrival.EpochStr);
ToF = arrival.time - departure.time;

% Initial position vector [km]
x1 = cspice_spkezr(departure.Label,departure.time,'J2000','NONE',primary.label);
r1 = x1(1:3);
v1_0 = x1(4:6);

% Final position vector [km]
r2 = cspice_spkpos(arrival.Label,arrival.time,'J2000','NONE',primary.label);      

% Initial guess as Hohmann transfer velocity
E2M.objective = @(x) objectiveFun(x,departure.time,ToF,primary.GM,r2,r1);
E2M.x0 = v1_0;
E2M.solver = 'fsolve';

% Set options for using finite difference method
optsFinDiff = optimoptions('fsolve');
optsFinDiff.Display='iter';
optsFinDiff.FiniteDifferenceType = 'central';

E2M.options = optsFinDiff;
clc
solFD=fsolve(E2M);

optsJacobian = optimoptions('fsolve');
optsJacobian.Display='iter';
optsJacobian.SpecifyObjectiveGradient = true;

E2M.options = optsJacobian;
solSTM=fsolve(E2M);


[A,P,E,ERROR,VI,VF,TPAR,THETA] = lambertMR(r1, r2, ToF, primary.GM,0,0,0 );

disp('The difference between the analytical lambert solver and the shooting one is:')
norm(solSTM - VI')
norm(solFD - VI')

%% Point 3 - delta_v minimizer (lambert solver) %%%%%%%%%
clear E2MT primary departure arrival

departure.Label = 'Earth';
arrival.Label = 'Mars';
frame = 'ECLIPJ2000';

primary.Label = 'Sun';

% Time conversion constant
DAY2SECS = 24*3600;

departure.lb = cspice_str2et('2022-Jan-1 00:00:00.0000 TDB')/DAY2SECS;
departure.ub = cspice_str2et('2023-Jan-1 00:00:00.0000 TDB')/DAY2SECS;

arrival.lb = cspice_str2et('2022-Mar-1 00:00:00.0000 TDB')/DAY2SECS;
arrival.ub = cspice_str2et('2023-Dec-31 00:00:00.0000 TDB')/DAY2SECS;

minToF = 300;

% Define inital guess for the state (from computing DV grid) 
departureguess = cspice_str2et('2022-Jan-20 00:00:00.0000 TDB')/DAY2SECS;
arrivalguess = departureguess + minToF;

y0 = [departureguess/DAY2SECS; arrivalguess];

% Define optimization problem structure
E2MT.objective = @(y) lambertcostFcn(y); 
E2MT.x0 = y0; % Initial guess for the state
E2MT.A = [1 -1];                    % Linear inequality constr for time
E2MT.B = 0;
E2MT.lb = [departure.lb; 
            arrival.lb];
        
E2MT.ub = [departure.ub;
            arrival.ub];
        
E2MT.solver='fmincon';

% Define and set options
opts = optimoptions('fmincon');
opts.ConstraintTolerance = 0.1;
opts.Display = 'iter';
opts.PlotFcn = 'optimplotfval';
opts.ScaleProblem = true;

E2MT.options = opts;

% Compute solution
[sol,DeltaV] = fmincon(E2MT);
[DV,xx] = lambertcostFcn(sol);

% Checks
% End state position error
departureTime = sol(1)*DAY2SECS;
arrivalTime = sol(2)*DAY2SECS;

planetDeparture = cspice_spkpos(departure.Label,departureTime,frame,'NONE',primary.Label);
planetArrival = cspice_spkpos(arrival.Label,arrivalTime,frame,'NONE',primary.Label);

% Plot SC trajectory
figure()
plot3(xx(:,1),xx(:,2),xx(:,3),'LineWidth',2)
hold on
grid on
axis equal
xlabel('[km]')
ylabel('[km]')
zlabel('[km]')

% Plot planets at arrival and departure
plot3(planetDeparture(1),planetDeparture(2),planetDeparture(3),...
    '-o','Color','b','MarkerSize',10,'MarkerFaceColor','#D9FFFF');
plot3(planetArrival(1),planetArrival(2),planetArrival(3),...
     '-o','Color','r','MarkerSize',10,'MarkerFaceColor','#D95319');
legend('Transfer arc', 'Earth @ departure','Mars @ arrival')

% Print departure and arrival dates and time of flight
departureDate = cspice_et2utc(departureTime,'C',0)
arrivalDate = cspice_et2utc(arrivalTime,'C',0)
ToF_DAYS = diff(sol)

%% delta_v minimizer (free initial state x1) %%%%%%%%%
clear primary departure arrival

primary.Label = 'Sun';
primary.GM = cspice_bodvrd(primary.Label,'GM',1);   % Gravitational param

departure.Label = 'Earth';
arrival.Label = 'Mars';
frame = 'ECLIPJ2000';

% Time conversion constant
DAY2SECS = 24*3600;

departure.lb = cspice_str2et('2022-Jan-1 00:00:00.0000 TDB')/DAY2SECS;
departure.ub = cspice_str2et('2023-Jan-1 00:00:00.0000 TDB')/DAY2SECS;

arrival.lb = cspice_str2et('2022-Mar-1 00:00:00.0000 TDB')/DAY2SECS;
arrival.ub = cspice_str2et('2023-Dec-31 00:00:00.0000 TDB')/DAY2SECS;


% Define inital guess for the state (from computing DV grid) 
departureguess = cspice_str2et('2022-Sep-09 00:00:00.0000 TDB');

ToF = 300; %  ToF in days
arrivalguess = departureguess/DAY2SECS + ToF;

stateDep = cspice_spkezr(departure.Label,departureguess,frame,'NONE',primary.Label);
stateDep(4:6) = stateDep(4:6) + cross([0;0;1],stateDep(1:3)/norm(stateDep(1:3)))*7;

y0 = [stateDep; departureguess/DAY2SECS; arrivalguess];

% Define optimization problem structure
E2MT.objective = @(y) costFcn(y,departure,arrival,primary,frame); 
E2MT.x0 = y0; % Initial guess for the state
E2MT.A = [0 0 0 0 0 0 1 -1];                    % Linear inequality constr for time
E2MT.B = 0;
E2MT.lb = [-5e7; -5e7; -5e4; 
            -28; -28; -10;
            departure.lb; 
            arrival.lb];
        
E2MT.ub = [3e8; 3e8; 5e4;
            40; 40; 10;
            departure.ub;
            arrival.ub];
        
E2MT.nonlcon = @(y) enforcePositions(y,departure,arrival,primary,frame);
E2MT.solver='fmincon';

% Define and set options
opts = optimoptions('fmincon');
opts.ConstraintTolerance = 0.1;
opts.Display = 'iter';
opts.PlotFcn = 'optimplotfval';
opts.ScaleProblem = true;

E2MT.options = opts;

% Compute solution
[sol,DeltaV] = fmincon(E2MT);

% Checks
% End state position error
departureTime = sol(7)*DAY2SECS;
arrivalTime = sol(8)*DAY2SECS;
departureState = sol(1:6);

positionError = flow2BP(departureTime,departureState,arrivalTime,primary.GM,1)...
    - cspice_spkpos(arrival.Label,arrivalTime,frame,'NONE',primary.Label)
departureDate = cspice_et2utc(departureTime,'C',0)
arrivalDate = cspice_et2utc(arrivalTime,'C',0)

x0 = sol(1:6);

% Set options for ODE solver
opts = odeset('reltol', 1e-12, 'abstol', 1e-12);

[tt,xx] = ode113(@twobodyode,[departureTime arrivalTime],x0,opts,primary.GM);
planetDeparture = cspice_spkpos(departure.Label,departureTime,frame,'NONE',primary.Label);
planetArrival = cspice_spkpos(arrival.Label,arrivalTime,frame,'NONE',primary.Label);

figure()
plot3(xx(:,1),xx(:,2),xx(:,3),'LineWidth',2)
hold on
grid on
axis equal
xlabel('[km]')
ylabel('[km]')
zlabel('[km]')

plot3(planetDeparture(1),planetDeparture(2),planetDeparture(3),...
    '-o','Color','b','MarkerSize',10,'MarkerFaceColor','#D9FFFF');
plot3(planetArrival(1),planetArrival(2),planetArrival(3),...
     '-o','Color','r','MarkerSize',10,'MarkerFaceColor','#D95319');
legend('Transfer arc', 'Earth @ departure','Mars @ arrival','Location','best')

 % Clear kernel pool
cspice_kclear();

%% Ex 2
clearvars; close all; clc

% Load kernels
cspice_furnsh('assignment01.tm');

% Define list of celestial bodies:
labels = {'Earth';
          'Sun';
          'Moon'};

% Initialize propagation data
bodies = nbody_init(labels);

% Adimensionalizing constants
DAY2SECS = 24*3600;

% select integration frame string (SPICE naming convention)
frame = 'J2000';

dep_ep1 = '2022-Jan-01 00:00:00.0000 TDB';
dep_ep2 = '2022-Dec-31 00:00:00.0000 TDB';

t1_start = cspice_str2et(dep_ep1);
t1_end = cspice_str2et(dep_ep2);

% [y0,xxG1,xxG2] = biellipticGuess(t1_end,bodies,frame);
% figure()
% plot3(xxG1(1,:),xxG1(2,:),xxG1(3,:))
% hold on
% plot3(xxG2(1,:),xxG2(2,:),xxG2(3,:))

[y0,xxG,DVguess] = initialGuess2(t1_end,bodies,frame);
figure()
plot3(xxG(1,:),xxG(2,:),xxG(3,:))
hold on

grid on
axis equal
moonR = cspice_spkpos('moon',y0(end)*DAY2SECS,frame,'NONE','EARTH');
rL2 = moonR/norm(moonR)*(norm(moonR) + 60e3);

% Plot L2 target position
plot3(rL2(1),rL2(2),rL2(3),'*r')

% Plot initial guess positions
plot3(y0(1),y0(2),y0(3),'*k')
plot3(y0(7),y0(8),y0(9),'*k')
plot3(y0(13),y0(14),y0(15),'*k')
plot3(y0(19),y0(20),y0(21),'*k')

times = linspace(t1_start,t1_end,24);
for i=1:length(times)
    t1 = times(i);
%     [y0,xxG1,xxG2] = biellipticGuess(t1,bodies,frame);
    [y0,xxG1,xxG2] = initialGuess2(t1,bodies,frame);

    [DVguess1,DV1,DV2,DV3] = totalDV(y0,bodies);
    
    % Create problem structure
    
    L2.objective = @(y) totalDV(y,bodies);
    L2.x0 = y0;
    
    L2.A = [zeros(1,24) 1 0 -1;
        zeros(1,24) 1 -1 0;
        zeros(1,24) 0 1 -1];
    L2.B = [0;0;0];
    
    L2.Aeq = [zeros(1,24) 1 0 0;
        0 0 0*1 zeros(1,24)];
    L2.Beq = [t1/DAY2SECS; 0];
    
    L2.nonlcon = @(y) constraints(y,bodies);
    L2.solver = 'fmincon';
    
    % Define and set options
    opts = optimoptions('fmincon');
    opts.ConstraintTolerance = 0.01;
    opts.StepTolerance = 0.001;
    opts.Display = 'iter';
    opts.PlotFcn = 'optimplotfval';
    opts.Algorithm = 'active-set';
%     opts.OutputFcn = @(x,optimOptions,state) plotTrajectory(x,bodies);

    opts.MaxFunctionEvaluations = 50e3;
    opts.ScaleProblem = true;
    L2.options = opts;
    
    [sol,DeltaV(i),~,output,lambda,gradient,H] = fmincon(L2);
    figure(i)
    [tt,yy] = plotTrajectory(sol,bodies);
    [~,DV1,DV2,DV3] = totalDV(y0,bodies);
end

figure()
plot(times,DeltaV,'LineWidth',2)
ylabel('DV [km/s]')
grid on

cspice_kclear();

%% Ex 3 - OKAY
clearvars; close all; clc

% Load kernels
cspice_furnsh('assignment01.tm');

% Define boundary conditions
r0 = [0 -29597.43 0]';
v0 = [1.8349 0.0002 3.1783]';
m0 = 735;

rf = [0 -29617.43 0]';
vf = [1.8371 0.0002 3.1755]';

% Define S/C and initial orbit parameters
Tmax = 100*1e-6;
Isp = 3000;
t0 = 0;
GM = 398600;
epsilon0 = 0.5*norm(v0)^2 - GM/norm(r0);
semiMajorAxis = - GM/(2*epsilon0);
T = 2*pi*sqrt(semiMajorAxis^3/GM);

DAY2SECS = 24*60*60;

% Initial guess for the unknowns
Lambda0_guess = [ 5 5 5 1e4 1e4 1e4 1 T/DAY2SECS]';

% fsolve options
options = optimoptions('fsolve');
options.Display = 'iter';
options.MaxFunctionEvaluations = 3000;
options.MaxIterations = 500;

% Initial state vector
x0 = [r0; v0; m0];

% Solve
[Z,fval] = fsolve(@(x) shootingFcn(x,x0,t0,rf,vf,GM,Tmax,Isp),Lambda0_guess,options);
[nope, SOL] = shootingFcn(Z,x0,t0,rf,vf,GM,Tmax,Isp);

% Plot powered orbit of S/C transfer

plot3(SOL.y(1,:),SOL.y(2,:),SOL.y(3,:),'LineWidth',2)
hold on
plot3(SOL.y(1,1),SOL.y(2,1),SOL.y(3,1),'*r')
plot3(SOL.y(1,end),SOL.y(2,end),SOL.y(3,end),'*k')
xlabel('[km]')
ylabel('[km]')
zlabel('[km]')

plot3(0,0,0,'-o','Color','b','MarkerSize',10,'MarkerFaceColor','#D9FFFF')
axis equal padded
grid on

% Plot lambdaM (should always be >0 for time-optimal problems)
figure()
plot(SOL.x,SOL.y(14,:),'LineWidth',2)
% title('$\lambda_m$','Interpreter','latex')
xlabel('$Time [s]$','Interpreter','latex')
ylabel('$\lambda_m$','Interpreter','latex')
grid on

% Plot control action
figure()
u=zeros(size(SOL.x));
st=u;
for j = 1:length(SOL.x)
    m = SOL.y(7,j);
    LambdaV = SOL.y(11:13,j);
    LambdaM = SOL.y(14,j);
    [ui,sti] = thrustMagRatio(LambdaV,LambdaM,Isp,m);
    u(j)=ui;
    st(j)=sti;
end

plot(SOL.x,u,'LineWidth',2)
grid on
xlabel('$Time [s]$','Interpreter','latex')
ylabel('$Thrust\; u [-]$','Interpreter','latex')

% Compute the solution endpoint defect
norm(SOL.y(1:6,end)- [rf;vf],inf)

% Compute t_f = f(Tmax)
Tmax = (100:50:300)*1e-6;
for i = 1:length(Tmax)
    [Z(:,i),fval] = fsolve(@(x) shootingFcn(x,x0,t0,rf,vf,GM,Tmax(i),Isp),Lambda0_guess,options);
end
times = Z(end,:);

% Plot the transfer time
figure()
plot(Tmax*1e3,times*24,'LineWidth',2)
grid on
xlabel('$Thrust [N]$','Interpreter','latex')
ylabel('$Time of Flight [h]$','Interpreter','latex')

% Clear kernel pool
cspice_kclear();

%% Functions
function dy = twobodyode(~,y,mu)
%TWOBODYODE Restricted two body problem ODE function
% 
% PROTOTYPE:
%   dy = twobodyode(t,y)
% 
% INPUT:
%   t[1]       Time                                             [T]
%   y[6x1]     State of the system [r v]
%              (position and velocity vectors)                  [L, L/T]
%   mu[1]      Planetary constant
%
% OUTPUT:
%   dy[6x1]         Derivative of the state [L, L/T]
%
% CONTRIBUTORS:
%   Davide Iafrate
%
% VERSIONS
%   2020-09-24: First version
%

% Calculate radius
r = norm(y(1:3));

% Set the derivatives of the state
dy = [ y(4); 
       y(5); 
       y(6);
        - mu*y(1)/r^3;
        - mu*y(2)/r^3;
        - mu*y(3)/r^3; ];
end

function [xt,traj] = flow2BP(t0,x0,t,mu,posOutOnly)
%FLOW2BP Restricted two body problem flow function
% 
% PROTOTYPE:
%   xt=flow2BP(x0,t0,t,mu)
% 
% INPUT:
%   x0[6x1]     Initial state of the system
%   t0[1]       Initial time
%   t[1]        Final time                                           
%   y[6x1]      State of the system [r v]
%               (position and velocity vectors)                  
%   mu[1]       Planetary constant
%
% OUTPUT:
%   xt[6x1]         Final state
%
% CONTRIBUTORS:
%   Davide Iafrate
%
% VERSIONS
%   2020-09-24: First version
%

tspan = [t0 t];

% Set options for ODE solver
opts = odeset('reltol', 1e-12, 'abstol', 1e-12);

% Solve the ode
[~,xx] = ode113(@twobodyode,tspan,x0,opts,mu);

if nargin>4
    if posOutOnly == 1
        xt = xx(end,1:3)';
    else
        xt = xx(end,:)';
    end
else
    xt = xx(end,:)';
end

if nargout>1
    traj = xx;
end

end

function [g,J] = objectiveFun(vi,t1,ToF,mu,r2,r1)
    g = flow2BP(t1,[r1; vi],t1+ToF,mu,1) - r2;
  
    % Output the jacobian if required
    
    if nargout > 1
          x0 = [r1;vi];
          STM = stateTransitionMatrix(t1,x0,t1+ToF,mu);
          PHIrv = STM(1:3,4:6);
          
          J = PHIrv;
    
    end
end

function DV = costFcn(y,departureObj,arrivalObj,primary,frame)
%COSTFCN function that computes the total DELTA-V required by the
%trajectory

% Time conversion constant
DAY2SECS = 24*3600;

r1 = y(1:3);
v1 = y(4:6);

x1 = [r1; v1];

t1 = y(7)*DAY2SECS;
t2 = y(8)*DAY2SECS;

x2 = flow2BP(t1,x1,t2,primary.GM);

v2 = x2(4:6);

stateDep = cspice_spkezr(departureObj.Label,t1,frame,'NONE',primary.Label);
stateArrival = cspice_spkezr(arrivalObj.Label,t2,frame,'NONE',primary.Label);

vDep = stateDep(4:6);
vArr = stateArrival(4:6);

dv1 = norm(v1 - vDep);
dv2 = norm(v2 - vArr);
DV = dv1 + dv2;
end

function [DV,xx] = lambertcostFcn(y)
%COSTFCN function that computes the total DELTA-V required by the
%trajectory solving the Lambert

% Time conversion constant
DAY2SECS = 24*3600;

t1 = y(1)*DAY2SECS;
t2 = y(2)*DAY2SECS;

primary.label = 'Sun';
primary.GM = cspice_bodvrd(primary.label,'GM',1);   % Gravitational param
            
% select frame string (SPICE naming convention)
frame = 'ECLIPJ2000';

% select departure and arrival strings
departure.Label = 'Earth';
arrival.Label = 'Mars';

ToF = t2 - t1;

% Initial position vector [km]
x1 = cspice_spkezr(departure.Label,t1,'ECLIPJ2000','NONE',primary.label);
r1 = x1(1:3);
v1_0 = x1(4:6);

% Final position vector [km]
r2 = cspice_spkpos(arrival.Label,t2,'ECLIPJ2000','NONE',primary.label);      

% Initial guess as planet velocity
E2M.objective = @(x) objectiveFun(x,t1,ToF,primary.GM,r2,r1);
E2M.x0 = v1_0;
E2M.solver = 'fsolve';

% Set options using STM
opts = optimoptions('fsolve');
opts.Display='iter';
opts.SpecifyObjectiveGradient = true;
opts.FunctionTolerance = 0.01;
E2M.options = opts;
v1=fsolve(E2M);

x1 = [r1; v1];

% Propagate to get the final solution
[xt,xx] = flow2BP(t1,x1,t2,primary.GM);

stateDep = cspice_spkezr(departure.Label,t1,frame,'NONE',primary.label);
stateArrival = cspice_spkezr(arrival.Label,t2,frame,'NONE',primary.label);

vDep = stateDep(4:6);
vArr = stateArrival(4:6);

v2 = xt(4:6);

dv1 = norm(v1 - vDep);
dv2 = norm(v2 - vArr);
DV = dv1 + dv2;

end

function [c, ceq] = enforcePositions(y,departureObj,arrivalObj,primary,frame)
% Time conversion constant
DAY2SECS = 24*3600;

r1 = y(1:3);
v1 = y(4:6);

t1 = y(7)*DAY2SECS;
t2 = y(8)*DAY2SECS;

radiusDep = cspice_spkpos(departureObj.Label,t1,frame,'NONE',primary.Label);
radiusArrival = cspice_spkpos(arrivalObj.Label,t2,frame,'NONE',primary.Label);

x2 = flow2BP(t1,[r1; v1],t2,primary.GM);
r2 = x2(1:3);

% Compute constraint violations
c1 = norm(r1-radiusDep);
c2 = norm(r2-radiusArrival);

ceq = [c1; c2];
c = [];
end

function STM = stateTransitionMatrix(t0,x0,t,mu)
%STM Summary of this function goes here
%   Detailed explanation goes here

stateLength = length(x0);

function dy = stmDynamics(t,Y)
    X = Y(1:6);
    x = X(1);
    y = X(2);
    z = X(3);
    r = norm(X(1:3));
    
    subMatrix = [x^2 x*y x*z;
                 x*y,y^2,y*z;
                 x*z,y*z,z^2];
    
    A = [zeros(3) eye(3);
        3*mu/r^5*subMatrix-mu/r^3*eye(3) zeros(3)];
    
    dX = twobodyode(t,X,mu);
    PHIvec = Y(stateLength+1:stateLength+stateLength^2);
    PHI = reshape(PHIvec,[stateLength stateLength]);
    dPHI = A*PHI;
    dPHIvec = dPHI(:);
    dy = [dX; dPHIvec];
end

STM0 = eye(stateLength);
y0 = [x0;STM0(:)];
opts = odeset('Reltol',1e-13,'AbsTol',1e-14);

sol = ode113(@stmDynamics,[t0 t],y0,opts);
endState = sol.y(:,end);
STMvec = endState(stateLength+1:stateLength+stateLength^2);
STM = reshape(STMvec,[stateLength stateLength]);

end

function [eqConstr, sol] = shootingFcn(unknown,x0,t0,rf,vf,mu,Tmax,Isp)
%SHOOTINGFCN    shooting function for the Time optimal TPBVP for continous
%               thrust guidance
% 
% INPUT:
%   unknown[8x1]    unknowns of the problem
%                       |-> Costate initial value Lambda0
%                       |-> final time tf
% 
%   x0[7x1]         initial state of the S/C
%   rf[3x1]         target position
%   vf[3x1]         target velocity
%   mu[1]           gravitational parameter of the primary
%   Tmax[1]         maximum thrust of the S/C
%   Isp[1]          specific impulse of the S/C
% 
% OUTPUT:
%   eqConstr[8x1]   constraint values
DAY2SECS = 24*60*60;
Lambda0 = unknown(1:7);
tf = unknown(8);            % IN DAYS!
tf = tf*DAY2SECS;
g0 = 9.80665*1e-3;

% Integrate dynamics
Y0 = [x0; Lambda0];
options = odeset('reltol', 3e-14, 'abstol', 1e-14);

SOL = ode113(@(t,Y) TPBVP_2BP(t,Y,mu,Tmax,Isp),[t0 tf],Y0,options);

% Retrieve the solution at final time
YF=SOL.y(:,end);
XF = YF(1:7);
LambdaF = YF(8:14);

r = XF(1:3);     LambdaR = LambdaF(1:3);
v = XF(4:6);     LambdaV = LambdaF(4:6);
m = XF(7);       LambdaM = LambdaF(7);

% Compute the Hamiltonian at final time
[u,St] = thrustMagRatio(LambdaV,LambdaM,Isp,m);
Hf = 1 + dot(LambdaR,v) - mu/norm(r)^3*dot(LambdaV,r) +...
        + Tmax*u*St/Isp/g0;

% Compute the constraints to be zeroed
errorH = Hf;
errorPosition = XF(1:3) - rf;
errorVelocity = XF(4:6) - vf;

eqConstr = [errorPosition; errorVelocity; LambdaM*100; errorH];

% Output solution of ODE if requested
% if nargout>2
    sol = SOL;
% end

end

function dY = TPBVP_2BP(t,Y,mu,Tmax,Isp)
%TPBVP_2BP  Two point BVP differential equations for the Low Thrust
%           controlled S/C problem
% 
% PROTOTYPE:
%   dY = TPBVP_2BP(Y,mu,Tmax,Isp)
% 
% INPUT:
%   Y[14x1]     S/C state and costate 
%   mu[1]       gravitational constant of the primary
%   Tmax[1]     maximum thrust of the S/C
%   Isp[1]      specific impulse of the S/C
% 
% OUTPUT:
%   dY[14x1]    dynamics of S/C state and of the costate

% Extract the state and costate
X = Y(1:7);
Lambda = Y(8:14);

r = X(1:3);     LambdaR = Lambda(1:3);
v = X(4:6);     LambdaV = Lambda(4:6);
m = X(7);       LambdaM = Lambda(7);

[u,~] = thrustMagRatio(LambdaV,LambdaM,Isp,m);

% Compute the unitary PRIMER VECTOR
alfa = - LambdaV/norm(LambdaV);

% Compute the derivatives of the S/C state
dX = EOM(t,X,u,alfa,mu,Tmax,Isp);

% Compute derivatives of the costate vector
dLambdaR = -3*mu/norm(r)^5*dot(r,LambdaV)*r + mu/norm(r)^3*LambdaV;
dLambdaV = -LambdaR;
dLambdaM = -u*norm(LambdaV)*Tmax/m^2;

% Output the derivative of the state and costate
dY = [dX; dLambdaR; dLambdaV; dLambdaM];
end

function [u,St] = thrustMagRatio(LambdaV,LambdaM,Isp,m)
% Generate the control input magnitude and direction
g0 = 9.80665*1e-3;  % [km/s^2]

% Switching function computation
St = - norm(LambdaV)*Isp*g0/m - norm(LambdaM);

if St >= 0
    u = 0;
else
    u = 1;
end

end

function dX = EOM(t,X,u,alfa,mu,Tmax,Isp)
%EOM equations of motion for the controlled 2BP
% 
% PROTOTYPE:
%   dX = EOM(t,X,u,params)
% 
% INPUT:
%   t[1]        time instant
%   X[7x1]      state vector [r(3); v(3); m(1)]'
%   u[1x1]      thrust magnitude control
%   alfa[3x1]   unit vector of thrust direction 
%   mu[1]       gravitational constant of the primary
%   Tmax[1]     maximum thrust of the S/C
%   Isp[1]      specific impulse of the S/C
% 
% OUTPUT:
%   dX[7x1] derivative of the state

% Define standard gravitational acceleratio
g0 = 9.80665*1e-3;       % [km/s^2] IMPROVEMENT: RETRIEVE FROM SPICE

% extract the state parameters
r = X(1:3);
v = X(4:6);
m = X(7);

% Thrust pointing unitary vector

dr = v;
dv = -mu/norm(r)^3*r + u*Tmax/m*alfa;      % Convert thrust term to km/s^2
dm = -norm(u)*Tmax/(Isp*g0);

dX = [dr; dv; dm];
end

function [bodies] = nbody_init(labels)
%NBODY_INIT Initialize planetary data for n-body propagation
%   Given a set of labels of planets and/or barycentres, returns a
%   cell array populated with structures containing the body label and the
%   associated gravitational constant.
%
%
% Author
%   Name: ALESSANDRO 
%   Surname: MORSELLI
%   Research group: DART
%   Department: DAER
%   University: Politecnico di Milano 
%   Creation: 26/09/2021
%   Contact: alessandro.morselli@polimi.it
%   Copyright: (c) 2021 A. Morselli, Politecnico di Milano. 
%                  All rights reserved.
%
%
% Notes:
%   This material was prepared to support the course 'Satellite Guidance
%   and Navigation', AY 2021/2022.
%
%
% Inputs:
%   labels : [1,n] cell-array with object labels
%
% Outputs:
%   bodies : [1,n] cell-array with struct elements containing the following
%                  fields
%                  |
%                  |--bodies{i}.name -> body label
%                  |--bodies{i}.GM   -> gravitational constant [km**3/s**2]
%
%
% Prerequisites:
%   - MICE (Matlab SPICE)
%   - Populated kernel pool (PCK kernels)
%

% Initialize output
bodies = cell(size(labels));

% Loop over labels
for i = 1:length(labels)
    % Store body label
    bodies{i}.name = labels{i};
    % Store body gravitational constant
    bodies{i}.GM   = cspice_bodvrd(labels{i}, 'GM', 1);
end

end

function [dxdt] = nbody_rhs(t, x, bodies, frame)
%NBODY_RHS Evaluates the right-hand-side of a N-body propagator
%   Evaluates the right-hand-side of a newtonian N-body propagator.
%   The integration centre is the Earth 
%
%
% Author
%   Name: ALESSANDRO 
%   Surname: MORSELLI
%   Research group: DART
%   Department: DAER
%   University: Politecnico di Milano 
%   Creation: 26/09/2021
%   Contact: alessandro.morselli@polimi.it
%   Copyright: (c) 2021 A. Morselli, Politecnico di Milano. 
%                  All rights reserved.
%
%
% Notes:
%   This material was prepared to support the course 'Satellite Guidance
%   and Navigation', AY 2021/2022.
%
%
% Inputs:
%   t      : [1,1] ephemeris time (ET SPICE), seconds past J2000 (TDB)
%   x      : [6,1] cartesian state vector wrt Earth-Barycentre
%   bodies : [1,6] cell-array created with function nbody_init
%
% Outputs:
%   dxdt   : [6,1] RHS, newtonian gravitational acceleration only
%
% Prerequisites:
%   - MICE (Matlab SPICE)
%   - Populated kernel pool (SPK, LSK, PCK kernels)
%

% Initialize right-hand-side
dxdt = zeros(6,1);

% Position derivative is object's velocity
dxdt(1:3) = x(4:6);

% Extract the object position at state x
r = x(1:3);

% Compute the primary acceleration
normRR2 = dot(r,r);
normRR = sqrt(normRR2);

dxdt(4:6) = -bodies{1}.GM*r/(normRR*normRR2);

% Compute perturbing accelerations from the other bodies

for i=2:length(bodies)

    % Retrieve position and velocity of i-th celestial body wrt Earth
    % Barycentre in an inertial frame
    rho_body_es = cspice_spkezr(bodies{i}.name, t, frame, 'NONE', 'EARTH');
    
    % Extract i-th body pos wrt Earth
    rho_body = rho_body_es(1:3);
    rho2 = dot(rho_body,rho_body);
    rho = sqrt(rho2);
    
    % Extract object position wrt. i-th celestial body
    dd_sc_body = r - rho_body;
    
    % Compute square distance and distance
    dist2 = dot(dd_sc_body, dd_sc_body);
    dist = sqrt(dist2);

    % Compute the gravitational acceleration using Newton's law
    aa_grav =  - bodies{i}.GM * (dd_sc_body /(dist*dist2) + rho_body/(rho*rho2));

    % Sum up acceleration to right-hand-side
    dxdt(4:6) = dxdt(4:6) + aa_grav;

end

end

function xt = flowNBP(t0,x0,t,bodies,frame,posOutOnly)
%FLOWNBP Restricted N body problem flow function
% 
% PROTOTYPE:
%   xt=flow2BP(x0,t0,t,mu)
% 
% INPUT:
%   x0[6x1]     Initial state of the system
%   t0[1]       Initial time
%   t[1]        Final time                                           
%   y[6x1]      State of the system [r v]
%               (position and velocity vectors)                  
%   mu[1]       Planetary constant
%
% OUTPUT:
%   xt[6x1]         Final state
%
% CONTRIBUTORS:
%   Davide Iafrate
%
% VERSIONS
%   2020-09-24: First version
%

tspan = [t0 t];

options = odeset('reltol', 1e-8, 'abstol', 1e-8);
[tt, xx] = ode113(@(t,x) nbody_rhs(t,x,bodies,frame), tspan, x0, options);


if nargin > 5
    if posOutOnly == 1
        xt = xx(end,1:3)';
    else
        xt = xx(end,:)';
    end
else
    xt = xx(end,:)';
end

end

function [c, ceq] = constraints(y,bodies)
%CONSTRAINTS function producing the equality constraints and the
%inequality constraints vectors
% 
% INPUT:
%   y[29x1] :       state vector of optimization variables
%   bodies[1,6] :   cell-array created with function nbody_init
% 
% OUTPUT:
%   c[4x1]  :       nonlinear inequality constraints
%   ceq[7x1] :      nonlinear equality constraints

% set reference frame
frame = 'J2000';

% Extract the state vectors
x1 = y(1:6);    r1 = x1(1:3);   v1 = x1(4:6);
x2 = y(7:12);   r2 = x2(1:3);   v2 = x2(4:6);
x3 = y(13:18);  r3 = x3(1:3);   v3 = x3(4:6);
x4 = y(19:24);  r4 = x4(1:3);   v4 = x4(4:6);

% Parameters used in calculations
muE = bodies{1}.GM;
radiiE = cspice_bodvrd( 'EARTH', 'RADII', 3);
R_e = radiiE(1);
DAY2SECS = 24*3600;

% First burn time
t1 = y(25)*DAY2SECS;

% Second burn time
t3 = y(26)*DAY2SECS;

% Last burn times
tN = y(27)*DAY2SECS;

% Perform time discretization of each arc in 2 segments
N = 3;
temp = linspace(t3,tN,N);
t = [linspace(t1,t3,N) temp(2:end)];


% Retrieve final orbital position of moon from ephemeris
xMoon_f = cspice_spkezr('moon',tN,frame,'NONE','EARTH');

% Find L2 final state, same velocity as the moon but radius 60000km greater
rL2_f = xMoon_f(1:3)/norm(xMoon_f(1:3))*(norm(xMoon_f(1:3)) + 60000);

% Compute the constraints
ceq = [norm(r1) - (200+R_e);
%        dot(v1,r1)*100;
       flowNBP(t(1),x1,t(2),bodies,frame) - x2;
       flowNBP(t(2),x2,t3,bodies,frame,1) - r3;
       flowNBP(t3,x3,t(4),bodies,frame) - x4;
       flowNBP(t(4),x4,tN,bodies,frame,1) - rL2_f];
   
c = [500e3 - norm(r3);
    norm(r3)- 1e6;
    norm(v1) - sqrt(2*muE/norm(r1))
    sqrt(muE/norm(r1)) - norm(v1)];

end

function [DV,DV1,DV2,DV3] = totalDV(y,bodies)
%TOTALDV cost function computing the total DeltaV required for the three
%maneuvers
% 
% INPUT:
%   y[29x1] :       state vector of optimization variables
%   bodies[1,6] :   cell-array created with function nbody_init
% 
% OUTPUT:
%   DV[1]   :         total DeltaV required by the three maneuvers
%   DV1[1]  :         DeltaV required by the first maneuver
%   DV2[1]  :         DeltaV required by the second maneuver
%   DV3[1]  :         DeltaV required by the third maneuver

frame = 'J2000';
% Extract the state vectors
x1 = y(1:6);    r1 = x1(1:3);   v1 = x1(4:6);
x2 = y(7:12);   r2 = x2(1:3);   v2 = x2(4:6);
x3 = y(13:18);  r3 = x3(1:3);   v3 = x3(4:6);
x4 = y(19:24);  r4 = x4(1:3);   v4 = x4(4:6);

% Parameters used in calculations
GM = bodies{1}.GM;
DAY2SECS = 24*3600;

% First burn time
t1 = y(25)*DAY2SECS;

% Second burn time
t3 = y(26)*DAY2SECS;

% Last burn times
t5 = y(27)*DAY2SECS;

% Perform time discretization of each arc in 2 segments
N = 3;
tvec1 = linspace(t1,t3,N);
tvec2 = linspace(t3,t5,N);

t = [tvec1 tvec2(2:end)];

% Retrieve final orbital position of moon from ephemeris
xMoon_f = cspice_spkezr('moon',t5,'J2000','NONE','EARTH');

% Find L2 final state, same velocity as the moon but radius 60000km greater
vL2_f = xMoon_f(4:6);

% Compute the flows
x3f = flowNBP(t(2),x2,t3,bodies,frame);
x5f = flowNBP(t(4),x4,t5,bodies,frame);

% Compute the circular velocity on the initial orbit
v0 = cross([0;0;1],r1/norm(r1))*sqrt(GM/norm(r1));

% Compute the DVs
DV1 = norm(v1 - v0);
DV2 = norm(v3 - x3f(4:6));
DV3 = norm(vL2_f - x5f(4:6));

DV = DV1 + DV2 + DV3;
end

function [y0,xxG1,xxG2] = biellipticGuess(t1,bodies,frame)
% Initial guess definition (using a 2BP bielliptic transfer as reference)
% 
%   INPUT:
%       t1[1]

opts = odeset('reltol', 1e-12, 'abstol', 1e-12);

% Parameters used in calculations
GM = bodies{1}.GM;
radiiE = cspice_bodvrd( 'EARTH', 'RADII', 3);
R_e = radiiE(1);
DAY2SECS = 24*3600;
r0 = R_e + 200;

% Get final position we need to reach
rp1 = 1/(2/(r0)-1/(500e3));

a1 = (rp1 + 500e3)/2;
e1 = (500e3-rp1)/(2*a1);
a2 = (500e3 + 400e3)/2;
T1 = 2*pi*sqrt(a1^3/GM);
T2 = 2*pi*sqrt(a2^3/GM);

% Get the two eccentric anomalies
th1 = pi/2;
th2 = pi;

E1 = 2*atan2(tan(th1/2),sqrt((1+e1)/(1-e1)));
E2 = 2*atan2(tan(th2/2),sqrt((1+e1)/(1-e1)));

DeltaT = sqrt(a1^3/GM)*(E2-e1*sin(E2)-E1+e1*sin(E1));
% Compute Delta t
t3 = t1 + DeltaT;
t5 = t3 + T2/2;

moonX = cspice_spkezr('moon',t5,frame,'NONE','EARTH');
moonR = moonX(1:3);
moonV = moonX(4:6);

[aM, eM, iM, RAANM, omegaM, fM] = car2kep(moonR,moonV,GM);
rL2 = moonR/norm(moonR)*(norm(moonR)+60e3);
r3 = -moonR/norm(moonR)*500e3;


% Orbital periods of the two transfer orbits
th1 = pi/2;
th2 = pi;
RAAN1 = RAANM;
om1 = fM+omegaM;


[r1,v1] = kep2car(a1,e1,iM,RAAN1,om1,th1,GM);
[r2,v2] = kep2car(a1,e1,iM,RAAN1,om1,th2,GM);

RAAN2 = RAANM;
om2 = om1;
a2 = (norm(r3)+norm(rL2)+100e3)/2;
e2 = -(norm(r3)-100e3-norm(rL2))/(2*a2);

th3 = pi;
th4 = 3*pi/2;

[r3,v3] = kep2car(a2,e2,iM,RAAN2,om2,th3,GM);
[r4,v4] = kep2car(a2,e2,iM,RAAN2,om2,th4,GM);

% Orbital period
T2 = 2*pi*sqrt(a2^3/GM);

% Propagation times
t3 = t1 + DeltaT;
t5 = t3 + T2/2;

x1 = [r1';v1']; x2_test = [r2';v2'];
x3 = [r3';v3']; x4_test = [r4';v4'];

% Propagate first orbit
% sol1 = ode113(@(t,x) nbody_rhs(t,x,bodies,frame),[t1 t3],x1,opts);

% Propagate second orbit
% sol2 = ode113(@(t,x) nbody_rhs(t,x,bodies,frame),[t3 t5],x3,opts);

% Propagate first orbit
sol1 = ode113(@(t,x) twobodyode(t,x,GM),[t1 t3],x1,opts);

% Propagate second orbit
sol2 = ode113(@(t,x) twobodyode(t,x,GM),[t3 t5],x3,opts);

% Extract points
x2 = deval(sol1,(t1+t3)/2);
x4 = deval(sol2,(t3+t5)/2);

% Output inital guess for the unknowns vector
y0 = [x1; 
      x2;
      x3;
      x4;
      t1/DAY2SECS;
      t3/DAY2SECS;
      t5/DAY2SECS];
  
 if nargout > 1
     xxG1 = sol1.y;
     xxG2 = sol2.y;
%      DV = norm(v0-VI1)+norm(VF1-VI2)+norm(VF2-moonV);
 end
 
end

function [y0,xxG,DV] = initialGuess(t1,bodies,frame)
% Initial guess definition (using a 2BP bielliptic transfer as reference)
% 
%   INPUT:
%       t1[1]

opts = odeset('reltol', 1e-12, 'abstol', 1e-12);

% Parameters used in calculations
GM = bodies{1}.GM;
radiiE = cspice_bodvrd( 'EARTH', 'RADII', 3);
R_e = radiiE(1);
DAY2SECS = 24*3600;

% Semi-major axis of the two transfer orbits
a1 = (200 + 500e3 + R_e)/2;
a2 = (500e3 + (378+60)*1e3 )/2;

% Orbital periods of the two transfer orbits
T1 = 2*pi*sqrt(a1^3/GM);
T2 = 2*pi*sqrt(a2^3/GM);

t3 = t1 + T1/2 + 2*DAY2SECS;
t5 = t3 + T2/2;

% Get final position we need to reach
moonX = cspice_spkezr('moon',t5,frame,'NONE','EARTH');
moonR = moonX(1:3);
moonV = moonX(4:6);

% Initial and final position of the first arc in J2000 frame
r1 = (6378+200)*moonR/norm(moonR);
r3 = (500e3 + 60*1e3 + norm(moonR))*-(moonR+[10;10;10])/norm(moonR);
r5 = moonR/norm(moonR)*(norm(moonR)+60e3);


% Get initial and final velocity for the first arc
[A,P,E,ERROR,VI1,VF1,TPAR,THETA] = lambertMR(r1, r3, t3-t1 , GM,0,0,0 );
x1 = [r1; VI1'];

% Circular orbit velocity
v0 = VI1/norm(VI1)*sqrt(GM/norm(r1));

% Get initial and final velocities for the second arc
[A,P,E,ERROR,VI2,VF2,TPAR,THETA] = lambertMR(r3, r5, t5-t3 , GM,0,0,0 );
x3 = [r3; VI2'];

% Propagate the first orbit for half a period
segmentsNumber = 2;
times1 = linspace(t1,t3,segmentsNumber+1);
times2 = linspace(t3,t5 - 1*DAY2SECS,segmentsNumber+1);

sol1 = ode113(@(t,x) nbody_rhs(t,x,bodies,frame),[t1 t3],x1,opts);

% Propagate second orbit
sol2 = ode113(@(t,x) nbody_rhs(t,x,bodies,frame),[t3 t5],x3,opts);

% Extract points
x2 = deval(sol1,(t1+t3)/2);
x4 = deval(sol2,(t3+t5)/2);
% Output inital guess for the unknowns vector
y0 = [x1; 
      x2;
      x3;
      x4;
      t1/DAY2SECS;
      t3/DAY2SECS;
      t5/DAY2SECS];
  
 if nargout > 1
     xxG = [sol1.y sol2.y];
     DV = norm(v0-VI1)+norm(VF1-VI2)+norm(VF2-moonV);
 end
 
end

function [y0,xxG,DV] = initialGuess2(t1,bodies,frame)
% Initial guess definition (using a 2BP bielliptic transfer as reference)
% 
%   INPUT:
%       t1[1]

opts = odeset('reltol', 1e-12, 'abstol', 1e-12);

% Parameters used in calculations
GM = bodies{1}.GM;
radiiE = cspice_bodvrd( 'EARTH', 'RADII', 3);
R_e = radiiE(1);
DAY2SECS = 24*3600;

% Semi-major axis of the two transfer orbits
a1 = (200 + 505e3 + R_e)/2;
a2 = (505e3 + (378+60)*1e3 )/2;

% Orbital periods of the two transfer orbits
T1 = 2*pi*sqrt(a1^3/GM);
T2 = 2*pi*sqrt(a2^3/GM);

t3 = t1 + T1/2;
t5 = t3 + T2/2;

% Get final position we need to reach
moonX = cspice_spkezr('moon',t5,frame,'NONE','EARTH');
moonR = moonX(1:3);
moonV = moonX(4:6);
r5 = moonR/norm(moonR)*(norm(moonR)+60e3);

% Initial and final position of the first arc in J2000 frame
r1dir = [moonR(1:2); 0]/norm([moonR(1:2); 0]);
r1norm = 6378+200;
r1 = r1norm*r1dir;
v1 = sqrt(GM*(2/r1norm - 1/a1))*cross([0;0;1],r1dir)/norm(cross([0;0;1],r1dir));
x1 = [r1;v1];

sol1 = ode113(@(t,x) twobodyode(t,x,GM),[t1 t3],x1,opts);

r3 = sol1.y(1:3,end);

normal=cross(r3,r5);
normal = normal/norm(normal);

if normal(3)>0
    v3 =sqrt(GM*(2/505e3 - 2/(norm(r5)+505e3)))*cross(normal,r3/norm(r3));
else
    v3 =sqrt(GM*(2/505e3 - 2/(norm(r5)+505e3)))*cross(-normal,r3/norm(r3));
end
x3 = [r3;v3];
% Propagate second orbit
sol2 = ode113(@(t,x) twobodyode(t,x,GM),[t3 t5],x3,opts);

% Extract points
x2 = deval(sol1,(t1+t3)/2);
x4 = deval(sol2,(t3+t5)/2);
% Output inital guess for the unknowns vector
y0 = [x1; 
      x2;
      x3;
      x4;
      t1/DAY2SECS;
      t3/DAY2SECS;
      t5/DAY2SECS];
  
 if nargout > 1
     xxG = [sol1.y sol2.y];
     DV = 10;
%      DV = norm(v0-VI1)+norm(VF1-VI2)+norm(VF2-moonV);
 end
 
end

function [stop,tt,yy] = plotTrajectory(y,bodies)
% set reference frame
frame = 'J2000';
opts = odeset('reltol', 1e-12, 'abstol', 1e-12);

% Extract the state vectors
x1 = y(1:6);    r1 = x1(1:3);   v1 = x1(4:6);
x2 = y(7:12);   r2 = x2(1:3);   v2 = x2(4:6);
x3 = y(13:18);  r3 = x3(1:3);   v3 = x3(4:6);
x4 = y(19:24);  r4 = x4(1:3);   v4 = x4(4:6);

DAY2SECS = 24*3600;
% First burn time
t1 = y(25)*DAY2SECS;

% Second burn time
t3 = y(26)*DAY2SECS;

% Last burn times
tN = y(27)*DAY2SECS;

% Parameters used in calculations
[tt1,xx1] = ode113(@(t,x) nbody_rhs(t,x,bodies,frame),[t1 t3],x1,opts);
[tt2,xx2] = ode113(@(t,x) nbody_rhs(t,x,bodies,frame),[t3 tN],x3,opts);

yy = [xx1; xx2];
tt=[tt1;tt2];

stop = 0;

% plot solution
figure()
plot3(yy(:,1),yy(:,2),yy(:,3))
hold on
moonR = cspice_spkpos('moon',tt(end),frame,'NONE','EARTH');
rL2 = moonR/norm(moonR)*(norm(moonR)+60e3);
plot3(rL2(1),rL2(2),rL2(3),'*r')
grid on
axis equal
end

function [a, e, i, RAAN, omega, f] = car2kep(r_vec,v_vec,mu)
% Convert from state vector in cartesian coordinates to keplerian coordinates
% SPECIAL CASES: 
%
%       [i = 0]
%           the NODE LINE is taken as the reference direction Gamma (vernal
%           equinox
%           RAAN is 0
%           argument of periapsis omega is undefined (NaN)
%       
%       [e = 0 & i != 0]
%           f is calculated as the argument of latitude, angular
%           distance between the line of nodes and the position vector
%           omega is 0
%       [e = 0 & i = 0]
%           f is calculated as the true longitude, angular
%           distance between the gamma direction and the position vector
% PROTOTYPE:
%   [a, e, i, RAAN, omega, M] = kep2car(r_vec,v_vec) 
%   
% INPUT:
%
%   r_vec[3]   position vector     [ km/s ]
%   v_vec[3]   velocity vector     [ km/s ]
%   mu[1]      gravitational paramer [ km^3/s^2 ]
%
% OUTPUT:
%
%   a[1]        semi-major axis     [ km ]
%   e[1]        eccentricity        [ - ]
%   i[1]        inclination         [ rad ]
%   RAAN[1]     right ascension of the ascending node [ rad ]
%   omega[1]    argument of perigee [ rad ]
%   f[1]        true anomaly        [ rad ]
%
% CONTRIBUTORS:
%   Davide Iafrate
%
% VERSIONS
%   2020-10-12: First version
%

eps = 1.e-10;

% calculate magnitude of r_vec, v_vec

r = norm(r_vec);
v = norm(v_vec);

%% compute the radial component of velocity

v_radial = dot(v_vec,r_vec)/r;

%% compute the angular momentum vector

h_vec = cross(r_vec,v_vec);

%% compute the magnitude of the angular momentum vector

h = norm(h_vec);

%% compute the inclination of the orbit
i = acos(h_vec(3)/h);

%% compute the node line

N_vec = cross([0 0 1], h_vec);

n = norm(N_vec);

%% compute the RAAN     % solving the case i = 0

if i < eps
    RAAN = 0;
    
elseif (N_vec(2)>= 0)
    RAAN = acos(N_vec(1)/n);
    
else
    RAAN = 2*pi - acos(N_vec(1)/n);

end

%% compute the eccentricity vector and its magnitude

e_vec = (r_vec.*(v^2 - mu/r) - r.*v_radial.*v_vec)./mu;
e = norm(e_vec);


%% compute the argument of perigee omega

if e < eps
    
    omega = 0;

elseif (e_vec(3) >= 0)
    omega = acos(dot(N_vec,e_vec)/(n*e));

else
    omega = 2*pi - acos(dot(N_vec,e_vec)/(n*e));
    
end



%% compute the true anomaly f

% REFERENCE:
% https://en.wikipedia.org/wiki/True_anomaly
% used for special cases

if e < eps     
    if i < eps               % Output the true longitude, omega is set to 0:
        f = acos(r_vec(1) / r);
        
        if v_vec(1) > 0
            
            f = 2*pi - acos(r_vec(1) / r);
        
        end
        
    else
        
        % Output the argument of latitude
        % https://en.wikipedia.org/wiki/True_anomaly#Circular_orbit
        % omega is assumed as 0
        if r_vec(3) < 0
            
            f = 2*pi - acos(dot(N_vec , r_vec) / (n*r));
            
        else
            f = acos(dot(N_vec , r_vec) / (n*r));
        end
    end
elseif (v_radial >= 0)
    f = acos(dot(e_vec,r_vec)/(e*r));
else
    f = 2*pi - acos(dot(e_vec,r_vec)/(e*r));
end

%% semi-major axis (scattered among various lectures)

if e < eps % Circular orbits
    a = r;
elseif e == 1 % Parabolic orbit
    a = 'infinite';
else % Elliptical/Hyperbolic orbit
    a = (h^2/mu)/(1 - e^2);
end
end

function [r_vec,v_vec]=kep2car(a,e,i,RAAN,omega,f,mu)

% Solution to: [a,e,f,RAAN,omega,i,mu]~~>[r,v]
%
% ATTENTION: 
% 1)    Not sure it works with parabolic and hyperbolic orbits;
% 1)    If i=0, RAAN loses physical meaning while omega represents the angle
%       between X end the eccentricity vector. The algorithm gives strange 
%       results for any other value. The reason being that while the axis Z
%       is the same despite the value of RAAN & omega, X & Y axis vary;
% 2)    when e=0, omega has no physical meaning. The reason is the same as
%       before.
% PROTOTYPE:
%    [r_vec,v_vec] = kep2car(a, e, i, RAAN, omega, f, mu) 
%   
% INPUT:
%   a[1]        semi-major axis     [ km ]
%   e[1]        eccentricity        [ - ]
%   i[1]        inclination         [ rad ]
%   RAAN[1]     right ascension of the ascending node [ rad ]
%   omega[1]    argument of perigee [ rad ]
%   f[1]        true anomaly        [ rad ]
%   mu[1]      gravitational paramer [ km^3/s^2 ]
%
% OUTPUT:
%   r_vec[3]   position vector     [ km ]
%   v_vec[3]   velocity vector     [ km/s ]
%   
% CONTRIBUTORS:
%   Davide Demartini
%   Davide Iafrate
%
% VERSIONS
%   2020-10-12: First version
%

if i==0 && e==0
    f = f + RAAN + omega;
    RAAN = 0;
    omega = 0;


elseif i==0
    omega = omega + RAAN;
    RAAN = 0;


elseif e==0
    f = f + omega;
    omega = 0;
end


% Rotation Matrixes

% rotation around K of an angle RAAN

R_RAAN=[cos(RAAN)   sin(RAAN) 0
        -sin(RAAN)  cos(RAAN) 0
         0          0         1];

% rotation around the line of nodes of an angle i
R_i=[1 0 0;
    0 cos(i) sin(i)
    0 -sin(i) cos(i)]; 

% rotation around h of amount of an angle omega
R_omega=[cos(omega) sin(omega) 0
    -sin(omega) cos(omega) 0
    0 0 1]; % Perigee anomaly

% Trasformation matrix RF Perifocal ~~> Geocentric
% R = R_omega*R_i*R_RAAN  GE~~>PF
% T = R'                  PF~~>GE
T=(R_omega*R_i*R_RAAN)';

% Parameter p
p=a*(1-e^2);

% Compute position components in PF RF
r_pf_x=p*cos(f)/(1+e*cos(f));
r_pf_y=p*sin(f)/(1+e*cos(f));
r_pf_vect=[r_pf_x; r_pf_y; 0];

% Compute velocity components in PF RF
co=sqrt(mu/p);

v_pf_vect= co * [-sin(f); (e+cos(f)); 0];

% Compute position and velocity components in Geocentric RF
% to convert vectors from the PF ref frame to the geocentric one
% we multiply the vectors in the PF frame by the matrix T
r_vec = T*r_pf_vect;
v_vec = T*v_pf_vect;

r_vec = r_vec';
v_vec = v_vec';

end