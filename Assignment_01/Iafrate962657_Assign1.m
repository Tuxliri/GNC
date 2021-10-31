% Spacecraft Guidance and Navigation (2021/2022)
% Assignment # 1
% Author: Davide Iafrate

%% Ex 1 - PRETTY OKAY
clearvars; close all; clc

% Load kernels
cspice_furnsh('assignment01.tm');

%%%%%%%%% #1 PROPAGATOR VALIDATION --- OK! %%%%%%%%%%
primary.label = 'Earth';
primary.GM = cspice_bodvrd(primary.label,'GM',1);   % Gravitational param
                                                    %   [km^3/s^2]
r0 = [ -7128.137, 0, 0 ]';   % Initial radius vector             [ km ]
v0 = [ 0, -9.781, 0 ]';      % Initial velocity vector           [ km/s ]

% Specific energy (constant)                                    [ km^2/s^2 ]
epsilon0 = 0.5*norm(v0)^2 - primary.GM/norm(r0);

a = -primary.GM/(2*epsilon0);        % Semi-major axis                   [ km ]
T = 2*pi*sqrt(a^3/primary.GM);       % Period                            [ s ]

% Initial conditions
x0 = [r0;v0];
t0 = 0;

xend = flow2BP(t0,x0,T,primary.GM);

% If error = 0 after one period our propagator is working correctly
error = norm(xend-x0);

%%%%%%%%% LAMBERT PROBLEM SHOOTING SOLVER - Validated %%%%%%%%%
%(Earth centered solution) - Validated
% Primary parameters 
primary.label = 'Earth';
primary.GM = cspice_bodvrd(primary.label,'GM',1);   % Gravitational param
                                                    
% select frame string (SPICE naming convention)
frame = 'J2000';

% select departure and arrival strings
departure.EpochStr = '1996-Nov-7 00:00:00.0000 TDB';
arrival.EpochStr = '1997-Sep-12 00:00:00.0000 TDB';

% Convert strings to J2000 times
departure.time = cspice_str2et(departure.EpochStr);
arrival.time = cspice_str2et(arrival.EpochStr);
ToF = arrival.time - departure.time;

% Retrieve initial and final positions of the planets
% Initial position vector [km]
r1 = [ 7000; 0; 1];
v1_0 = [0; 8.5; 0];
t1 = 0;
ToF = 2000;

% Final position vector [km]
r2 = [200; 7000; -100];

% Initial guess as Hohmann transfer velocity
E2M.objective = @(x) objectiveFun(x,0,ToF,primary.GM,r2,r1);
E2M.x0 = v1_0;
E2M.solver = 'fsolve';

% Set options for using finite difference method
optsFinDiff = optimoptions('fsolve');
optsFinDiff.Display='iter';
optsFinDiff.FiniteDifferenceType = 'central';

E2M.options = optsFinDiff;

sol=fsolve(E2M);

% Set options for using the Jacobian computed through the STM
optsJacobian = optimoptions('fsolve');
optsJacobian.Display='iter';
optsJacobian.SpecifyObjectiveGradient = true;

E2M.options = optsJacobian;

sol2 = fsolve(E2M);

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

sol=fsolve(E2M);

[A,P,E,ERROR,VI,VF,TPAR,THETA] = lambertMR(r1, r2, ToF, primary.GM,0,0,0 );

disp('The difference between the analytical lambert solver and the shooting one is:')
norm(sol - VI)

%%%%%%%%% delta_v minimizer %%%%%%%%%
clear departure arrival
departure.Label = 'Earth';
arrival.Label = 'Mars';


% Time conversion constant
DAY2SECS = 24*3600;

departure.lb = cspice_str2et('2022-Jan-1 00:00:00.0000 TDB')/DAY2SECS;
departure.ub = cspice_str2et('2023-Jan-1 00:00:00.0000 TDB')/DAY2SECS;

arrival.lb = cspice_str2et('2022-Mar-1 00:00:00.0000 TDB')/DAY2SECS;
arrival.ub = cspice_str2et('2023-Dec-31 00:00:00.0000 TDB')/DAY2SECS;


% Define inital guess for the state (from computing DV grid) 
departureguess = cspice_str2et('2022-Sep-9 00:00:00.0000 TDB');

minToF = 300; % minimum ToF in days
arrivalguess = departureguess/DAY2SECS + minToF;

stateDep = cspice_spkezr(departure.Label,departureguess,frame,'NONE',primary.label);
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
            38; 38; 10;
            departure.ub;
            arrival.ub];
        
E2MT.nonlcon = @(y) enforcePositions(y,departure,arrival,primary,frame);
E2MT.solver='fmincon';

% Define and set options
opts = optimoptions('fmincon');
opts.ConstraintTolerance = 1;
opts.Display = 'iter';
opts.PlotFcn = 'optimplotfval';
opts.ScaleProblem = true;

E2MT.options = opts;

% Compute solution (why doesn't it output the minimum DeltaV found??)
[sol,DeltaV,~,output,lambda,gradient,H] = fmincon(E2MT);

% Checks
% End state position error
departureTime = sol(7)*DAY2SECS;
arrivalTime = sol(8)*DAY2SECS;
departureState = sol(1:6);
positionError = flow2BP(departureTime,departureState,arrivalTime,primary.GM,1)...
    - cspice_spkpos(arrival.Label,arrivalTime,frame,'NONE',primary.label)
departureDate = cspice_et2utc(departureTime,'C',0)
arrivalDate = cspice_et2utc(arrivalTime,'C',0)

% Clear kernel pool
cspice_kclear();

%% Ex 2
clearvars; close all; clc

% Load kernels
cspice_furnsh('assignment01.tm');

% Define list of celestial bodies:
labels = {'Sun';
          'Earth';
          'Moon'};

% Initialize propagation data
bodies = nbody_init(labels);

% select integration frame string (SPICE naming convention)
frame = 'J2000';


figure()
hold on
plot3(xx(:,1),xx(:,2),xx(:,3))
hold on
plot3(xx1(:,1),xx1(:,2),xx1(:,3))
axis equal

% Create problem structure
test = totalDV(y0,bodies);
L2.objective = @(y) totalDV(y,bodies);
L2.x0 = y0;
L2.nonlcon = @(y) constraints(y,bodies);
L2.solver = 'fmincon';

% Define and set options
opts = optimoptions('fmincon');
opts.ConstraintTolerance = 1;
opts.Display = 'iter';
opts.PlotFcn = 'optimplotfval';
opts.Algorithm = 'sqp';
L2.options = opts;

[sol,DeltaV,~,output,lambda,gradient,H] = fmincon(L2);

cspice_kclear();
%% Ex 3 - PRETTY OKAY
clearvars; close all; clc

% Load kernels
cspice_furnsh('assignment01.tm');

% Define boundary conditions
r0 = [0 -29597.43 0]';
v0 = [1.8349 0.0002 3.1783]';
m0 = 735;
rf = [0 -29617.43 0]';
vf = [1.8371 0.0002 3.1755]';

% Define S/C parameters
Tmax = 300*1e-6;
Isp = 3000;
t0 = 0;
GM = 398600;
epsilon0 = 0.5*norm(v0)^2 - GM/norm(r0);
semiMajorAxis = - GM/(2*epsilon0);
T = 2*pi*sqrt(semiMajorAxis^3/GM);

DAY2SECS = 24*60*60;
% Initial guess for the unknown
% Lambda0_guess = [.1 .1 .1 .1 .1 .1 .1 1e4]'; %(Wrong time)?
Lambda0_guess = [ 5 5 5 1e4 1e4 1e4 1 T/DAY2SECS]';  % By professor
options = optimoptions('fsolve');
options.Display = 'iter';
options.MaxFunctionEvaluations = 3000;
options.MaxIterations = 500;

x0 = [r0; v0; m0];

[Z,fval] = fsolve(@(x) shootingFcn(x,x0,t0,rf,vf,GM,Tmax,Isp),Lambda0_guess,options);
[nope, SOL] = shootingFcn(Z,x0,t0,rf,vf,GM,Tmax,Isp);

% Plot powered orbit of S/C transfer
plot3(SOL.y(1,:),SOL.y(2,:),SOL.y(3,:))
axis equal
grid on

% Plot lambdaM (should always be >0 for time-optimal problems)
figure()
plot(SOL.x,SOL.y(14,:),'LineWidth',2)
title('$\lambda_m$','Interpreter','latex')
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
title('Control magnitude')
grid on

% Compute the solution endpoint defect
norm(SOL.y(1:6,end)- [rf;vf],inf)

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

function xt = flow2BP(t0,x0,t,mu,posOutOnly)
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

end

function [g,J] = objectiveFun(vi,t1,ToF,mu,r2,r1)
    g = flow2BP(t1,[r1; vi],t1+ToF,mu,1) - r2;
  
    % Use this to check if the jacobian is correct
%     https://www.mathworks.com/help/releases/R2021a/optim/ug/checking-validity-of-gradients-or-jacobians.html
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

stateDep = cspice_spkezr(departureObj.Label,t1,frame,'NONE',primary.label);
stateArrival = cspice_spkezr(arrivalObj.Label,t2,frame,'NONE',primary.label);

vDep = stateDep(4:6);
vArr = stateArrival(4:6);

dv1 = norm(v1 - vDep);
dv2 = norm(v2 - vArr);
DV = dv1 + dv2;
end

function [c, ceq] = enforcePositions(y,departureObj,arrivalObj,primary,frame)
% How can I avoid computing the flow AGAIN?

% Time conversion constant
DAY2SECS = 24*3600;

r1 = y(1:3);
v1 = y(4:6);

t1 = y(7)*DAY2SECS;
t2 = y(8)*DAY2SECS;

radiusDep = cspice_spkpos(departureObj.Label,t1,frame,'NONE',primary.label);
radiusArrival = cspice_spkpos(arrivalObj.Label,t2,frame,'NONE',primary.label);

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

eqConstr = [errorPosition; errorVelocity; LambdaM*1000; errorH];

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
%   The integration centre is the Solar-System-Barycentre (SSB) and only
%   Newtonian gravitational accelerations are considered.
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

if not( strcmpi(frame, 'ECLIPJ2000') || strcmpi(frame, 'J2000') )
    msg = 'Invalid integration reference frame, select either J2000 or ECLIPJ2000';
    error(msg);
end

% Initialize right-hand-side
dxdt = zeros(6,1);

% Position detivative is object's velocity
dxdt(1:3) = x(4:6);
% Extract the object position from state x
rr_es_obj = x(1:3);

for i=1:length(bodies)

    % Retrieve position and velocity of i-th celestial body wrt Earth
    % Barycentre in an inertial frame
    rv_es_body = cspice_spkezr(bodies{i}.name, t, frame, 'NONE', 'EARTH');
    
    % Extract object position wrt. i-th celestial body
    dd_sc_body = rr_es_obj - rv_es_body(1:3);
    
    % Extract i-th body pos wrt central
    rr_body_prim = rv_es_body(1:3);
    rho2 = max(dot(rr_body_prim,rr_body_prim),0.01);
    rho = sqrt(rho2);
    
    % Compute square distance and distance
    dist2 = dot(dd_sc_body, dd_sc_body);
    dist = sqrt(dist2);

    % Compute the gravitational acceleration using Newton's law
    aa_grav =  - bodies{i}.GM * (dd_sc_body /(dist*dist2) + rr_body_prim/(rho*rho2));

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

options = odeset('reltol', 1e-12, 'abstol', 1e-12);
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
muE = bodies{2}.GM;
radiiE = cspice_bodvrd( 'EARTH', 'RADII', 3);
R_e = radiiE(1);

% Second burn time
t3 = y(25);

%First and last burn times
t1 = y(26);     tN = y(27);

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
       norm(v1) - sqrt(muE/norm(r1));
       dot(v1,r1);
       flowNBP(t(1),x1,t(2),bodies,frame) - x2;
       flowNBP(t(2),x2,t3,bodies,frame,1) - r3;
       flowNBP(t3,x3,t(4),bodies,frame) - x4;
       flowNBP(t(4),x4,tN,bodies,frame,1) - rL2_f];
   
c = [5e5 - norm(r3);
    t3 - tN;
    t1 - t3;
    t3 - tN];

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
muE = bodies{2}.GM;

% Second burn time
t3 = y(2);

%First and last burn times
t1 = y(26);     tN = y(27);

% Perform time discretization of each arc in 2 segments
N = 3;
t = [linspace(t1,t3,N) linspace(t3,tN,N)];

% Retrieve final orbital position of moon from ephemeris
xMoon_f = cspice_spkezr('moon',tN,'J2000','NONE','EARTH');

% Find L2 final state, same velocity as the moon but radius 60000km greater
vL2_f = xMoon_f(4:6);

% Compute the flows
x3f = flowNBP(t(2),x2,t3,bodies,frame);
x5f = flowNBP(t(4),x4,tN,bodies,frame);

% Compute the velocities before and after each impulse
v0 = cross([0 0 1],r1/norm(r1))*sqrt(muE/norm(r1));

% Compute the DVs
DV1 = norm(v1 - v0);
DV2 = norm(v3 - x3f(4:6));
DV3 = norm(vL2_f - x5f(4:6));

DV = DV1 + DV2 + DV3;
end

function y0 = initialGuess()
% Initial guess definition (use 2BP bielliptic transfer as reference)
opts = odeset('reltol', 1e-12, 'abstol', 1e-12);

GM = 398600;
a1 = (200 + 500000 + 6378)/2;
a2 = (500e3 + (378+60)*1e3 )/2;

T1 = 2*pi*sqrt(a1^3/GM);
T2 = 2*pi*sqrt(a2^3/GM);

ri = [1;0;0]*(6378+200);
vi = sqrt(2*GM/norm(ri) - GM/a1)*[0;1;0];
t0 = 0;

xi = [ri; vi];
[tt,xx] = ode113(@twobodyode,[t0 T1/2],xi,opts,GM);
xf = xx(end,:)';
v2 = sqrt(2*GM/norm(xf(1:3)) - GM/a2)*[0; -1; 0];
xi2 = [xf(1:3); v2];
[tt1,xx1] = ode113(@twobodyode,[t0 T2/2],xi2,opts,GM);

Npoints = length(tt);
Npoints1 = length(tt1);

y0 = [xi; xx(round(Npoints/2),:)'; xx1(1,:)'; xx(round(Npoints1/2),:)'; t0; t0+T1/2; t0+T1/2+T2/2 ];
