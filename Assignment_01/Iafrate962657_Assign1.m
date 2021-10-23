% Spacecraft Guidance and Navigation (2021/2022)
% Assignment # 1
% Author: Davide Iafrate

%% Ex 1
clearvars; close all; clc

% Load kernels
cspice_furnsh('assignment01.tm');

% #1 PROPAGATOR VALIDATION
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

xend = flow2BP(x0,t0,T,primary.GM);

% If xend==x0 after one period our propagator is working correctly
norm(xend-x0);

% STM Validation
dr0 = [100 10 10]';
dv0 = [0.5 0 0]';

% Perturbed state vector
dx0 = [dr0; dv0];
t=0.1;

dxend = flow2BP(x0+dx0,t0,t,primary.GM) - xend;
STM = stateTransitionMatrix(t0,t,x0,primary.GM);
dxend1 = STM*dx0;

%% #2 LAMBERT PROBLEM SHOOTING SOLVER
% VALIDATE WITH EXAMPLE 8.8 of Curtis

% Primary parameters
primary.label = 'Sun';
primary.GM = cspice_bodvrd(primary.label,'GM',1);   % Gravitational param
                                                    
% select frame string (SPICE naming convention)
frame = 'ECLIPJ2000';

% Set departure/arrival data
departure.EpochStr = '1996-Nov-7 00:00:00.0000 TDB';
arrival.EpochStr = '1997-Sep-12 00:00:00.0000 TDB';
departure.Label = 'Earth';
arrival.Label = 'Mars';

departure.time = cspice_str2et(departure.EpochStr);
arrival.time = cspice_str2et(arrival.EpochStr);
ToF = arrival.time - departure.time;

% Retrieve initial and final positions of the planets
% Initial position vector [km]
r1 = cspice_spkpos(departure.Label,departure.time,frame,'NONE',primary.label);      
% Final position vector [km]
r2 = cspice_spkpos(arrival.Label,arrival.time,frame,'NONE',primary.label);      

% TEST
% r1=[1.0499e8;1.0465e8;716.93];
% r2 = [-2.0858e7;-2.1842e8;4.06244e6];

% Initial guess as Hohmann transfer velocity
E2M.objective = @(x) objectiveFun(x,departure.time,ToF,primary.GM,r2,r1);

% problem.x0 = [sqrt(primary.GM*(2/norm(r1) - 2/(norm(r1)+norm(r2))));0; 0];
E2M.x0 = [-24.429;21.782;0.9481];
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
% optsJacobian.CheckGradients = true;  % CHECK IF THE JACOBIAN IS CORRECT

E2M.options = optsJacobian;

sol2 = fsolve(E2M);

%% #3 delta_v minimizer
departure.Label = 'Earth';
arrival.Label = 'Mars';
departure.lb = cspice_str2et('2019-Jun-1 00:00:00.0000 TDB');
departure.ub = cspice_str2et('2021-Jun-1 00:00:00.0000 TDB');
arrival.lb = cspice_str2et('2020-Jan-1 00:00:00.0000 TDB');
arrival.ub = cspice_str2et('2022-Jan-1 00:00:00.0000 TDB');

% Define inital guess for the state (from computing DV grid) 
stateDep = cspice_spkezr(departure.Label,7511*24*3600,frame,'NONE',primary.label);
y0 = [stateDep; 7511*24*3600; 7716*24*3600];

% Define problem structure
E2M.objective = @(y) costFcn(y,departure,arrival,primary,frame); 
E2M.x0 = y0; % Initial guess for the state
E2M.lb = [-Inf*ones(6,1); departure.lb; arrival.lb];
E2M.ub = [+Inf*ones(6,1); departure.ub; arrival.ub];
E2M.nonlcon = @(y) enforcePositions(y,departure,arrival,primary,frame);
E2M.solver='fmincon';

% Define and set options
opts = optimoptions('fmincon');
opts.ConstraintTolerance = 1;
opts.Display = 'iter';
opts.PlotFcn = 'optimplotfval';
opts.MaxFunctionEvaluations = 1000;
E2M.options = opts;

% Compute solution
[sol,DeltaV,~,output,lambda,gradient,H] = fmincon(E2M);

cspice_kclear();

%% Ex 2
% clearvars; close all; clc

% Load kernels
% cspice_furnsh('assignment01.tm');

%% Functions

function xt = flow2BP(x0,t0,t,mu,posOutOnly)
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
opts = odeset('Reltol',1e-13,'AbsTol',1e-14,'Stats','off');

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
    g = flow2BP([r1; vi],t1,ToF,mu,1) - r2;
  
    % Use this to check if the jacobian is correct
%     https://www.mathworks.com/help/releases/R2021a/optim/ug/checking-validity-of-gradients-or-jacobians.html
    if nargout > 1
          x0 = [r1;vi];
          STM = stateTransitionMatrix(t1,t1+ToF,x0,mu);
          PHIrv = STM(1:3,4:6);
          
          J = PHIrv;
    
    end
end

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

function DV = costFcn(y,departureObj,arrivalObj,primary,frame)
%COSTFCN function that computes the total DELTA-V required by the
%trajectory

r1 = y(1:3);
v1 = y(4:6);
t1 = y(7);
t2 = y(8);

x2 = flow2BP([r1; v1],t1,t2,primary.GM,0);
v2 = x2(4:6);

stateDep = cspice_spkezr(departureObj.Label,t1,frame,'NONE',primary.label);
stateArrival = cspice_spkezr(arrivalObj.Label,t2,frame,'NONE',primary.label);

vDep = stateDep(4:6);
vArr = stateArrival(4:6);

dv1 = v1 - vDep;
dv2 = v2 - vArr;
DV = norm(dv1) + norm(dv2);
end

function [c, ceq] = enforcePositions(y,departureObj,arrivalObj,primary,frame)
% How can I avoid computing the flow AGAIN?
r1 = y(1:3);
v1 = y(4:6);
t1 = y(7);
t2 = y(8);

radiusDep = cspice_spkpos(departureObj.Label,t1,frame,'NONE',primary.label);
radiusArrival = cspice_spkpos(arrivalObj.Label,t2,frame,'NONE',primary.label);

x2 = flow2BP([r1; v1],t1,t2,primary.GM);
r2 = x2(1:3);

% Compute constraint violations
c1 = norm(r1-radiusDep);
c2 = norm(r2-radiusArrival);
ToF = t2-t1;

ceq = [c1; c2;];
c = ToF;
end

function STM = stateTransitionMatrix(t0,t,x0,mu)
%STM Summary of this function goes here
%   Detailed explanation goes here

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
    PHIvec = Y(7:42);
    PHI = reshape(PHIvec,[6 6]);
    dPHI = A*PHI;
    dPHIvec = dPHI(:);
    dy = [dX; dPHIvec];
end

STM0 = eye(6);
y0 = [x0;STM0(:)];
opts = odeset('Reltol',1e-13,'AbsTol',1e-14);

sol = ode113(@stmDynamics,[t0 t],y0,opts);
endState = sol.y(:,end);
STMvec = endState(7:42);
STM = reshape(STMvec,[6 6]);

end

% function dX = EOM(t,X,u,