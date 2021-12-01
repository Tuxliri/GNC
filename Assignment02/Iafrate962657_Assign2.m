% Spacecraft Guidance and Navigation (2021/2022)
% Assignment # 1
% Author: Davide Iafrate

%% Ex 1
clearvars; close all; clc

% Load kernels
cspice_furnsh('assignment02.tm');

% Get parameters needed
mu = cspice_bodvrd('Earth','GM',1);

station(1).name = 'Milano';
station(1).lat = 45.50122;          % station latitude [rad]
station(1).lon = 9.15461;          % station longitude [rad]
station(1).alt = 20;                        % station altitude [m]
station(1).frameName = 'MILANO_TOPO';       % Topocentric frame name [-]
station(1).minimumElevation = deg2rad(30);  % minimum elevation of the
                                            % tracked object [rad]
station(1).R = diag([deg2rad(100e-3) deg2rad(100e-3) 0.01e3]);  % Noise meaurement matrix
                                            % for Azimuth, Elevation, Range

station(2).name = 'Wellington';
station(2).lat = -41.28437;          % station latitude [rad]
station(2).lon = 174.76697;          % station longitude [rad]
station(2).alt = 117;                        % station altitude [m]
station(2).frameName = 'WELLINGTON_TOPO';       % Topocentric frame name [-]
station(2).minimumElevation = deg2rad(20);  % minimum elevation of the
                                            % tracked object [rad]
                              
station(3).name = 'La Silla';
station(3).lat = -29.26117;          % station latitude [rad]
station(3).lon = -70.73133;          % station longitude [rad]
station(3).alt = 2400;                        % station altitude [m]
station(3).frameName = 'LA-SILLA_TOPO';       % Topocentric frame name [-]
station(3).minimumElevation = deg2rad(10);  % minimum elevation of the
                                            % tracked object [rad]

% Initial et time
t0 = cspice_str2et('2021-Nov-19 14:25:39.6520 TDB');

% Final et time
tf = cspice_str2et('2021-Nov-23 14:00:00.0000 TDB');

% Define timespan vector
tspan = t0:60:tf;

% Define initial state vector and covariance
r0mean = [23721.610; 17903.673; -49.918];
v0mean = [-1.150987; 1.529718;  3.122389];

x0mean = [r0mean; v0mean];

P0 = [+2.6e-2 +1.4e-2 -1.8e-3 0 0 0;
      +1.4e-2 +1.8e-2 +2.3e-3 0 0 0;
      -1.8e-3 +2.3e-3 +1.0e-2 0 0 0;
       0 0 0 +1.6e-7 0 0;
       0 0 0 0 +1.6e-7 0;
       0 0 0 0 0 +1.6e-7];

% Set options for ODE solver
opts = odeset('reltol', 1e-12, 'abstol', 1e-12);

% Solve the ode
sol = ode113(@twobodyode,tspan,x0mean,opts,mu);

xx = deval(sol,tspan);
tt = tspan;
% Point 1

for i = 1:3
    
    % Latitude and longitude in degrees, altitude in meters
    lat = station(i).lat;
    lon = station(i).lon;
    alt = station(i).alt;
    
    minimumElevation = station(i).minimumElevation;
    VIS1(:,i) = visibilityOld(xx(1:3,:),tt,lat,lon,alt,minimumElevation);
    
%     VIS(:,i) = visibility(xx,tt,lat,lon,alt,minimumElevation);
    
    % Set first and last elements to 0
    VIS1(1,i) = 0; 
    VIS1(end,i) = 0;
    
    % Detect visibility windows
    riseIndexes = strfind(VIS1(:,i)',[0 1]);
    setIndexes = strfind(VIS1(:,i)',[1 0]);
    station(i).tRise = tt(riseIndexes);
    station(i).tSet = tt(setIndexes);
    station(i).visibilityMatrix = [station(i).tRise' station(i).tSet'];
    station(i).visibilityIndexes = [riseIndexes' setIndexes'];
end

% Point 2

tSet = sort([station(1).tSet station(2).tSet station(3).tSet]);    % assemble from each station

linearPt = zeros(6,6,length(tSet));

% Linear method
for i=1:length(tSet)
    tf = tSet(i);
    PHI = stateTransitionMatrix(t0,x0mean,tf,mu);
    linearPt(:,:,i) = PHI*P0*PHI';
end

% Unscented transform
n = 6;
unscentedPt = zeros(6,6,length(tSet));

% Generate sigma points
matrixRoot = sqrtm(n*P0);

chi0 = x0mean;

for i = 1:n
    chi(:,i) = x0mean + matrixRoot(i,:)';
%     chi(:,i) = x0 + matrixRoot(:,i);

end

for i = (n+1):2*n
    chi(:,i) = x0mean - matrixRoot(i-n,:)';
%     chi(:,i) = x0 - matrixRoot(:,i-n);

end

% Transform the sigma points using the nonlinear transform (flow of 2BP)
%%%% CHECK THE RESULTS, they seem to be too big %%%%
for i=1:length(tSet)
    tf = tSet(i);
    Y0 = flow2BP(t0,chi0,tf,mu);
    
    for j = 1:2*n
        x0mean = chi(:,j);
        Y(:,j) = flow2BP(t0,x0mean,tf,mu);
    end
    
    % Global matrix
    Ytot = [Y0 Y]';
    
    % Compute the mean
    Ymean = mean(Ytot);
    
    % Compute the covariance

    unscentedPt(:,:,i) = cov(Ytot);
end

% Point 3: Montecarlo analysis

% Draw 100 samples of the state from initial covariance P0
n = 100;    Sigma = P0;
R = mvnrnd(x0mean,Sigma,n);	% Each column of R contains a state to propagate

tf = tSet(end);

for i = 1:n
	xi = R(i,:);
	xt = flow2BP(t0,xi,tf,mu);
	B(i,:) = xt';
end

% Compute the mean of propagated states 
x_f_mean = mean(B);

% Compute the covariance of the propagated states
Pfinal = cov(B);

%% Exercise 2

% Replace the path below with your own installation path if not present as
% a subfolder
addpath('sgp4');

% Set SGP4 parameters
typerun    = 'u';  % user-provided inputs to SGP4 Matlab function
opsmode    = 'a';  % afspc approach ('air force space command')
whichconst =  72;  % WGS72 constants (radius, gravitational parameter)

% Improvement: READ FROM FILE
longstr1 = '1 28922U 05051A   21323.60115338 -.00000088  00000-0  00000+0 0  9998';
longstr2 = '2 28922  58.3917  37.3090 0004589  75.2965 284.7978  1.69421077 98469';

satrec = twoline2rv(longstr1, longstr2, typerun,'e',...
                            opsmode, whichconst);

% Get TLE epoch
[year,mon,day,hr,min,sec] = invjday(satrec.jdsatepoch, satrec.jdsatepochf);
sat_epoch_str = sprintf('%d-%02d-%02dT%02d:%02d:%02.6f', [year,mon,day,hr,min,sec]);
sat_epoch_et = cspice_str2et(sat_epoch_str);

fprintf('Satellite num ID: %d\n', satrec.satnum);
fprintf('TLE reference epoch: UTC %s', sat_epoch_str);

% Evaluate the TLE at its reference epoch
[satrec,rteme,vteme] = sgp4(satrec,  0.0);

% For each station
for k=1:1
    
    % Maybe it can be useful to extract also the time instants of the
    % observations
    stationMeasurements = measurementFun(satrec,station(k));
    station(k).measurements = stationMeasurements;
end

tf = cspice_str2et('2021-Nov-23 14:00:00.0000 TDB');
objective = @(x) costFcn(x,station(1),t0,tf,mu);
x0 = x0mean;
options = optimoptions('lsqnonlin');
options.Display = 'iter';

[x, resnorm, residual, exitflag, ~, ~, jac] =...
    lsqnonlin(objective,x0,[],[],options);

Jac =full(jac);
P_ls = resnorm/(length(residual)-length(x0)).*inv(Jac'*Jac);

% clean kernel pool
cspice_kclear();
%% Functions

function [azimuth, elevation, range, range_rate] = antenna_pointing(stationName, et, rv_target_eci)

%station_visibility Summary of this function goes here
%   Detailed explanation goes here

topoFrame = [stationName, '_TOPO'];

% Get from SPK kernel the station state in ECI reference frame
rv_station_eci = cspice_spkezr(stationName, et, 'J2000', 'NONE', 'EARTH');

% Compute state vector of satellite as seen from the station in J2000
rv_station_sat_eci = rv_target_eci - rv_station_eci;

% Get state rotation matrix from ECI (J2000) to TOPOCENTRIC
ROT_ECI2TOPO = cspice_sxform('J2000',topoFrame, et);

% Convert target ECI state into topocentric
rv_station_sat_topo = zeros(size(rv_station_sat_eci));
for i = 1:size(rv_station_sat_eci,2)
    rv_station_sat_topo(:,i) = ROT_ECI2TOPO(:,:,i)*rv_station_sat_eci(:,i);
end

% Compute range and range-rate
% Compute euclidean norm along direction 1 (column-wise)
range      = vecnorm(rv_station_sat_topo(1:3,:),2,1);

% Range rate is the scalar product of the velocity and a unit vector in the
% range direction
range_rate = dot(rv_station_sat_topo(4:6,:), rv_station_sat_topo(1:3,:)./range);

% Compute azimuth and elevation
rll_station_sat = cspice_xfmsta(rv_station_sat_topo,'RECTANGULAR','LATITUDINAL','EARTH');

%fprintf('%.5f %.5f\n',sat_range-rll_station_sat(1,:))
%fprintf('%.5f %.5f\n',sat_range_rate-rll_station_sat(4,:))

azimuth = rll_station_sat(2,:);
elevation = rll_station_sat(3,:);

end

function VIS = visibilityOld(r,et,lat,lon,alt,minimumElevation)
%VISIBILITY function that checks for visibility of the given object from a
%certain ground station
%
% INPUT:
%   r       position vector of the object in ECI coordinates
%   et      observation epoch
%   lat     station latitude
%   lon     station longitude
%   alt     station altitude
%
% OUTPUT:
%   VIS     logic visibility conidition
%               |-> 1 if visible
%               |-> 0 if not visible

% Get Earth radii (equatorial and polar)
% (requires pck00010.tpc kernel)

radii = cspice_bodvrd('EARTH', 'RADII', 3);
re = radii(1); rp = radii(3);

% Compute flattening
flat = (re - rp) / re;

% Convert to radians and km
lat_rad = deg2rad(lat); lon_rad = deg2rad(lon);
alt_km = alt / 1000;

% Compute station pos wrt Earth center (ECEF)
a = cspice_pgrrec(...
'EARTH', lon_rad, lat_rad, alt_km, re, flat);

% Compute ECEF2TOPO rotation matrix
ECEF2TOPO = cspice_eul2m(...
lat_rad - pi, pi - lon_rad, pi / 2, 2, 1, 2);

% Compute the east, north, up vectors in ECEF coordinates


for i = 1:length(et)
    
    % Convert object position from ECI to ECEF coordinates
    ROT_ECI2ECEF = cspice_pxform('J2000','IAU_EARTH',et(i));
    rECEF = ROT_ECI2ECEF*r(:,i);
    
    % Compute radius from station to satellite in ECEF
    rho_ECEF = rECEF-a;
    
    rho_TOPO = ECEF2TOPO*rho_ECEF;
    
    % compute radius in topocentric
    El(i) = asin(rho_TOPO(3)/norm(rho_TOPO));
end

VIS = (El > minimumElevation)';

% Plot something
end

function VIS = visibility(xECI,et,station)
%VISIBILITY function that checks for visibility of the given object from a
%certain ground station
%
% INPUT:
%   r       position vector of the object in ECI coordinates
%   et      observation epoch
%   lat     station latitude
%   lon     station longitude
%   alt     station altitude
%
% OUTPUT:
%   VIS     logic visibility conidition
%               |-> 1 if visible
%               |-> 0 if not visible

% Get Earth radii (equatorial and polar)
% (requires pck00010.tpc kernel)
stationName = station.name;

for i = 1:length(et)
    tt = et(i);
    rv_target = xECI(i,:);
    [Az, El(i), rho, drho] = antenna_pointing(stationName, tt, rv_target);
end

VIS = (El > station.minimumElevation)';

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

function STM = stateTransitionMatrix(t0,x0,t,mu)
%STATETRANSITIONMATRIX function computing the state transition matrix

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

function measurements = measurementFun(satrec,station)
    
% INPUT:
%
%   tt[1xN]         time vector of the propagated trajectory
%   xx[6xN]         state vector of the propagated trajectoy
%   station         struct with the station parameters


% Constant for arcseconds to radians conversions
arcsec2rad = pi / (180*3600);

visibilityMatrix = station.visibilityMatrix;

Sigma = station.R;

[rows,~] = size(visibilityMatrix);

tt = [];
yNoisy = [];

% For each visibility window
for window = 1:rows
  
    % JD dates of the visibility window
    tspan = station.visibilityMatrix(window,1):60:...
            station.visibilityMatrix(window,2);
    tt = [tt tspan];
    
    % Get reference epoch
%     jdsatepoch = satrec.jdsatepoch;
%     jdsatepochf = satrec.jdsatepochf;
%     referenceJD = jdsatepoch + jdsatepochf;
    
    % Get TLE epoch
    [year,mon,day,hr,min,sec] = invjday(satrec.jdsatepoch, satrec.jdsatepochf);
    sat_epoch_str = sprintf('%d-%02d-%02dT%02d:%02d:%02.6f', [year,mon,day,hr,min,sec]);
    epoch_et = cspice_str2et(sat_epoch_str);

    % corrective coefficients taken from EOP-Last5Years.txt
    ddpsi = -0.112684*arcsec2rad; %  [rad]
    ddeps = -0.006350*arcsec2rad; %  [rad]
    
    
    % for each time instant in the visibility windows of the station
    for jj = 1:length(tspan)
        
        % Propagate the states using SGP4
        minutesSinceRef = (tspan(jj)-epoch_et)/60;
        [satrec,rteme,vteme] = sgp4(satrec,  minutesSinceRef);
        
        et = tspan(jj);
        
        ttt = cspice_unitim(et, 'ET', 'TDT')/cspice_jyear()/100;
        
        % ---------------------- eci transformations ----------------------
        % -------- teme2eci    - transform teme to eci vectors
        ateme = [0;0;0];
        [reci, veci, ~] = teme2eci(rteme, vteme, ateme, ttt, ddpsi, ddeps);
        
        rv_target_eci = [reci; veci];
        
        stationName = station.name;
        % NB: for Milan the measurements are azimuth and elevation, for
        % the other 2 stations are right ascension and declination!!!
        % IMPROVEMENT: use array of function handles f(et,state)
        switch stationName
            case 'Milano'
                [azimuth, elevation, range, ~] =...
                    antenna_pointing(stationName, et, rv_target_eci);
                yMean(jj,1:3) = [rad2deg(azimuth) rad2deg(elevation) range];
            otherwise
                
        end
    end
    
    % Generate the noisy measurement and output the corresponding
    % measurement times
    measurements(window).yNoisy = mvnrnd(yMean,Sigma)';
    measurements(window).measurementTimes = tt;
end
end

function costVect = costFcn(initialStateGuess,station,t0,tf,mu)

tspan = t0:60:tf;

% Nx3 matrix of measurements
measurements = [station(1).measurements(1).yNoisy';
                station(1).measurements(2).yNoisy';
                station(1).measurements(3).yNoisy';
                station(1).measurements(4).yNoisy'];
            
% Column vector of measurement times
measurementTimes = station(1).measurements(1).measurementTimes';
   
% Set options for ODE solver
opts = odeset('reltol', 1e-12, 'abstol', 1e-12);

% Propagate trajectory on tspan
sol = ode113(@twobodyode,tspan,initialStateGuess,opts,mu);

xx = deval(sol,tspan);

% Compute the vector of measurements to be used
usedMeasurementsIndexes = measurementTimes<=tf;
usedMeasurements = measurements(usedMeasurementsIndexes,1:3);

stationName = station.name;

costVect = zeros(3*length(usedMeasurementsIndexes),1);
for i=1:length(usedMeasurementsIndexes)
    index = (measurementTimes(i)-tspan(1))/60 +1;
    
    rv_target_eci = xx(:,index);
    et = measurementTimes(i);
    
    [azimuth, elevation, range, ~] =...
                    antenna_pointing(stationName, et, rv_target_eci);
    estimatedMeasurement = [rad2deg(azimuth) rad2deg(elevation) range];
    
    residual = estimatedMeasurement'-usedMeasurements(i,1:3)';
    residual(3) = residual(3);
    costVect((3*i:(3*i+2))-2) = residual;
end
end