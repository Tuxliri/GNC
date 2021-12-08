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
station(1).R = diag([100e-3 100e-3 0.01e3].^2);  % Noise meaurement matrix
                                            % for Azimuth, Elevation, Range
station(1).FoV = 6;                  % Antenna field of view [deg]


station(2).name = 'Wellington';
station(2).lat = -41.28437;          % station latitude [rad]
station(2).lon = 174.76697;          % station longitude [rad]
station(2).alt = 117;                        % station altitude [m]
station(2).frameName = 'WELLINGTON_TOPO';       % Topocentric frame name [-]
station(2).minimumElevation = deg2rad(20);  % minimum elevation of the
                                            % tracked object [rad]
station(2).R = diag(rad2deg([0.0005 0.0005].^2));  % Noise meaurement matrix for
                                       % Right Ascension and Declination
station(2).FoV = 2;                  % Antenna field of view [deg]
                              
station(3).name = 'La-Silla';
station(3).lat = -29.26117;          % station latitude [rad]
station(3).lon = -70.73133;          % station longitude [rad]
station(3).alt = 2400;                        % station altitude [m]
station(3).frameName = 'LA-SILLA_TOPO';       % Topocentric frame name [-]
station(3).minimumElevation = deg2rad(10);  % minimum elevation of the
                                            % tracked object [rad]
station(3).R = diag(rad2deg([0.001 0.001].^2));  % Noise meaurement matrix
                                     % for Right Ascension and Declination
station(3).FoV = 1;                  % Antenna field of view [deg]

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

% -----Point 1-----

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

% -----Point 2-----

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

% -----Point 3: Montecarlo analysis-----

% Draw 100 samples of the state from initial covariance P0
n = 100;    Sigma = P0;
R = mvnrnd(x0mean,Sigma,n);	% Each column of R contains a state to propagate

tLastVisibility = tSet(end);
xxLastVisibility = deval(sol,tLastVisibility);

% LAST STATION IS LA-SILLA, find a way to automate this maybe

% The nominal center of the FoV of the station is in direction:
[azimuth, elevation, range, ~] = antenna_pointing(station(3).name,...
                                        tLastVisibility, xxLastVisibility);

antennaCenterVersor = [cos(azimuth)*cos(elevation) -sin(azimuth)*cos(elevation) sin(elevation)];

for i = 1:n
	xi = R(i,:);
	xt = flow2BP(t0,xi,tLastVisibility,mu);
    
    % Compute the satellite position versor
    [azimuth, elevation, range, ~] = antenna_pointing(station(3).name,...
                                            tLastVisibility, xt);
    satVersor(i,:) = [cos(azimuth)*cos(elevation) -sin(azimuth)*cos(elevation) sin(elevation)];
    
	B(i,:) = xt';
    
    % Compute how many points are inside the FoV of the station
end

% Compute the sample mean of propagated states 
x_f_mean = mean(B);

% Compute the sample covariance of the propagated states
Pfinal = cov(B);

% The nominal center of the FoV of the station is in direction:
[azimuth, elevation, range, ~] = antenna_pointing(station(3).name,...
                                        tLastVisibility, x_f_mean');

antennaCenterVersor = [cos(azimuth)*cos(elevation) -sin(azimuth)*cos(elevation) sin(elevation)];

for i = 1:n
    % Compute the angle B/W the satellite and the antenna center
    angle(i) = acos(dot(antennaCenterVersor,satVersor(i,:)));
end
    
% Plot the % of points inside the station FoV
insideFoVPoints = sum(rad2deg(angle)<(station(3).FoV)/2);
pie([insideFoVPoints 100-insideFoVPoints],{'Points inside FoV', 'Points outside FoV'})

sqrt(trace(unscentedPt(:,:,end)))
sqrt(trace(linearPt(:,:,end)))

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

% Generate the measurements for each station
for k=1:3
    
    % Maybe it can be useful to extract also the time instants of the
    % observations
    station(k).measurements = measurementFun(satrec,station(k));
    
end

tf = cspice_str2et('2021-Nov-23 14:00:00.0000 TDB');

x0=sol.y(:,1);
[residual0,parout]=costFcn(x0,station,t0,tf,mu,0);

objective = @(x) costFcn(x,station,t0,tf,mu,0);
x0 = x0mean;
options = optimoptions('lsqnonlin','Algorithm',...
                        'levenberg-marquardt','Display','iter');
                    
[x, resnorm, RES, exitflag, ~, ~, jac] =...
    lsqnonlin(objective,x0,[],[],options);

Jac =full(jac);
P_ls = resnorm/(length(RES)-length(x0)).*inv(Jac'*Jac);

% clean kernel pool
cspice_kclear();

%% Exercise 3

addpath(genpath('tdmReader'))
epoch = cspice_str2et('2016-10-10 00:00:00.000');

for i = 1:length(rdm)
    [xk(:,i),Pk(:,:,i)] = UKF(xk(:,i-1),Pk(:,:,i-1), yk, Rk,tk_1,tk,mu,stations);
end

%% Functions
function [sat_azimuth, sat_elevation, sat_range, sat_range_rate] =...
                        antenna_pointing(stationName, et, rv_target_eci)

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
sat_range      = vecnorm(rv_station_sat_topo(1:3,:),2,1);

% Range rate is the scalar product of the velocity and a unit vector in the
% range direction
sat_range_rate = dot(rv_station_sat_topo(4:6,:), rv_station_sat_topo(1:3,:)./sat_range);

% Compute azimuth and elevation
rll_station_sat = cspice_xfmsta(rv_station_sat_topo,'RECTANGULAR','LATITUDINAL','EARTH');

%fprintf('%.5f %.5f\n',sat_range-rll_station_sat(1,:))
%fprintf('%.5f %.5f\n',sat_range_rate-rll_station_sat(4,:))

sat_azimuth = rll_station_sat(2,:);
sat_elevation = rll_station_sat(3,:);

end

function varargout = measure(stationName, et, rv_target_eci)
%OUTPUT IN RADIANS!
% Get from SPK kernel the station state in ECI reference frame
rv_station_eci = cspice_spkezr(stationName, et, 'J2000', 'NONE', 'EARTH');
   
% Compute state vector of satellite as seen from the station in J2000
rv_station_sat_eci = rv_target_eci - rv_station_eci;

switch stationName
    case 'Milano'
        topoFrame = [stationName, '_TOPO'];

        % Get state rotation matrix from ECI (J2000) to TOPOCENTRIC
        ROT_ECI2TOPO = cspice_sxform('J2000',topoFrame, et);

        % Convert target ECI state into topocentric
        rv_station_sat_topo = zeros(size(rv_station_sat_eci));
        for i = 1:size(rv_station_sat_eci,2)
            rv_station_sat_topo(:,i) = ROT_ECI2TOPO(:,:,i)*rv_station_sat_eci(:,i);
        end

        % Compute range and range-rate
        % Compute euclidean norm along direction 1 (column-wise)
        range = vecnorm(rv_station_sat_topo(1:3,:),2,1);

        % Range rate is the scalar product of the velocity and a unit vector in the
        % range direction
        range_rate = dot(rv_station_sat_topo(4:6,:), rv_station_sat_topo(1:3,:)./range);

        % Compute azimuth and elevation
        rll_station_sat = cspice_xfmsta(rv_station_sat_topo,...
            'RECTANGULAR','LATITUDINAL','EARTH');

        azimuth = rll_station_sat(2,:);
        elevation = rll_station_sat(3,:);
        
        varargout = {azimuth, elevation, range, range_rate};
        
    case {'Wellington','La-Silla'}
        
        % Get position vector of the satellite as seen from the station in
        % ECI coordinates (J2000)
        rho = rv_station_sat_eci(1:3);
        
        % Compute the local right ascension and declination of the
        % satellite
        RA = atan2(rho(2),rho(1));
        DECL = asin(rho(3)/norm(rho));
        
        varargout = {RA, DECL};
        
    otherwise
        error('Station not part of the network!')
end


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
    [Az, El(i), rho, drho] = measure(stationName, tt, rv_target);
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

function dy = twobodyodeJ2(~,y,mu)
% Restructed two body problem ODE function
% 
% PROTOTYPE:
%   dy = twobodyode(t,y)
% 
% INPUT:
%   t[1]       Time                                             [T]
%   y[6x1]     State of the oscillator 
%              (position and velocity vectors)                  [L, L/T]
%   mu[1]      Planetary constant                               [L^3/T^2]
%   J2[1]      Second zonal harmonic                            [-]
%   R_e[1]     Equatorial radius of Earth                       [km]
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

J2 = 0.00108263;                % second zonal armonic of earth       [ - ]
R_e = cspice_bodvrd( 'EARTH', 'RADII', 3);
R_e = R_e(1);

% Calculate radius
r = norm(y(1:3));

% Set the derivatives of the state
dy = [ y(4);
       y(5);
       y(6);
        - mu*y(1)/r^3 + ((1.5*J2*mu*R_e^2)/r^4) * y(1) / r * (5*y(3)^2/r^2 - 1);
        - mu*y(2)/r^3 + ((1.5*J2*mu*R_e^2)/r^4) * y(2) / r * (5*y(3)^2/r^2 - 1);
        - mu*y(3)/r^3 + ((1.5*J2*mu*R_e^2)/r^4) * y(3) / r * (5*y(3)^2/r^2 - 3)];
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

% Extract station data
visibilityMatrix = station.visibilityMatrix;
stationName = station.name;
Sigma = station.R;

[rows,~] = size(visibilityMatrix);

% For each visibility window
for window = 1:rows
    
    yMean = [];
        
    % JD dates of the visibility window
    tStart = station.visibilityMatrix(window,1);
    tEnd = station.visibilityMatrix(window,2);
    
    tSpanWindow = tStart:60:tEnd;
    
    % Get TLE epoch
    [year,mon,day,hr,min,sec] = invjday(satrec.jdsatepoch, satrec.jdsatepochf);
    sat_epoch_str = sprintf('%d-%02d-%02dT%02d:%02d:%02.6f', [year,mon,day,hr,min,sec]);
    epoch_et = cspice_str2et(sat_epoch_str);

    % corrective coefficients taken from EOP-Last5Years.txt
    ddpsi = -0.112684*arcsec2rad; %  [rad]
    ddeps = -0.006350*arcsec2rad; %  [rad]
    
    
    % for each time instant in the visibility windows of the station
    for jj = 1:length(tSpanWindow)
        
        % Propagate the states using SGP4
        minutesSinceRef = (tSpanWindow(jj)-epoch_et)/60;
        [satrec,rteme,vteme] = sgp4(satrec,  minutesSinceRef);
        
        % Get the J2000 time instant
        et = tSpanWindow(jj);
        
        % Compute the time for the TEME2ECI conversion
        ttt = cspice_unitim(et, 'ET', 'TDT')/cspice_jyear()/100;
        
        % ---------------------- eci transformations ----------------------
        % -------- teme2eci    - transform teme to eci vectors
        ateme = [0;0;0];
        [reci, veci, ~] = teme2eci(rteme, vteme, ateme, ttt, ddpsi, ddeps);
        
        % Satellite 'real' state
        rv_target_eci = [reci; veci];
        
        % Generate reference measurements from reference state
        switch stationName
            case 'Milano'
                [azimuth, elevation, range, ~] = ...
                    measure(stationName, et, rv_target_eci);
                yMean(jj,1:3) = [rad2deg(azimuth) rad2deg(elevation) range];
            case {'Wellington','La-Silla'}
                [RA, DECL] = ...
                    measure(stationName, et, rv_target_eci);
                yMean(jj,1:2) = [rad2deg(RA) rad2deg(DECL)];
                
        end
    end
    
    % Generate the noisy measurement and output the corresponding
    % measurement times
    measurements(window).yNoisy = mvnrnd(yMean,Sigma)';
    measurements(window).measurementTimes = tSpanWindow;
end

end

function [residual,estMeasurements] = costFcn(initialStateGuess,stations,t0,tf,mu,J2Flag)

tspan = t0:60:tf;

% Set options for ODE solver
opts = odeset('reltol', 1e-12, 'abstol', 1e-12);

if J2Flag == 0
    % Propagate trajectory on tspan
    sol = ode113(@twobodyode,tspan,initialStateGuess,opts,mu);
else
    sol = ode113(@twobodyodeJ2,tspan,initialStateGuess,opts,mu);
end

xx = deval(sol,tspan);

residual = [];

for k=1:1%length(stations)
    T = struct2table(stations(k).measurements);
    
    % Nx3 or Nx2 matrix of measurements
    measurements = cell2mat(T{:,1}')';

    % Column vector of measurement times
    measurementTimes = cell2mat(T{:,2}');

    [res,parout] = difference(sol,stations(k),measurementTimes,measurements,tspan);
    estMeasurements = parout;
    residual = [residual; res];

end

end

function [weightedRes, parout] = difference(sol,station,measurementTimes,measurements,tspan)

    tf = tspan(end);

    % Compute the vector of measurements to be used
    usedMeasurementsIndexes = measurementTimes<=tf;
    usedMeasurements = measurements(usedMeasurementsIndexes,:);

    stationName = station.name;
    
    % Get the covariances of the measurements of the stations
    sigma_meas = diag(station.R);
    
    % Compute the weights
    W_m = diag(1 ./ sigma_meas);
    
%     costVect = zeros(3*length(usedMeasurementsIndexes),1);
    
    weightedRes = [];
    parout = [];
    
    for i=1:length(usedMeasurementsIndexes)
%         index = (measurementTimes(i)-tspan(1))/60 +1;

        rv_target_eci = deval(sol,measurementTimes(i));
        et = measurementTimes(i);
        
        switch stationName
            case 'Milano'
                [azimuth, elevation, range, ~] =...
                                measure(stationName, et, rv_target_eci);
                estimatedMeasurement = [rad2deg(azimuth) rad2deg(elevation) range];
                
            case {'Wellington','La-Silla'}
                [RA, DECL] =...
                        measure(stationName, et, rv_target_eci);
                estimatedMeasurement = [rad2deg(RA) rad2deg(DECL)];
        end
%         parout(:,i) = [parout estimatedMeasurement'];
        weightedRes = [weightedRes;
                       W_m*(estimatedMeasurement'-usedMeasurements(i,:)')];
        
    end
end

function [xk,Pk] = UKF(xk_1,Pk_1, yk, Rk,tk_1,tk,mu,stations)

n = length(xk_1);
l = length(yk);
%------ Compute the SIGMA points ------
matrixRoot = sqrtm(n*Pk_1);

chi_0_k_1 = xk_1;

chi_k_1 = zeros(n,2*n);

% TRY TO AVOID FOR LOOPS:
for i = 1:n
    chi_k_1(:,i) = xk_1 + matrixRoot(i,:)';
end

for i = (n+1):2*n
    chi_k_1(:,i) = xk_1 - matrixRoot(i-n,:)';
end

% Weights computation
alfa = 1e-3;
beta = 2;
k = 0;
lambda = alfa^2*(n+k)-n;

W0m = lambda/(n+lambda);
W0c = lambda/(n+lambda)+(1-alfa^2+beta);

Wi_mean = 1/(2*(n+lambda))*eye(2*n);
Wi_cov = Wi_mean;

W_mean = diag([W0m; diag(Wi_mean)]);
W_cov = diag([W0c; diag(Wi_cov)]);

%------ PREDICTION STEP ------
% Propagate sigma points to tk
chi0K = flow2BP(tk_1,chi_0_k_1,tk,mu);

chiK = zeros(n,2*n);

for i = 1:2*n
    x0 = chiK(:,i);
    chiK(:,i) = flow2BP(tk_1,x0,tk,mu);
end
    
% Global matrix
chiKtot = [chi0K chiK]';
    
% Compute the weighted mean
xk_apriori_mean = mean(repmat(diag(W_mean)',2*n+1,1).*chiKtot);
    
% Compute the wieghted covariance

Pk_apriori = weightedcov(chiKtot,W_mean,W_cov); % CHECK IF THE WEIGHTS ARE ADDED CORRECTLY

% Compute the sigma points of the predicted measurements
Y0k1 = measure(chi0K,tk,stations(4));
Y0k2 = measure(chi0K,tk,stations(5));

for i = 1:2*n
    x = chiK(:,i);
    stationA = measure(x0,tk,stations(4));
    stationB = measure(x0,tk,stations(5));
    Y_i_k(:,i) = [stationA; stationB];
end

Yk = [[Y0k1; Y0k2] Y_i_k];

Yk_mean_apriori = mean(Wi_mean.*Yk,2);

%------ UPDATE STEP ------
% Compute measurements covariance and cross-covariance
Pyy_k = cov(Wi_cov.*Yk) + Rk;   % NOT THE CORRECT WAY TO ADD WEIGHTS!
Pxy_k = xcov(chiK,Wi_cov.*Yk);

% Compute the gain
Kk = Pxy_k*inv(Pyy_k);

% Compute the a posteriori mean and covariance
xk = xk_apriori_mean + Kk*(yk - Yk_mean_apriori);
Pk = Pk_apriori - Kk*Pyy_k*(Kk');
end

function C = weightedcov(Y, Wm, Wc)
%   Weighted Covariance Matrix
%
%   WEIGHTEDCOV returns a symmetric matrix C of weighted covariances
%   calculated from an input T-by-N matrix Y whose rows are
%   observations and whose columns are variables and an input T-by-1 vector
%   w of weights for the observations. This function may be a valid
%   alternative to COV if observations are not all equally relevant
%   and need to be weighted according to some theoretical hypothesis or
%   knowledge.
%
%   C = WEIGHTEDCOV(Y, w) returns a positive semidefinite matrix C, i.e. all its
%   eigenvalues are non-negative.
%
%   If w = ones(size(Y, 1), 1), no difference exists between
%   WEIGHTEDCOV(Y, w) and COV(Y, 1).
%
%   REFERENCE: mathematical formulas in matrix notation are available in
%   F. Pozzi, T. Di Matteo, T. Aste,
%   "Exponential smoothing weighted correlations",
%   The European Physical Journal B, Volume 85, Issue 6, 2012.
%   DOI:10.1140/epjb/e2012-20697-x. 
%
% % ======================================================================
%
% % COMPUTE WEIGHTED COVARIANCE MATRIX
%   c = weightedcov(Y, w);                                                        % Weighted Covariance Matrix
%
% % ======================================================================
%
%   See also CORRCOEF, COV, STD, MEAN.
%   Check also WEIGHTEDCORRS (FE 20846) and KENDALLTAU (FE 27361)

% Check input
ctrl = isvector(w) & isreal(w) & ~any(isnan(w)) & ~any(isinf(w)) & all(w > 0);
if ~ctrl
  error('Check w: it needs be a vector of real positive numbers with no infinite or nan values!')
end
ctrl = isreal(Y) & ~any(isnan(Y)) & ~any(isinf(Y)) & (size(size(Y), 2) == 2);
if ~ctrl
  error('Check Y: it needs be a 2D matrix of real numbers with no infinite or nan values!')
end
ctrl = length(w) == size(Y, 2);
if ~ctrl
  error('size(Y, 1) has to be equal to length(w)!')
end

wm = diag(Wm);
wc = diag(Wc);

[T, N] = size(Y);                   % T: number of observations;
                                    % N: number of variables
C = Y - repmat(wm' * Y, T, 1);       % Remove mean (which is, also, weighted)
C = C' * (C .* repmat(w, 1, N));    % Weighted Covariance Matrix
C = 0.5 * (C + C');                 % Must be exactly symmetric

end