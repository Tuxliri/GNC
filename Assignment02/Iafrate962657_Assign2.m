% Spacecraft Guidance and Navigation (2021/2022)
% Assignment # 1
% Author: Davide Iafrate

%% Exercise 1
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

n = length(x0mean);

% Set options for ODE solver
opts = odeset('reltol', 1e-12, 'abstol', 1e-12);

% Solve the ode
sol = ode113(@twobodyode,tspan,x0mean,opts,mu);

xx = deval(sol,tspan);
tt = tspan;

% -----Point 1-----
figure()
hold on
for i = 1:3
    
    % Latitude and longitude in degrees, altitude in meters
    lat = station(i).lat;
    lon = station(i).lon;
    alt = station(i).alt;
    
    legendVisibility{i} = station(i).name;
    
    minimumElevation = station(i).minimumElevation;
    [VIS1(:,i),El(:,i)] = visibilityOld(xx(1:3,:),tt,lat,lon,alt,minimumElevation);
    
%     VIS(:,i) = visibility(xx,tt,lat,lon,alt,minimumElevation);
    
    % Set first and last elements to 0
    VIS1(1,i) = 0; 
    VIS1(end,i) = 0;
    Elevations = El(:,i).*(VIS1(:,i)>0);
    Elevations(Elevations==0) = NaN;
    
    plot(tt,rad2deg(Elevations),'LineWidth',1.5)
    
    % Detect visibility windows
    riseIndexes = strfind(VIS1(:,i)',[0 1]);
    setIndexes = strfind(VIS1(:,i)',[1 0]);
    station(i).tRise = tt(riseIndexes);
    station(i).tSet = tt(setIndexes);
    station(i).visibilityMatrix = [station(i).tRise' station(i).tSet'];
    station(i).visibilityIndexes = [riseIndexes' setIndexes'];
end

% Plot minimum elevation lines
yline(30,'--','30??')
yline(20,'--','20??')
yline(10,'--','10??')

grid on
grid minor
legend(legendVisibility{:},'','','','Location','best')
ylabel('Elevation [$\deg$]','Interpreter','latex')
xLabels = cspice_et2utc(xticks(gca),'C',0);
xLabels = xLabels(:,1:17);
xticklabels(xLabels)
ylim([0 90])

% -----Point 2-----

tSet = sort([station(1).tSet station(2).tSet station(3).tSet]);    % assemble from each station

linearPt = zeros(6,6,length(tSet));

% Linear method (LinCOV)
for i=1:length(tSet)
    tf = tSet(i);
    PHI = stateTransitionMatrix(t0,x0mean,tf,mu);
    linearPt(:,:,i) = PHI*P0*PHI';
    sqrtLinCOVPrr(i) = sqrt(trace(linearPt(1:3,1:3,i)));
    sqrtLinCOVPvv(i) = sqrt(trace(linearPt(4:6,4:6,i)));

end

% Unscented transform method

for i = 1:length(tSet)
    tf = tSet(i);
    fun = @(x0) flow2BP(t0,x0,tf,mu);
    [xf_mean,Pt_i] = unscentedTransform(x0mean,P0,fun);
    unscentedPt(:,:,i) = Pt_i;
    sqrtUTPrr(i) = sqrt(trace(Pt_i(1:3,1:3)));
    sqrtUTPvv(i) = sqrt(trace(Pt_i(4:6,4:6)));
end

figure()
hold on
scatter([t0 tSet],[sqrt(trace(P0(1:3,1:3))) sqrtLinCOVPrr],'filled')
scatter(tSet, sqrtUTPrr,'d')
grid on
ylabel('$\sqrt{tr(P_{rr})}$ [km]','Interpreter','latex')
grid minor
xLabels = cspice_et2utc(xticks(gca),'C',0);
xLabels = xLabels(:,1:17);
xticklabels(xLabels)
legend('linCOV','UT method')

figure()
hold on
scatter([t0 tSet],[sqrt(trace(P0(4:6,4:6))) sqrtLinCOVPvv],'filled')
scatter(tSet, sqrtUTPvv,'d')
grid on
ylabel('$\sqrt{tr(P_{vv})}$ [km/s]','Interpreter','latex')
grid minor
xLabels = cspice_et2utc(xticks(gca),'C',0);
xLabels = xLabels(:,1:17);
xticklabels(xLabels)
legend('linCOV','UT method')

% -----Plot the uncertainy ellipsoid -----
confidenceLevel = 95/100;   % Confidence level [%]

% Plot the linear Transform uncertainty ellipsoid
figure()
error_ellipse(linearPt(1:3,1:3,end),[0;0;0],'conf',confidenceLevel)

axis equal
grid on
grid minor

xlabel('x [km]')
ylabel('y [km]')
zlabel('z [km]')
view(-31,20)
alpha(.8)

% Plot the unscented Transform uncertainty ellipsoid
figure()
error_ellipse(unscentedPt(1:3,1:3,end),[0;0;0],'conf',confidenceLevel)
axis equal
grid on
grid minor

xlabel('x [km]')
ylabel('y [km]')
zlabel('z [km]')
view(-31,20)
alpha(.9)



% -----Point 3: Montecarlo analysis-----

% Draw 300 samples of the state from initial covariance P0
n = 300;    Sigma = P0;
R = mvnrnd(x0mean,Sigma,n);	% Each column of R contains a state to propagate

tLastVisibility = tSet(end);
xxLastVisibility = deval(sol,tLastVisibility);


% LAST STATION IS LA-SILLA

% The nominal center of the FoV of the station is in direction:
[azimuth, elevation, ~, ~] = antenna_pointing(station(3).name,...
                                        tLastVisibility, xxLastVisibility);

% antenna center pointing towards the propagated mean state
antennaCenterVersor = [cos(azimuth)*cos(elevation) -sin(azimuth)*cos(elevation) sin(elevation)];

for i = 1:n
	xi = R(i,:);
	xt = flow2BP(t0,xi,tLastVisibility,mu);
    
    % Compute the satellite position versor
    [azimuth, elevation, ~, ~] = antenna_pointing(station(3).name,...
                                            tLastVisibility, xt);
    satVersor(i,:) = [cos(azimuth)*cos(elevation) -sin(azimuth)*cos(elevation) sin(elevation)];
    
	B(i,:) = xt';
    
    % Compute how many points are inside the FoV of the station
end

% Compute the sample mean of propagated states 
x_f_mean = mean(B);

% Compute the sample covariance of the propagated states
Pfinal = cov(B);

% Plot the montecarlo position distribution
figure()
hold on
axis equal
scatter3(B(:,1),B(:,2),B(:,3),'filled')
error_ellipse(Pfinal(1:3,1:3),x_f_mean(1:3)','conf',confidenceLevel)
alpha(0.5)
grid on
grid minor

xlabel('x [m]')
ylabel('y [m]')
zlabel('z [m]')

% The nominal center of the FoV of the station is in direction:
[azimuth, elevation, ~, ~] = antenna_pointing(station(3).name,...
                                        tLastVisibility, x_f_mean');

% antenna center pointing towards the mean of montecarlo realizationa
antennaCenterVersor = [cos(azimuth)*cos(elevation) -sin(azimuth)*cos(elevation) sin(elevation)];

for i = 1:n
    % Compute the angle B/W the satellite and the antenna center
    angle(i) = acos(dot(antennaCenterVersor,satVersor(i,:)));
end

figure()
% Plot the % of points inside the station FoV
insideFoVPoints = sum(rad2deg(angle)<(station(3).FoV)/2);
outsideFoVPoints = n-insideFoVPoints;
labels = {strcat(num2str(insideFoVPoints/n*100,3),'% Points inside FoV'),...
          strcat(num2str(outsideFoVPoints/n*100,3),'% Points outside FoV')};
pie([insideFoVPoints outsideFoVPoints],labels)

%% Exercise 2
close all

% Replace the path below with your own installation path if not present as
% a subfolder
addpath('sgp4');

mu = cspice_bodvrd('Earth','GM',1);

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
    station(k).measurements = measurementFun(satrec,station(k));
end

tfList = [cspice_str2et('2021-Nov-21 14:00:00.0000 TDB');
      cspice_str2et('2021-Nov-22 14:00:00.0000 TDB');
      cspice_str2et('2021-Nov-23 14:00:00.0000 TDB')];

P_ls = zeros(6,6,3);
  
USEJ2ORBIT = 1;
x0 = x0mean;
    
  
options = optimoptions('lsqnonlin','Algorithm',...
        'levenberg-marquardt','Display','iter','StepTolerance',1e-10);
    
for jj = 1:3
    tf = tfList(jj);
    objective = @(x) costFcn(x,station,sat_epoch_et,tf,mu,USEJ2ORBIT);

    [x(:,jj), resnorm, RES, ~, ~, ~, jac] =...
        lsqnonlin(objective,x0,[],[],options);
    
    Jac = full(jac);
    P_ls(:,:,jj) = resnorm/(length(RES)-length(x0)).*inv(Jac'*Jac);
end

faceNumber = [20 20 20];
spacing = (0.1:0.1:0.3)/2;
colours = {'#EDB120','#D95319','#4DBEEE'};
legendEntries = {'2021-Nov-21 14:00', '2021-Nov-22 14:00', '2021-Nov-23 14:00'};

% Plot the three final position uncertainty ellipsoids
for jj = 1:3
    figure()
    mean_center = real(x(1:3,jj)');
    error_ellipse(real(P_ls(1:3,1:3,jj)),[0,0,0])
    axis equal
    xt = xticks;
    yt = yticks;
    zt = zticks;
    xticklabels(num2cell(xt*1e3))
    yticklabels(num2cell(yt*1e3))
    zticklabels(num2cell(zt*1e3))
    grid on
    grid minor
    
    xlabel('x [m]')
    ylabel('y [m]')
    zlabel('z [m]')
    legend(legendEntries{jj},'Location','best')
end

% ----- Point (2b) ---- With a priori information
tf = tfList(1);

objective = @(x) costFcn(x,station,sat_epoch_et,tf,mu,1,P0,x0mean);
[x_apriori, resnorm, RES, exitflag, ~, ~, jac] =...
    lsqnonlin(objective,x0,[],[],options);

Jac_apriori =full(jac);
P_ls_apriori = resnorm/(length(RES)-length(x0)).*inv(Jac_apriori'*Jac_apriori);

% clean kernel pool
cspice_kclear();

%% Exercise 3 

% Load kernels
cspice_furnsh('assignment02.tm');
cspice_furnsh('kernels\exomars.bsp');

% Get parameters needed
mu = cspice_bodvrd('Sun','GM',1);

% Read all the measurements
fileStruct    = dir(fullfile('tdm','*.tdm'));
fileNames     = fullfile({fileStruct.folder}, {fileStruct.name});

nFiles       = length(fileNames);
successFlags = true(1,nFiles);
meCA         = cell(1,nFiles);

measurementsVect = [];
measurementTimes = [];
STATION_ID = [];

for iFile = 1:nFiles
   try
       % READ THE DATA HERE
       fileName = fileNames{iFile};
       [data,et,metaData] = readTDM(fileName);
       switch metaData.PARTICIPANT_1
           
           case 'MALARGUE'
               STATION_ID = [STATION_ID;
                             4*ones(size(et))'];
           case 'NEW_NORCIA'
               STATION_ID = [STATION_ID;
                             5*ones(size(et))'];

       end
       
       measurementsVect = [measurementsVect; data];
       measurementTimes = [measurementTimes; et'];
   catch ME
       meCA{iFile} = ME;
       successFlags(iFile) = false;
   end
end

% Sort chronologically the measurements and measurement times

[measurementTimes,order] = sort(measurementTimes);
measurementsVect = measurementsVect(order,:);
STATION_ID = STATION_ID(order);

% Build the Rk matrix
sigma_Az = deg2rad(1.5e-3);
sigma_El = deg2rad(1.3e-3);
sigma_rng = 75e-3;
rho_Az_El = 0.1;
rho_Az_rng = 0;
rho_El_rng = 0;

R = [sigma_Az^2                     rho_Az_El*sigma_Az*sigma_El 0;
     rho_Az_El*sigma_Az*sigma_El    sigma_El^2                  0;
     0                                      0         sigma_rng^2];
 
station(4).name = 'MALARGUE';
station(4).lat = -35.77601;                 % station latitude [rad]
station(4).lon = -69.39819;                 % station longitude [rad]
station(4).alt = 1550;                      % station altitude [m]
station(4).frameName = 'Malargue_TOPO';     % Topocentric frame name [-]
station(4).minimumElevation = deg2rad(20);  % minimum elevation of the
                                            % tracked object [rad]


station(4).R = R;  % Noise meaurement matrix

station(5).name = 'NEW_NORCIA';
station(5).lat = -31.04823;                 % station latitude [rad]
station(5).lon = 116.19147;                 % station longitude [rad]
station(5).alt = 252;                       % station altitude [m]
station(5).frameName = 'NEW_NORCIA_TOPO';   % Topocentric frame name [-]
station(5).minimumElevation = deg2rad(20);  % minimum elevation of the
                                            % tracked object [rad]


station(5).R = R;                           % Noise meaurement matrix

addpath(genpath('tdm'))
epoch = cspice_str2et('2016-10-10T00:00:00.000');

r0mean = [+1.68467660241E+08 -1.07050886902E+08 -5.47243873455E+07];
v0mean = [+1.34362486580E+01 +1.68723391839E+01 +8.66147058235E+00];

x0mean = [r0mean'; v0mean'];

P0 = [+2.01e+04  -7.90e+00  -4.05e+00  -5.39e-03  +6.37e-06  +3.26e-06 ;
  -7.90e+00  +2.01e+04  +2.64e+00  +6.37e-06  -5.38e-03  -2.07e-06 ;
  -4.05e+00  +2.64e+00  +2.01e+04  +3.25e-06  -2.03e-06  -5.38e-03 ;
  -5.39e-03  +6.37e-06  +3.25e-06  +1.92e-07  -2.28e-09  -1.16e-09 ;
  +6.37e-06  -5.38e-03  -2.03e-06  -2.28e-09  +1.91e-07  +7.31e-10 ;
  +3.26e-06  -2.07e-06  -5.38e-03  -1.16e-09  +7.31e-10  +1.91e-07];

% Create a structure for all the measurements with fields:
%   - Measurements
%   - Measurement times
%   - Station ID

target = '-143';
frame = 'J2000';
observer = 'SUN';

TIMEVECTOR = [epoch; measurementTimes];
MEASUREMENTS_MATRIX = [0 0 0; measurementsVect]';
STATION_ID = [0; STATION_ID];
Pk(:,:,1) = P0;
xk(:,1) = x0mean;
rv_exomars0 = cspice_spkezr(target,epoch,frame,'NONE',observer);
errorTrueTraj(:,1) = xk(:,1) - rv_exomars0;

Rk = R;
tk_1 = TIMEVECTOR(1);

errorTrueTraj = zeros(6,length(TIMEVECTOR)-1);
uncertaintyCovariance = errorTrueTraj;

% [Az,El,range] = measure(station(STATION_ID(1)).name, TIMEVECTOR(1), rv_target_eci);

for i = 2:length(TIMEVECTOR)-1
    tk = TIMEVECTOR(i);
    yk = MEASUREMENTS_MATRIX(:,i);
    stat = station(STATION_ID(i));
    [xk(:,i),Pk(:,:,i)] = UKF(xk(:,i-1),Pk(:,:,i-1), yk, Rk,tk_1,tk,mu,stat);
    tk_1 = tk;
    
    % Check the error B/W the true trajectory provided and the covariance
    % matrix estimated by the filter
 
    et = tk;
    
    rv_exomars(:,i) = cspice_spkezr(target,et,frame,'NONE',observer);
    
    % Plot these two
    errorTrueTraj(:,i) = abs(xk(:,i)-rv_exomars(:,i));
    uncertaintyCovariance(:,i) = diag(sqrt(abs(Pk(:,:,i))));
end

% FIX LEGEND LABELS
figure()
hold on

% Estimated error position
plot(TIMEVECTOR(2:end),vecnorm(uncertaintyCovariance(1:3,:)),'LineWidth',1.4)

% Actual error position
plot(TIMEVECTOR(2:end),vecnorm(errorTrueTraj(1:3,:)),'LineWidth',1.4)

legendEntries = {'$\sigma_{rr,k}$',...
                '$||\vec{r}_k^{UKF}-\vec{r}_{true}||$'};
legend(legendEntries{:},'Interpreter','latex','FontSize',14)
grid on
grid minor
xLabels = cspice_et2utc(xticks(gca),'C',0);
xLabels = xLabels(:,1:17);
xticklabels(xLabels)
ylabel('Position estimation error [m]','Interpreter','latex')

figure()
hold on

% Estimated velocity error 
plot(TIMEVECTOR(2:end),vecnorm(uncertaintyCovariance(4:6,:)),'LineWidth',1.4)

% Actual velocity error
plot(TIMEVECTOR(2:end),vecnorm(errorTrueTraj(4:6,:)),'LineWidth',1.4)  

legendEntries = {'$\sigma_{vv,k}$',...
                '$||\vec{v}_k^{UKF}-\vec{v}_{true}||$'};
legend(legendEntries{:},'Interpreter','latex','FontSize',14)
grid on
grid minor
xLabels = cspice_et2utc(xticks(gca),'C',0);
xLabels = xLabels(:,1:17);
xticklabels(xLabels)
ylabel('Velocity estimation error [m/s]','Interpreter','latex')

% Plot the final position uncertainty ellipsoid
figure()
error_ellipse(real(Pk(1:3,1:3,end)),[0;0;0],'conf',0.95)
axis equal
xlabel('x [km]')
ylabel('y [km]')
zlabel('z [km]')

% ------ Point 3 -------

finalPositionUncertainty = sqrt(trace(Pk(1:3,1:3,end)));

finalVelocityUncertainty = sqrt(trace(Pk(4:6,4:6,end)));

%% Functions

function [Ymean,covariance, CHIk,Wm,Wc] = unscentedTransform(x0mean,P0,fun)
% Unscented transform
n = length(x0mean);

% Weights computation
alfa = 1e-2;
beta = 2;
k = 0;
lambda = alfa^2*(n+k)-n;

W0m = lambda/(n+lambda);
W0c = lambda/(n+lambda)+(1-alfa^2+beta);

Wi_mean = 1/(2*(n+lambda))*eye(2*n);
Wi_cov = Wi_mean;

Wm = diag([W0m; diag(Wi_mean)]);
Wc = diag([W0c; diag(Wi_cov)]);

% Generate sigma points
matrixRoot = sqrtm((n+lambda)*P0);

CHI0 = x0mean;
% CHI = zeros(n,2*n);

for i = 1:n
    CHI(:,i) = x0mean + matrixRoot(i,:)';
end

for i = (n+1):2*n
    CHI(:,i) = x0mean - matrixRoot(i-n,:)';
end

% Transform the sigma points using the nonlinear transform
    
Y0 = fun(CHI0);

% Y = zeros(length(Y0),2*n);

for j = 1:2*n
    x0 = CHI(:,j);
    Y(:,j) = fun(x0);
end

% Global matrix
Ytot = [Y0 Y]';     % 13x6 matrix containing on each row the propagation of
                    % one of the sigma points
                    
[covariance, Ymean] = weightedcov(Ytot, Wm, Wc);

CHIk = Ytot;
end

function [Ymean,covariance, CHIk,Wm,Wc] = unscentedTransform2(x0mean,P0,fun)
% Unscented transform
n = length(x0mean);

% Weights computation
alfa = 1e-2;
beta = 2;
k = 0;
lambda = alfa^2*(n+k)-n;

W0m = lambda/(n+lambda);
W0c = lambda/(n+lambda)+(1-alfa^2+beta);

Wi_mean = 1/(2*(n+lambda))*eye(2*n);
Wi_cov = Wi_mean;

Wm = diag([W0m; diag(Wi_mean)]);
Wc = diag([W0c; diag(Wi_cov)]);

% Generate sigma points
matrixRoot = sqrtm((n+lambda)*P0);

CHI0 = x0mean;

for i = 1:n
    CHI(:,i) = x0mean + matrixRoot(i,:)';
end

for i = (n+1):2*n
    CHI(:,i) = x0mean - matrixRoot(i-n,:)';
end

% Transform the sigma points using the nonlinear transform
    
Y0 = fun(CHI0);

for j = 1:2*n
    x0 = CHI(:,j);
    Y(:,j) = fun(x0);
end

% Global matrix
Ytot = [Y0 Y];     % 13x6 matrix containing on each row the propagation of
                    % one of the sigma points
                    
[covariance, Ymean] = weightedcov2(Ytot, Wm, Wc);

CHIk = Ytot;
end

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
    case {'Milano','MALARGUE','NEW_NORCIA'}
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

function [VIS,El] = visibilityOld(r,et,lat,lon,alt,minimumElevation)
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
%   VIS     logic visibility condition
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
El = El';
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
opts = odeset('reltol', 1e-12, 'abstol', 1e-14);

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

function [residual,estMeasurements] = costFcn(initialStateGuess,stations,t0,tf,mu,J2Flag,P0,x0)

tspan = t0:60:tf;

% Set options for ODE solver
opts = odeset('reltol', 1e-12, 'abstol', 1e-12);

if J2Flag == 0
    % Propagate trajectory on tspan
    sol = ode113(@twobodyode,tspan,initialStateGuess,opts,mu);
elseif J2Flag == 1
    sol = ode113(@twobodyodeJ2,tspan,initialStateGuess,opts,mu);
end

xx = deval(sol,tspan);

residual = [];

if nargin>6
    % If least squares with "a priori" append it to the residual
    W_ap = sqrtm(inv(P0));
    diff_apriori_weighted = W_ap * (initialStateGuess - x0);
    residual = [residual; diff_apriori_weighted(:)];
else
end

for k=1:3
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
    usedMeasurementsIndexes = measurementTimes < tf;
    usedMeasurements = measurements(usedMeasurementsIndexes,:);

    stationName = station.name;
    
    % Get the covariances of the measurements of the stations
    sigma_meas = sqrt(diag(station.R));
    
    % Compute the weights
    W_m = diag(1 ./ sigma_meas);
        
    weightedRes = [];
    parout = [];
    
    for i=1:length(usedMeasurements)

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
        weightedRes = [weightedRes;
                       W_m*(estimatedMeasurement'-usedMeasurements(i,:)')];
        
    end
end

function [xk,Pk] = UKF(xk_1,Pk_1, yk, Rk,tk_1,tk,mu, station)

% Transform measurements from deg to rad
yk = [wrapTo2Pi(deg2rad(yk(1)-360));
      deg2rad(yk(2));
      yk(3)];
  
%------ PREDICTION STEP ------
% Propagate sigma points with unscented tranform to tk
fun = @(x0) flow2BP(tk_1,x0,tk,mu);
[xk_apriori_mean,Pk_apriori,CHIk,Wm,Wc] = unscentedTransform2(xk_1,Pk_1,fun);

% Compute the predicted measurements of the sigma points
n = length(xk_1);

r_earth_eci = cspice_spkezr('EARTH',tk,'J2000','NONE','SUN');

for i=1:(2*n+1)
    CHIi_k = CHIk(:,i);
    
    stationName = station.name;
    rv_target_eci = CHIi_k - r_earth_eci;
    [Az,El,range] = measure(stationName, tk, rv_target_eci);
    Y_i_k = [wrapTo2Pi(Az); El; range];
    Y(:,i) = Y_i_k;
end

%------ UPDATE STEP ------
% Compute measurements covariance and cross-covariance

[Pyy_k ,Yk_mean_apriori] = weightedcov2(Y, Wm, Wc);
Pyy_k = Pyy_k + Rk;
Pxy_k = weightedxcov(Y,CHIk,Wm, Wc);

% Compute the gain
Kk = Pxy_k*inv(Pyy_k);

% Compute the a posteriori mean and covariance
% Convert the measurement to the counterclockwise azimuth convention


xk = xk_apriori_mean + (Kk*(yk - Yk_mean_apriori));
Pk = Pk_apriori - Kk*Pyy_k*(Kk');
end

function [C,Ymean] = weightedcov(Y, Wm, Wc)
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

% % Check input
% ctrl = isreal(Wm) & ~any(isnan(Wm)) & ~any(isinf(Wm));
% if ~ctrl
%   error('Check Wm: it needs be a matrix of real positive numbers with no infinite or nan values!')
% end
% 
% % Check input
% ctrl = isreal(Wc) & ~any(isnan(Wc)) & ~any(isinf(Wc));
% if ~ctrl
%   error('Check Wc: it needs be a vector of real positive numbers with no infinite or nan values!')
% end
% 
% ctrl = isreal(Y) & ~any(isnan(Y)) & ~any(isinf(Y)) & (size(size(Y), 2) == 2);
% if ~ctrl
%   error('Check Y: it needs be a 2D matrix of real numbers with no infinite or nan values!')
% end
% ctrl = length(diag(Wc)) == size(Y, 1);
% if ~ctrl
%   error('size(Y, 1) has to be equal to length(w)!')
% end

wm = diag(Wm);

[T, N] = size(Y);                   % T: number of observations;
                                    % N: number of variables
                                    

Ymean=zeros(1,N);
Py=zeros(N);

for i=1:T
Wi = Wm(i,i);
Ymean = Ymean + Wi*Y(i,:);
end

for i=1:T
Wi = Wc(i,i);
Py = Py + Wi*(Y(i,:)-Ymean)'*(Y(i,:)-Ymean);

end

C = (Py+Py')/2;
end

function [C,Xmean] = weightedcov2(X,Wm, Wc)
% INPUT:
% Y [3x13]
% X [6x13]

                                    
[N, L] = size(X);               % N: number of variables

Xmean = zeros(N,1);
Pyy = zeros(N);

for i=1:L
Wi = Wm(i,i);
Xmean = Xmean + Wi*X(:,i);
end

for i=1:L
Wi = Wc(i,i);
Pyy = Pyy + Wi*(X(:,i)-Xmean)*(X(:,i)-Xmean)';
end

C = nearestSPD((Pyy+Pyy')/2);

end

function C = weightedxcov(Y,X,Wm, Wc)
% INPUT:
% Y [3x13]
% X [6x13]

[T, L] = size(Y);                   % T: number of variables;
                                    % N: number of observations
[N, L] = size(X);

Ymean = zeros(T,1);
Xmean = zeros(N,1);
Pxy = zeros(N,T);

for i=1:L
Wi = Wm(i,i);
Ymean = Ymean + Wi*Y(:,i);
Xmean = Xmean + Wi*X(:,i);
end

for i=1:L
Wi = Wc(i,i);
Pxy = Pxy + Wi*(X(:,i)-Xmean)*(Y(:,i)-Ymean)';
end

C = Pxy;
end

function Ahat = nearestSPD(A)
% nearestSPD - the nearest (in Frobenius norm) Symmetric Positive Definite matrix to A
% usage: Ahat = nearestSPD(A)
%
% From Higham: "The nearest symmetric positive semidefinite matrix in the
% Frobenius norm to an arbitrary real matrix A is shown to be (B + H)/2,
% where H is the symmetric polar factor of B=(A + A')/2."
%
% http://www.sciencedirect.com/science/article/pii/0024379588902236
%
% arguments: (input)
%  A - square matrix, which will be converted to the nearest Symmetric
%    Positive Definite Matrix.
%
% Arguments: (output)
%  Ahat - The matrix chosen as the nearest SPD matrix to A.

if nargin ~= 1
    error('Exactly one argument must be provided.')
end

% test for a square matrix A
[r,c] = size(A);
if r ~= c
    error('A must be a square matrix.')
elseif (r == 1) && (A <= 0)
    % A was scalar and non-positive, so just return eps
    Ahat = eps;
    return
end

% symmetrize A into B
B = (A + A')/2;

% Compute the symmetric polar factor of B. Call it H.
% Clearly H is itself SPD.
[~,Sigma,V] = svd(B);
H = V*Sigma*V';

% get Ahat in the above formula
Ahat = (B+H)/2;

% ensure symmetry
Ahat = (Ahat + Ahat')/2;

% test that Ahat is in fact PD. if it is not so, then tweak it just a bit.
p = 1;
k = 0;
while p ~= 0
    [~,p] = chol(Ahat);
    k = k + 1;
    if p ~= 0
        % Ahat failed the chol test. It must have been just a hair off,
        % due to floating point trash, so it is simplest now just to
        % tweak by adding a tiny multiple of an identity matrix.
        mineig = min(eig(Ahat));
        Ahat = Ahat + (-mineig*k.^2 + eps(mineig))*eye(size(A));
    end
end
end

function [data,et,measurement] = readTDM(fileName)
% TDM files parser
% Reads the TDM file, outputing the data in a 3xN matrix, the times in a
% vector and the metadata in a measurement struct

    fid = fopen(fileName);
    tline = fgetl(fid);
    while ~strcmp(tline,'DATA_STOP')
        if strcmp(tline,'META_START')
            tline = fgetl(fid);
            while ~strcmp(tline,'META_STOP')
                msg = split(tline);
                measurement.(msg{1}) = msg{3};
                tline = fgetl(fid);
            end
        elseif strcmp(tline,'DATA_START')
            data = [];
            et = [];
            i=1;
            while ~strcmp(tline,'DATA_STOP')&&~strcmp(tline,'')
                % Get the strings from file
                fgetl(fid); % Skip empty line
                try
                    tline = fgetl(fid);
                    tline1 = fgetl(fid);
                    tline2 = fgetl(fid);
                catch
                    break
                end
                % Split the strings
                try
                    msg1 = split(tline);
                    msg2 = split(tline1);
                    msg3 = split(tline2);
                catch
                    break
                end
                % Get the three measurements
                AZ = str2double(msg1{end});
                EL = str2double(msg2{end});
                range = str2double(msg3{end});
                
                % Convert azimuth to counterclockwise convention by taking it
                % negative
                
                data(i,:) = [wrapTo360(-AZ) EL range];
                et(i) = cspice_str2et(msg1{3});
                
                i=i+1;
            end
        else
            tline = fgetl(fid);
        end
    end
end

function h = error_ellipse(varargin)
% ERROR_ELLIPSE - plot an error ellipse, or ellipsoid, defining confidence region
%    ERROR_ELLIPSE(C22) - Given a 2x2 covariance matrix, plot the
%    associated error ellipse, at the origin. It returns a graphics handle
%    of the ellipse that was drawn.
%
%    ERROR_ELLIPSE(C33) - Given a 3x3 covariance matrix, plot the
%    associated error ellipsoid, at the origin, as well as its projections
%    onto the three axes. Returns a vector of 4 graphics handles, for the
%    three ellipses (in the X-Y, Y-Z, and Z-X planes, respectively) and for
%    the ellipsoid.
%
%    ERROR_ELLIPSE(C,MU) - Plot the ellipse, or ellipsoid, centered at MU,
%    a vector whose length should match that of C (which is 2x2 or 3x3).
%
%    ERROR_ELLIPSE(...,'Property1',Value1,'Name2',Value2,...) sets the
%    values of specified properties, including:
%      'C' - Alternate method of specifying the covariance matrix
%      'mu' - Alternate method of specifying the ellipse (-oid) center
%      'conf' - A value betwen 0 and 1 specifying the confidence interval.
%        the default is 0.5 which is the 50% error ellipse.
%      'scale' - Allow the plot the be scaled to difference units.
%      'style' - A plotting style used to format ellipses.
%      'clip' - specifies a clipping radius. Portions of the ellipse, -oid,
%        outside the radius will not be shown.
%
%    NOTES: C must be positive definite for this function to work properly.

default_properties = struct(...
  'C', [], ... % The covaraince matrix (required)
  'mu', [], ... % Center of ellipse (optional)
  'conf', 0.5, ... % Percent confidence/100
  'scale', 1, ... % Scale factor, e.g. 1e-3 to plot m as km
  'style', '', ...  % Plot style
  'clip', inf); % Clipping radius

if length(varargin) >= 1 & isnumeric(varargin{1})
  default_properties.C = varargin{1};
  varargin(1) = [];
end

if length(varargin) >= 1 & isnumeric(varargin{1})
  default_properties.mu = varargin{1};
  varargin(1) = [];
end

if length(varargin) >= 1 & isnumeric(varargin{1})
  default_properties.conf = varargin{1};
  varargin(1) = [];
end

if length(varargin) >= 1 & isnumeric(varargin{1})
  default_properties.scale = varargin{1};
  varargin(1) = [];
end

if length(varargin) >= 1 & ~ischar(varargin{1})
  error('Invalid parameter/value pair arguments.') 
end

prop = getopt(default_properties, varargin{:});
C = prop.C;

if isempty(prop.mu)
  mu = zeros(length(C),1);
else
  mu = prop.mu;
end

conf = prop.conf;
scale = prop.scale;
style = prop.style;

if conf <= 0 | conf >= 1
  error('conf parameter must be in range 0 to 1, exclusive')
end

[r,c] = size(C);
if r ~= c | (r ~= 2 & r ~= 3)
  error(['Don''t know what to do with ',num2str(r),'x',num2str(c),' matrix'])
end

x0=mu(1);
y0=mu(2);

% Compute quantile for the desired percentile
k = sqrt(qchisq(conf,r)); % r is the number of dimensions (degrees of freedom)

hold_state = get(gca,'nextplot');

if r==3 & c==3
  z0=mu(3);
  
  % Make the matrix has positive eigenvalues - else it's not a valid covariance matrix!
  if any(eig(C) <=0)
    error('The covariance matrix must be positive definite (it has non-positive eigenvalues)')
  end

  % C is 3x3; extract the 2x2 matricies, and plot the associated error
  % ellipses. They are drawn in space, around the ellipsoid; it may be
  % preferable to draw them on the axes.
  Cxy = C(1:2,1:2);
  Cyz = C(2:3,2:3);
  Czx = C([3 1],[3 1]);

  [x,y,z] = getpoints(Cxy,prop.clip);
  h1=plot3(x0+k*x,y0+k*y,z0+k*z,prop.style);hold on
  [y,z,x] = getpoints(Cyz,prop.clip);
  h2=plot3(x0+k*x,y0+k*y,z0+k*z,prop.style);hold on
  [z,x,y] = getpoints(Czx,prop.clip);
  h3=plot3(x0+k*x,y0+k*y,z0+k*z,prop.style);hold on

  
  [eigvec,eigval] = eig(C);

  [X,Y,Z] = ellipsoid(0,0,0,1,1,1);
  XYZ = [X(:),Y(:),Z(:)]*sqrt(eigval)*eigvec';
  
  X(:) = scale*(k*XYZ(:,1)+x0);
  Y(:) = scale*(k*XYZ(:,2)+y0);
  Z(:) = scale*(k*XYZ(:,3)+z0);
  h4=surf(X,Y,Z);
  colormap gray
  alpha(0.3)
  camlight
  if nargout
    h=[h1 h2 h3 h4];
  end
elseif r==2 & c==2
  % Make the matrix has positive eigenvalues - else it's not a valid covariance matrix!
  if any(eig(C) <=0)
    error('The covariance matrix must be positive definite (it has non-positive eigenvalues)')
  end

  [x,y,z] = getpoints(C,prop.clip);
  h1=plot(scale*(x0+k*x),scale*(y0+k*y),prop.style);
  set(h1,'zdata',z+1)
  if nargout
    h=h1;
  end
else
  error('C (covaraince matrix) must be specified as a 2x2 or 3x3 matrix)')
end
%axis equal

set(gca,'nextplot',hold_state);
end
function [x,y,z] = getpoints(C,clipping_radius)
%---------------------------------------------------------------
% getpoints - Generate x and y points that define an ellipse, given a 2x2
%   covariance matrix, C. z, if requested, is all zeros with same shape as
%   x and y.

n=100; % Number of points around ellipse
p=0:pi/n:2*pi; % angles around a circle

[eigvec,eigval] = eig(C); % Compute eigen-stuff
xy = [cos(p'),sin(p')] * sqrt(eigval) * eigvec'; % Transformation
x = xy(:,1);
y = xy(:,2);
z = zeros(size(x));

% Clip data to a bounding radius
if nargin >= 2
  r = sqrt(sum(xy.^2,2)); % Euclidian distance (distance from center)
  x(r > clipping_radius) = nan;
  y(r > clipping_radius) = nan;
  z(r > clipping_radius) = nan;
end
end
%---------------------------------------------------------------
function x=qchisq(P,n)
% QCHISQ(P,N) - quantile of the chi-square distribution.
if nargin<2
  n=1;
end

s0 = P==0;
s1 = P==1;
s = P>0 & P<1;
x = 0.5*ones(size(P));
x(s0) = -inf;
x(s1) = inf;
x(~(s0|s1|s))=nan;

for ii=1:14
  dx = -(pchisq(x(s),n)-P(s))./dchisq(x(s),n);
  x(s) = x(s)+dx;
  if all(abs(dx) < 1e-6)
    break;
  end
end
end
%---------------------------------------------------------------
function F=pchisq(x,n)
% PCHISQ(X,N) - Probability function of the chi-square distribution.
if nargin<2
  n=1;
end
F=zeros(size(x));

if rem(n,2) == 0
  s = x>0;
  k = 0;
  for jj = 0:n/2-1;
    k = k + (x(s)/2).^jj/factorial(jj);
  end
  F(s) = 1-exp(-x(s)/2).*k;
else
  for ii=1:numel(x)
    if x(ii) > 0
      F(ii) = quadl(@dchisq,0,x(ii),1e-6,0,n);
    else
      F(ii) = 0;
    end
  end
end
end
%---------------------------------------------------------------
function f=dchisq(x,n)
% DCHISQ(X,N) - Density function of the chi-square distribution.
if nargin<2
  n=1;
end
f=zeros(size(x));
s = x>=0;
f(s) = x(s).^(n/2-1).*exp(-x(s)/2)./(2^(n/2)*gamma(n/2));
end
%---------------------------------------------------------------
function properties = getopt(properties,varargin)
%GETOPT - Process paired optional arguments as 'prop1',val1,'prop2',val2,...
%
%   getopt(properties,varargin) returns a modified properties structure,
%   given an initial properties structure, and a list of paired arguments.
%   Each argumnet pair should be of the form property_name,val where
%   property_name is the name of one of the field in properties, and val is
%   the value to be assigned to that structure field.
%
%   No validation of the values is performed.
%
% EXAMPLE:
%   properties = struct('zoom',1.0,'aspect',1.0,'gamma',1.0,'file',[],'bg',[]);
%   properties = getopt(properties,'aspect',0.76,'file','mydata.dat')
% would return:
%   properties = 
%         zoom: 1
%       aspect: 0.7600
%        gamma: 1
%         file: 'mydata.dat'
%           bg: []
%
% Typical usage in a function:
%   properties = getopt(properties,varargin{:})

% Process the properties (optional input arguments)
prop_names = fieldnames(properties);
TargetField = [];
for ii=1:length(varargin)
  arg = varargin{ii};
  if isempty(TargetField)
    if ~ischar(arg)
      error('Propery names must be character strings');
    end
    f = find(strcmp(prop_names, arg));
    if length(f) == 0
      error('%s ',['invalid property ''',arg,'''; must be one of:'],prop_names{:});
    end
    TargetField = arg;
  else
    % properties.(TargetField) = arg; % Ver 6.5 and later only
    properties = setfield(properties, TargetField, arg); % Ver 6.1 friendly
    TargetField = '';
  end
end
if ~isempty(TargetField)
  error('Property names and values must be specified in pairs.');
end
end