% lab02_nbody_propagator.m - N-body propagator with SPICE.
%
% Description:
%   The script integrates the motion of an asteroid in the Solar system.
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
% Prerequisites:
%   - MICE (Matlab SPICE)
%   - Meta-kernel in current directory
%   - Kernels in subfolder 'kernels' (PCK, SPK, LSK)
%
clc, clear all, close all

%% Initialization

% Load kernel
cspice_furnsh('sgn_nbody_propagator.tm');

% Define list of celestial bodies:
labels = {'Sun';
          'Mercury';
          'Venus';
          'Earth';
          'Moon';
          'Mars Barycenter';
          'Jupiter Barycenter';
          'Saturn Barycenter';
          'Uranus Barycenter';
          'Neptune Barycenter';
          'Pluto Barycenter'};

% Initialize propagation data
bodies = nbody_init(labels);

% select integration frame string (SPICE naming convention)
frame = 'J2000';

%% Define initial state of object
% Values taken from NASA Horizon
% (https://ssd.jpl.nasa.gov/horizons.cgi#top)
%
%*******************************************************************************
%Target body name: 65803 Didymos (1996 GT)         {source: JPL#181}
%Center body name: Solar System Barycenter (0)     {source: DE441}
%Center-site name: BODY CENTER
%*******************************************************************************
%Start time      : A.D. 2021-Sep-26 00:00:00.0000 TDB
%Stop  time      : A.D. 2021-Oct-26 00:00:00.0000 TDB
%Step-size       : 1440 minutes
%*******************************************************************************
%2459485.500000000 = A.D. 2021-Sep-28 00:00:00.0000 TDB 
% X =-2.905136253172431E+08 Y =-1.688522330052698E+08 Z =-5.856187743850798E+07
% VX= 7.991281819510085E+00 VY=-1.190752125042429E+01 VZ=-5.908031506561663E+00

ref_epoch_str = '2021-Sep-28 00:00:00.0000 TDB';

et0 = cspice_str2et(ref_epoch_str);
x0 = [-2.905136253172431E+08; -1.688522330052698E+08; -5.856187743850798E+07;
       7.991281819510085E+00; -1.190752125042429E+01; -5.908031506561663E+00];

%*******************************************************************************
%Target body name: 21 Lutetia (A852 VA)            {source: JPL#118}
%Center body name: Solar System Barycenter (0)     {source: DE441}
%Center-site name: BODY CENTER
%*******************************************************************************
%2459485.500000000 = A.D. 2021-Sep-28 00:00:00.0000 TDB 
% X =-4.159078677837126E+08 Y = 1.993906728979437E+07 Z = 3.283807596788065E+07
% VX= 9.072577190190879E-02 VY=-1.511736475258038E+01 VZ=-6.712713574158217E+00
% x0 = [-4.159078677837126E+08; 1.993906728979437E+07; 3.283807596788065E+07;
%       9.072577190190879E-02; -1.511736475258038E+01; -6.712713574158217E+00];

%% Perform integration
final_epoch_str = '2025-Sep-28 00:00:00.0000 TDB';
etf = cspice_str2et(final_epoch_str);

options = odeset('reltol', 1e-12, 'abstol', 1e-12);
[tt, xx] = ode113(@(t,x) nbody_rhs(t,x,bodies,frame), [et0 etf], x0, options);

%% Plots
figure(1)
plot3(xx(:,1), xx(:,2), xx(:,3),'color','#D95319');
xlabel('x [AU]')
ylabel('y [AU]')
zlabel('z [AU]')
title (['@SSB:',frame])
axis equal
hold on

others = {'Earth', 'Mars Barycenter', 'Venus Barycenter'};
planet_colors = {'#77AC30','#A2142F','#EDB120'};
for i = 1:length(others)
    rr_planet = zeros(3, length(tt));
    for j = 1:length(tt)
        rr_planet(:,j) = cspice_spkpos(others{i},tt(j),frame,'NONE','SSB');
    end
    plot3(rr_planet(1,:), rr_planet(2,:), rr_planet(3,:),'color',planet_colors{i});
end
legend(['65803 Didymos',others],'Location','best');

%% Plot Earth distance
dist = zeros(1,length(tt));

for j = 1:length(tt)
    rr_earth_obj = xx(j,1:3)' ...
                   - cspice_spkpos('Earth',tt(j),frame,'NONE','SSB');
               
    dist(j) = cspice_convrt(norm(rr_earth_obj),'km','au');
end

figure(2)
plot(tt/cspice_spd(), dist,'k')
title('Relative distance Earth-Dydimos')
xlabel('Epoch [MJD2000]')
ylabel('Distance [AU]')

%% Clear kernel
cspice_kclear()