function [sat_azimuth, sat_elevation, sat_range, sat_range_rate] = antenna_pointing(stationName, et, rv_target_eci)

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