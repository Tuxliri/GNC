clear all
clc

mu = 398600;                % Earth's gravitational parameter   [ km^3/s^2 ]
r0 = [ -7128.137, 0, 0 ]';   % Initial radius vector             [ km ]
v0 = [ 0, -9.781, 0 ]';      % Initial velocity vector           [ km/s ]

I = eye(6);
X0 = [r0; v0; I(:)];
t0 = 0;
t = 3.8515e4;

%% Period computation
opts = odeset('Reltol',1e-13,'AbsTol',1e-14);

% tspan = linspace(t0,t,round(t-t0));
tspan = 0:1:t;
[T,XT] = ode113(@(t,y) f(t,y,mu),tspan,X0,opts);

STM1=XT(end,7:42)';
STM1 = reshape(STM1,6,6);

STM2 = stm_2(X0(1:6),@phi,t0,t,mu);


% delta_x0 = 1e-3*[100 200 300 8 10 13]';
% 
% delta_XF = STM1*delta_x0;