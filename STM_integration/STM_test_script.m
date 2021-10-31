clear all
clc

mu = 398600;                % Earth's gravitational parameter   [ km^3/s^2 ]
r0 = [ -7128.137, 0, 0 ]';   % Initial radius vector             [ km ]
v0 = [ 0, -3.781, 0 ]';      % Initial velocity vector           [ km/s ]

I = eye(6);
x0 = [r0; v0];
t0 = 0;
t = 500;

deltax0 = [10 0 0 0.1 0 0]';

%% Period computation
xf = phi(x0 + deltax0,t0,t,mu);
deltax = xf-x0;

STM1= stateTransitionMatrix(t0,x0,t,mu);
STM2 = stmFD(t0,x0,t,mu);

deltax_stm1 = STM1*deltax0;
deltax_stm2 = STM2*deltax0;

figure()
hold on
plot3(deltax0(1),deltax0(2),deltax0(3),'*g')
plot3(deltax(1),deltax(2),deltax(3),'*r')
plot3(deltax_stm1(1),deltax_stm1(2),deltax_stm1(3),'*k')
% [deltax deltax_stm2 deltax_stm1]