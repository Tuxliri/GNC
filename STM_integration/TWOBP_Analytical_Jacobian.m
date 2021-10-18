syms x y z vx vy vz mu real
r = [x y z]';
v = [vx vy vz]';

X = [r; v];

f = [vx vy vz -mu*x/(x^2+y^2+z^2)^(3/2) -mu*y/(x^2+y^2+z^2)^(3/2) -mu*z/(x^2+y^2+z^2)^(3/2)]';
J = jacobian(f,X);

PHI = sym('PHI%d%d', [6 6]);
