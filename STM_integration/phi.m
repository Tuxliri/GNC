function X = phi(t0,X0,t_end,mu)

% Set options for ODE solver
opts = odeset('Reltol',1e-13,'AbsTol',1e-14,'Stats','off');

[t, y] = ode113(@(t,y) twobodyode(t,y,mu), (t0:1:t_end), X0, opts);
X = y(end,:);
end

