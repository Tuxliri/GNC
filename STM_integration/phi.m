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

