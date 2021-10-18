function PHI = stm_2(X0,f,t0,t,mu)

STM = zeros(size(X0,1));
phi0 = f(t0,X0,t,mu);
eps_i = sqrt(eps)*max(1,abs(X0(1:6)));

for i=1:size(X0)
    EPS_i = zeros(size(X0));
    EPS_i(i) = eps_i(i);
    phii = f(t0,X0 + EPS_i,t,mu);
    STM(:,i) = (phii - phi0)'/eps_i(i);
end

PHI = STM;

end