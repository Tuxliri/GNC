function STM = stmFD(t0,x0,t,mu)

STM = zeros(size(x0,1));
phi0 = phi(x0,t0,t,mu);
eps_i = sqrt(eps)*max(1,abs(x0(1:6)));


for i=1:size(x0)
    EPS_i = zeros(size(x0));
    EPS_i(i) = eps_i(i);
    phii = phi(x0 + EPS_i,t0,t,mu);
    STM(:,i) = (phii - phi0)/eps_i(i);
end

end