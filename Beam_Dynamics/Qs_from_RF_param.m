function Qs=Qs_from_RF_param(V,h,gamma,eta,phis,particle)

    % computes Qs (at zero amplitude) from RF parameters:
    % - V = RF voltage [V],
    % - h = RF harmonic number,
    % - gamma = relativistic mass factor,
    % - eta = slip factor = alphap - 1/gamma^2,
    % - phis = synchrotron phase [rad],
    % - particle -> 'proton' or 'electron'.
    
    [e,m0,c,E0]=particle_param(particle);
    beta=sqrt(1.-1./(gamma^2));
    p0=m0*beta*gamma*c;
    Qs=sqrt(e*V*eta*h*cos(phis)/(2*pi*beta*c*p0));

end