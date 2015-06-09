function dQdN=tune_shift_m0(Zeff,beta_func,gamma, R,sigma_z,f0,particle)
    % dQdN=tune_shift_m0(Zeff,beta_func,gamma,R,sigma_z,f0,particle)
    %
    % compute sacherer mode 0 tune shift for gaussian mode.
    % It takes into account the beta_func doing betaY/beta_av where beta_av is
    % the smooth approx one: beta_av=R/Q;
    
    [e,mp,c,E0]=particle_param(particle);
    T0=1/f0;
    dQdN=-Zeff*beta_func/((4*sqrt(pi)*gamma*R*mp*(2*pi)^2*sigma_z)/(e^2*T0));
end        