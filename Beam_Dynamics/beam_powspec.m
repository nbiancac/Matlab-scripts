function distr=beam_powspec(f,sigma_t,tune,chroma,gammatr,gamma,f_rev,typedistr)
% distr=beam_powspec(f,sigma_t,tune,chroma,gammatr,gamma,f_rev,typedistr)
    if strcmp(typedistr,'Gauss')
        eta=1/gammatr^2-1/gamma^2;
        omega=2*pi*f;
        omega_rev=2*pi*f_rev;
        omega_beta=omega_rev*tune;
        omega_xi=omega_beta*chroma/eta;
        distr=exp(-((omega-omega_xi)*sigma_t).^2);
    end
end