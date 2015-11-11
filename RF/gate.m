 function [pulse21_gated,S21_gated]=gate(space_pulse,pulse21,sigma_t,s1,s2)
%      function [pulse21_gated,S21_gated]=gate(space_pulse,pulse21,s1,s2)
    ind_gate=[space_pulse>s1 & space_pulse<s2];
    
    space_pulse_gated=space_pulse;
    pulse21_gated=zeros(size(pulse21));
    
    pulse21_gated(ind_gate)=pulse21(ind_gate);
    
    
    
    
    Nfull=length(space_pulse_gated);
    
    dt=diff(space_pulse_gated(1:2))/constants('clight');
    df=1/(dt*Nfull);
    freq_full=df.*[(-Nfull+1)/2:1:(Nfull-1)/2];
    wind_freq=exp(-0.5*(2*pi*freq_full.*sigma_t).^2);
    
    S21_gated_all=(fftshift(fft(fftshift(pulse21_gated'))))./wind_freq;
    
    S21_gated=S21_gated_all(freq_full>=0);
 end