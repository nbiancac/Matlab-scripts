function Zt=transverse_resonator(freq,Q,fr,R)
% Zt=transverse_resonator(freq,Q,fr,R)
    Zt=(R*fr./freq)./(1-1i*Q.*(fr./freq-freq./fr));
end
