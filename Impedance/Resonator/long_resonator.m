function Zl=long_resonator(freq,Q,fr,R)
% Zt=long_resonator(freq,Q,fr,R)
    Zl=(R)./(1-1i*Q.*(fr./freq-freq./fr));
end