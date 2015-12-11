function AdB_wire=wire_lossdB(f,rho_wire,d,mur,L,Zc)
% function AdB_wire=wire_lossdB(f,rho_wire,d,mur,L,Zc)
% Computes loss from inner conductor of a coax (wire for bench
% measurements) accounting for dependency on frequency due to skin depth
% effect. 
% f= fequency in Hz
% rho_wire=wire resistivity in ohm.m
% d= wire diameter
% mur=magnetic permeability
% L=wire length
% Zc=ch. impedance of the line

delta=skin_depth(f,mur,1/rho_wire);
R=rho_wire./(delta*pi*d); % resistance per meter
alpha=R/(2*Zc);
AdB_wire=20*alpha/log(10)*L;
end