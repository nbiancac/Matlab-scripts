function [e,m0,c,E0]=particle_param(particle)
% Output: e [C], mass [kg], speed of light [m/s], rest energy [J] '''
if strcmp(particle,'electron')
    e=1.602176487e-19; % elementary charge
    m0=9.10938e-31; % electron mass in kg
    c=299792458; % speed of light
    E0=m0*c^2; % rest energy
elseif strcmp(particle,'proton')
    e=1.602176487e-19; % elementary charge
    m0=1.6726216e-27; % proton mass in kg
    c=299792458; % speed of light
    E0=m0*c^2; % rest energy
    
    
end