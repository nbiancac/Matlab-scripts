function s=constants(var)
% various constants: e.g. s=constants('clight')

if strcmp(var,'clight')
    s=299792456.2;
elseif strcmp(var,'Z0')
    s=120*pi;
elseif strcmp(var,'mu0')
    s=4*pi*10^-7;
elseif strcmp(var,'epsilon0')
    s=(4*pi*10^-7)/(120*pi)^2;    
elseif strcmp(var,'mp')
    s=1.67262178e-27;
elseif strcmp(var,'e')
    s=1.60217657e-19;
end





end