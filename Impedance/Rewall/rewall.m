function [Zl,Zdip]=rewall(f,b,t,sigma,L, range)
% function [Zl,Zdip]=rewall(f,b,sigma,L, 'intermediate')
% intermediate frequency range
if strcmp(range,'intermediate')
    Zdip=constants('clight')./(2/pi./f)*(1+1i)./(pi*b^3*skin_depth(f,constants('mu0'),sigma)*sigma)*L;
    Zl=(1+1i)*L./(2*pi*b*sigma*skin_depth(f,constants('mu0'),sigma));
elseif strcmp(range,'low')
    Zdip=constants('Z0')*i1*t*2*pi*f*L/(2*pi*b*constants('clight'));
    Zl=constants('Z0')*i1*t*L/(pi*b^3);
end


end