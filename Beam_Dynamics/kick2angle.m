function Dxprime=kick2angle(ktloss,machine,y0)
% computes the angle kick from a transverse kick factor.
% Dxprime=kick2angle(ktloss,machine,y0)

[q,mp,~,~]=particle_param('proton');
c=constants('clight');

Dxprime=-machine.Nb*q^2*y0*ktloss/(mp*c^2*machine.beta^2);

end