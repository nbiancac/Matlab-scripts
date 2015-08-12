function [stab,S0]=stab_diagramLHC_diffFD(gamma,epsnormx,epsnormy,plane, a,current_octF,current_octD,distribution)

% computes the stability diagram for the LHC (from the octupoles)
% TO BE FULLY CHECKED
%
% Input parameters:
% - gamma: relativistic mass factor of the beam,
% - epsnorm: normalized emittance (e.g. 3.75e-6 m at 7TeV),
% - a: octupole matrix as [axxF, axxD; ayyF, ayyD; axyF, axyD];
% - current_octF: current in the focusing octupoles (max is supposed to be 550A at 7TeV i.e.  gamma=7.4605e+03),
% - current_octD: current in the defocusing octupoles (max is supposed to be 550A at 7TeV i.e.  gamma=7.4605e+03)
% (most beneficial situation when sign(current_octF)=-sign(current_octD) )

% In output:
% - stabx, staby: two columns : stability diagram (i.e. Re(DeltaQ) vs. -Im(DeltaQ) )


% from Elias Metral's Mathematica notebook "StabilityDiagramsForNicolas.nb",
% June 10th, 2010.

current_max=550; % at 7 TeV

beta=sqrt(1-1/gamma^2); % relativistic velocity factor
F=current_octF/current_max*(7460.5/gamma); % reduction factor for the foc. octupole current
D=current_octD/current_max*(7460.5/gamma); % reduction factor for the defoc. octupole current

eps1sigma = eval(['epsnorm',plane])/(beta*gamma);
sigma = sqrt(eps1sigma);
%cRB = -0.65; % It is the c defined by Ruggiero and Berg in their Landau damping theory
%aRB = 270440*F;
%S0 = abs(-5*aRB*sigma^2);
% aRB = 263127*F-7309*D;
% cRB = -(87709*F-87709*D)/aRB;
% S0 = -5*aRB*sigma^2;

% NB 25-3-2015
axxF=a(1,1); axxD=a(1,2);
ayyF=a(2,1); ayyD=a(2,2);
axyF=a(3,1); axyD=a(3,2);

aRB = eval(['a',plane,plane,'F'])*F+eval(['a',plane,plane,'D'])*D;
cRB = (axyF*F+axyD*D)/aRB;

switch distribution
    case 'oldparabolic'
    
        q=[-10:0.002:10];
        S0 = -5*aRB*sigma^2;    
        TransferFunction=(4*1j/cRB)*(-(3+q.*(1+2*cRB)/cRB).*(q.^2).*log(q)+(1/((1-cRB)^2))* ...
            (((cRB+q).^3).*log(cRB+q)/cRB+(cRB-1).*(q.^2-2*cRB*q.^2-2*cRB*q-cRB)+ ...
            cRB*((1+q).^2).*(2*cRB*q-3*q-cRB).*log(1+q)));

        deltaQcoh=S0*1j./TransferFunction;
        deltaQcohRe=real(deltaQcoh).';
        deltaQcohIm=imag(deltaQcoh).';

    case 'parabolic'

        rate=1e-15;
        % 2nd order bunch distribution (cut at 3.2 sigma - see F. Ruggiero et al)
        q=[-10:0.01:10]+1j*rate;

        n=2;
        b=(n+3)*eps1sigma;
        a=(n+1)*(n+2)/(b^2);

        S0=-5*aRB*eps1sigma;

        I1=-(((cRB+q).^3.*log(1+q)-(cRB+q).^3.*log(cRB+q)+(-1+cRB)*(cRB*(cRB+2*cRB*q+(-1+2*cRB)*q.^2) + ...
            (-1+cRB)*q.^2.*(3*cRB+q+2*cRB*q).*(log(q)-log(1+q))))/(6*(-1+cRB)^2*cRB^2));

        deltaQcoh=-(aRB/(n*a*b))./I1;
        deltaQcohRe=real(deltaQcoh).';
        deltaQcohIm=imag(deltaQcoh).';

    case 'gaussian'

	rate=1e-15;
    S0=-5*aRB*eps1sigma;
    % gaussian distribution (as in Headtail)
    q=10.^[-3:1/38:2];q=[-q(end:-1:1) q]+1j*rate;


    I1=(1-cRB-(q+cRB-cRB*q).*exp(q).* expint(q) + ...
        cRB*exp(q/cRB).* expint(q/cRB)) / ((1-cRB)^2);

    deltaQcoh=-aRB*eps1sigma./I1;
    deltaQcohRe=real(deltaQcoh).';
    deltaQcohIm=imag(deltaQcoh).';

end

stab=[deltaQcohRe sign(S0)*deltaQcohIm];
stab=stab(2:end-1,:);

end