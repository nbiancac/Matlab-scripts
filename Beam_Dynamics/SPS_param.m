function SPS=SPS_param(E0,E,optics)


if E==26e9
        SPS.Estr='26GeV';
        SPS.Nb=1.3e11;
        SPS.M=288;
        SPS.V=2e6;
        SPS.scenario='Injection';
        SPS.taub=3e-9/4;
elseif E==55e9

        SPS.Estr='55GeV';
        SPS.Nb=1.6e11;
        SPS.M=1;
        SPS.V=4.6e6;
        SPS.scenario='55GeV plateau';
        SPS.taub=3.5e-9/4;
elseif E==270e9

        SPS.Estr='270GeV';
        SPS.Nb=1.3e11;
        SPS.M=288;
        SPS.V=4.5e6;
        SPS.scenario='270GeV plateau';
        SPS.taub=1.6e-9/4;

elseif E==450e9

        SPS.Estr='450GeV';
        SPS.Nb=1.3e11;
        SPS.M=288;
        SPS.V=7e6;
        SPS.scenario='Extraction';
        SPS.taub=1.6e-9/4;
end


e=1.602176487e-19;
c=299792458;
SPS.machine='SPS';
SPS.h=4620; 
SPS.gamma=E*e/E0;
SPS.beta=sqrt(1.-1./(SPS.gamma^2));
SPS.sigmaz=SPS.taub*SPS.beta*c; % rms bunch length
SPS.circ=6.9e3; % total circumference in m
SPS.R=SPS.circ/(2.*pi); % machine radius

if strcmp(optics,'Q20')
    SPS.Qx=20.13;SPS.Qxfrac=SPS.Qx-floor(SPS.Qx);
    SPS.Qy=20.18;SPS.Qyfrac=SPS.Qy-floor(SPS.Qy);
    SPS.gammatr=18;
elseif strcmp(optics,'Q26')
    SPS.Qx=26.13;SPS.Qxfrac=SPS.Qx-floor(SPS.Qx);
    SPS.Qy=26.18;SPS.Qyfrac=SPS.Qy-floor(SPS.Qy);
    SPS.gammatr=22.8;
end    
SPS.alphap=1/SPS.gammatr^2;
SPS.eta=SPS.alphap-1./(SPS.gamma^2); % slip factor
SPS.phis=0;
SPS.Qs=Qs_from_RF_param(SPS.V,SPS.h,SPS.gamma,SPS.eta,SPS.phis,'proton');
SPS.f0=c*SPS.beta/SPS.circ; % rev. frequency
SPS.omega0=2.*pi*SPS.f0;
SPS.omegas=SPS.Qs*SPS.omega0;
SPS.dphase=0.; 
SPS.chromax=0;
SPS.chromay=0;
SPS.betavy=SPS.circ/floor(SPS.Qy)/2/pi;
SPS.betavx=SPS.circ/floor(SPS.Qx)/2/pi;
end