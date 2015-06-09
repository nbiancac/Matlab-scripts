function machine=PS_param(E0,E)

machine={};

if (E==2.32e9); machine.E=2.32; machine.V=2e5;machine.Estr='2.3GeV';taub_full=120e-9; % full length in s
elseif (E==2.93e9); machine.E=2.93; machine.V=2e5;machine.Estr='2.9GeV';taub_full=120e-9; % full length in s
elseif (E==8.19e9); machine.E=8.19; machine.V=2e5; machine.Estr='8.19GeV';taub_full=55e-9; % full length in s
elseif (E==14e9); machine.E=14; machine.V=2e5; machine.Estr='14GeV';taub_full=55e-9; % full length in s
elseif (E==26.41e9); machine.E=26.41; machine.V=2e5; machine.Estr='26GeV';taub_full=45e-9; % full length in s
else
        error('PS energy not in database');
end

e=constants('e');
c=constants('clight');
mp=constants('mp');
machine.name='PS';
machine.h=8; 
machine.gamma=E*e/E0;
machine.beta=sqrt(1.-1./(machine.gamma^2));
machine.mom=machine.gamma*mp*machine.beta*c/(1e9*e/c); %GeV/c
machine.Ekin=(E-E0/e)/1e9; %GeV
machine.sigmat=taub_full/4; % rms bunch length in s
machine.taub=machine.sigmat; % rms bunch length in s
machine.sigmaz=machine.sigmat*machine.beta*c; % rms bunch length
machine.circ=628.31853072; % total circumference in m
machine.R=machine.circ/(2.*pi); % machine radius
machine.Qx=6.13;machine.Qxfrac=machine.Qx-floor(machine.Qx);
machine.Qy=6.26;machine.Qyfrac=machine.Qy-floor(machine.Qy);
machine.alphap=0.02776; % momentum compaction factor
machine.eta=(machine.alphap-1./(machine.gamma^2)); % slip factor
machine.gammatr=sqrt(abs(1/machine.alphap));
machine.phis=0;
machine.Qs=Qs_from_RF_param(machine.V,machine.h,machine.gamma,abs(machine.eta),machine.phis,'proton');
machine.f0=c*machine.beta/machine.circ; % rev. frequency
machine.omega0=2.*pi*machine.f0;
machine.omegas=machine.Qs*machine.omega0;
machine.dphase=0.; 
machine.chromay=0;
machine.chromax=0;
machine.Nb=1e12;
machine.M=1;
end

