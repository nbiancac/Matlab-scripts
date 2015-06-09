MatlabDir='/afs/cern.ch/user/n/nbiancac/scratch0/Matlab/';
path(path,MatlabDir);
MMMDir='/afs/cern.ch/user/n/nbiancac/scratch0/FiniteLength/ModeMatching/';
path(path,MMMDir);
SaveDir2='/afs/cern.ch/user/n/nbiancac/scratch0/FiniteLength/Step_in+out/Results/';


%% cavity step in + out

beta=0.9525; % PS @ 2GeV
sigma=1e-12;
material={};
material.name=['Step_in+out'];
material.sigma=sigma;
material.e_r='1';
material.mu_r='1';
L=2;
b=0.035;
d=0.075;
t=d-b;
fin=[]; fstep=[]; fout=[];
fdiscrete=1e7:1e7:1e8; f_post=[];
P=5; S=25;
MMM_transverse(beta,b,t,L,material,fin,fstep,fout,fdiscrete,f_post,P,S,SaveDir2,0,0,0);
letto=dlmread([SaveDir2,'MMMdip_Re_L',num2str(L),'_Beta0.9525_b',num2str(b),'_t',num2str(t),'_Material_Step_in+out_fmin10000000_fmax100000000_P5_S',num2str(S),'.txt']);
freq=letto(:,1); re=letto(:,2);
letto=dlmread([SaveDir2,'MMMdip_Im_L',num2str(L),'_Beta0.9525_b',num2str(b),'_t',num2str(t),'_Material_Step_in+out_fmin10000000_fmax100000000_P5_S',num2str(S),'.txt']);
im=letto(:,2);
figure(39)
plot(freq,re,'r',freq,im,'k')

Z0=377;
Z_step=-Z0/2*pi*(d-b)/b^2*(d^2-b^2)/(d^2+b^2)
Z_step_NG=-Z0/2*pi*L/b^2*(d^2-b^2)/(d^2+b^2)
Z_step_MMM=im(1)