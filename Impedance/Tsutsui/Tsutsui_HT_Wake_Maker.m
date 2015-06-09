% Tsutsui_HT_Wake_Makes.m

%--------------------------------------
% Computes Z[x,y][s][driv,det]
%--------------------------------------

%% Set Kicker dimension

%[a,b,c,Lz,KickerName]=[half width,half height, half hight with ferrite, name ]

TsutsuiDir='/afs/cern.ch/user/n/nbiancac/scratch0/IMPEDANCE/Tsutsui/';
KickersWakeDir='/afs/cern.ch/user/n/nbiancac/scratch0/IMPEDANCE/Tsutsui/IMP4HT/Wakes_tsutsui_PS/';
KickersImpedanceDir='/afs/cern.ch/user/n/nbiancac/scratch0/IMPEDANCE/Tsutsui/IMP4HT/Impedance_tsutsui_PS/';

%% Loop over PS Kickers


cd(TsutsuiDir);
C=importdata('PSKickerTable.txt');

for ii=1:length(C.textdata)
    
KickerName=char(C.textdata(ii));
a=C.data(ii,1)
b=C.data(ii,2)
t=C.data(ii,3)
Lz=C.data(ii,4)

%% Tsutsui's formula
cd(TsutsuiDir);
ftsu=0.01:0.01:100;

if rotate==0
[Zlongtsu,Zxtsu,Zxtsudet,Zytsu,Zytsudet]=Tsutsui(a,b,b+t,Lz,ftsu'*1e9);
elseif rotate==1
[Zlongtsu,Zytsu,Zytsudet,Zxtsu,Zxtsudet]=Tsutsui(a,b,b+t,Lz,ftsu'*1e9);
disp('rotated');
end

ftsuforFFT=[-ftsu(end:-1:1) 0 ftsu];
ZytsudetforFFT=[-conj(Zytsudet(end:-1:1))  1i*imag(Zytsudet(1)) Zytsudet];
ZytsuforFFT=[-conj(Zytsu(end:-1:1))  1i*imag(Zytsu(1)) Zytsu];
ZxtsudetforFFT=[-conj(Zxtsudet(end:-1:1))  1i*imag(Zxtsudet(1)) Zxtsudet];
ZxtsuforFFT=[-conj(Zxtsu(end:-1:1))  1i*imag(Zxtsu(1)) Zxtsu];
ZlongtsuforFFT=[-conj(Zlongtsu(end:-1:1))  1i*imag(Zlongtsu(1)) Zlongtsu];

c=299792458; 

Ntsu=size(ftsuforFFT,2);
dftsu=1e9*(ftsuforFFT(2)-ftsuforFFT(1));
dttsu=1/(Ntsu*dftsu);% sec
time=(0:1:Ntsu-1)*dttsu;
Wytsudetk=[time*c*1000;-dftsu*1e-12*Ntsu*imag(ifft(ifftshift(ZytsudetforFFT)))]';% (mm,V/Cm)
Wytsuk=[time*c*1000;-dftsu*1e-12*Ntsu*imag(ifft(ifftshift(ZytsuforFFT)))]';
Wxtsudetk=[time*c*1000;-dftsu*1e-12*Ntsu*imag(ifft(ifftshift(ZxtsudetforFFT)))]';
Wxtsuk=[time*c*1000;-dftsu*1e-12*Ntsu*imag(ifft(ifftshift(ZxtsuforFFT)))]';
Wlongtsuk=[time*c*1000;-dftsu*1e-12*Ntsu*imag(ifft(ifftshift(ZlongtsuforFFT)))]';
        
% prepare for Headtail (sign change, V/pCmm and nsec)

cd(KickersWakeDir)
       
WytsuHDTLk=[Wytsuk(:,1)/1000*1e9/c -Wytsuk(:,2)/1000];
dlmwrite([KickerName,'.Ydip.dat'],WytsuHDTLk,'\t');
WxtsuHDTLk=[Wxtsuk(:,1)/1000*1e9/c -Wxtsuk(:,2)/1000];
dlmwrite([KickerName,'.Xdip.dat'],WxtsuHDTLk,'\t');
WytsudetHDTLk=[Wytsudetk(:,1)/1000*1e9/c -Wytsudetk(:,2)/1000];
dlmwrite([KickerName,'.Yqua.dat'],WytsudetHDTLk,'\t');
WxtsudetHDTLk=[Wxtsudetk(:,1)/1000*1e9/c -Wxtsudetk(:,2)/1000];
dlmwrite([KickerName,'.Xqua.dat'],WxtsudetHDTLk,'\t');
WlongtsuHDTLk=[Wlongtsuk(:,1)/1000*1e9/c -Wlongtsuk(:,2)];

% type 4 wake for HDTLattice
wakestsu4HDTLk=[WxtsuHDTLk(:,1),WxtsuHDTLk(:,2),WytsuHDTLk(:,2),WxtsudetHDTLk(:,2),WytsudetHDTLk(:,2)];
dlmwrite([KickerName,'.wake'],wakestsu4HDTLk,'delimiter','\t');

% save impedance
cd(KickersImpedanceDir);

dlmwrite([KickerName,'.Xdip.Impedance.table.dat'],[ftsu',real(Zxtsu)', imag(Zxtsu)'],'delimiter','\t');
dlmwrite([KickerName,'.Ydip.Impedance.table.dat'],[ftsu',real(Zytsu)', imag(Zytsu)'],'delimiter','\t');
dlmwrite([KickerName,'.Xquad.Impedance.table.dat'],[ftsu',real(Zxtsudet)', imag(Zxtsudet)'],'delimiter','\t');
dlmwrite([KickerName,'.Yquad.Impedance.table.dat'],[ftsu',real(Zytsudet)', imag(Zytsudet)'],'delimiter','\t');



   
%% Plottery
%%
figure(1)
set(gca,'fontsize',12);
plot(ftsuforFFT,imag(ZytsudetforFFT),'b','LineWidth',2)
hold on;
plot(ftsuforFFT,imag(ZxtsudetforFFT),'--b','LineWidth',2)
plot(ftsuforFFT,imag(ZytsuforFFT),'r','LineWidth',2)
plot(ftsuforFFT,imag(ZxtsuforFFT),'--r','LineWidth',2)
grid on;
title([KickerName,' Transverse Impedance']);
legend('Im(Zy det)', 'Im(Zx det)', 'Im(Zy dip)','Im(Zx dip)');
xlim([0 4]);
xlabel('Frequency in GHz')
ylabel('Impedance in \Omega /m')
%saveas(gcf,['Zt_',KickerName,'.fig'],'fig');
%saveas(gcf,['Zt_',KickerName,'.pdf'],'pdf');

figure(2)
plot(ftsu,real(Zlongtsu),'LineWidth',2)
hold on
plot(ftsu,imag(Zlongtsu),'r','LineWidth',2)
legend('Tsutsui: Re(Z_{//})','Tsutsui: Im(Z_{//})');
%xlim([0 4]);
title([KickerName,' Longitudinal Impedance Zs'])
xlabel('Frequency in GHz')
ylabel('Impedance in \Omega')
grid on
% saveas(gcf,['Zs_',KickerName,'.fig'],'fig');
% saveas(gcf,['Zs_',KickerName,'.pdf'],'pdf');

figure(3)
time=WxtsudetHDTLk(:,1);
plot(time,WxtsudetHDTLk(:,2),'--r','LineWidth',2)
hold on
plot(time,WxtsuHDTLk(:,2),'r','LineWidth',2)
plot(time,WytsudetHDTLk(:,2),'--b','LineWidth',2)
plot(time,WytsuHDTLk(:,2),'b','LineWidth',2)
legend('W_{x det}','W_{x driv}', 'W_{y det}','W_{y driv}');
%xlim([0 4]);
title([KickerName,' Transverse wakes'])
xlabel('Time  [ns]')
ylabel('Wake Field  [V/(mm pC)]')
grid on
% saveas(gcf,['Wt_',KickerName,'.fig'],'fig');
% saveas(gcf,['Wt_',KickerName,'.pdf'],'pdf');

figure(4)
plot(time, WlongtsuHDTLk(:,2),'b','LineWidth',2)
legend('W_{long}');
%xlim([0 4]);
title([KickerName,' Longitudinal wake'])
xlabel('Time  [ns]')
ylabel('Wake Field  [V/(mm pC)]')
grid on
% saveas(gcf,['Ws_',KickerName,'.fig'],'fig');
% saveas(gcf,['Ws_',KickerName,'.pdf'],'pdf');
%%
cd(TsutsuiDir);

end