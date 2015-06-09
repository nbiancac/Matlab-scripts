function [a,qsec]=octupole(DataDir)

% computes the maximum detuning coefficients from the octupoles and maximum
% Q'' due to the octupoles, at mom [in GeV/c] (these are inversely proportional to the energy).

% a is a matrix 3*2: [axF  axD;
%                     ayF  ayD;
%                     axyF axyD]
% Those are the detuning coefficient with maximum current (550A) in the focusing
% octupoles and zero in the defocusing ones (with letter F), or zero in the focusing ones and maximum (550A) 
% in the defocusing octupoles (with letter D).
% components ax are multiplied by Jx and detune Qx
% components ay are multiplied by Jy and detune Qy
% components axy are multiplied by Jy (resp. Jx) and detune Qx (resp. Qy)

% qsec is a matrix 2*2: [Q''xF Q''xD;
%                        Q''yF Q''yD]
% Those are the Q'' at 7TeV, with maximum current (550A) in the focusing
% octupoles and zero in the defocusing ones (with letter F), or zero in the focusing ones and maximum (550A)
% in the defocusing octupoles (with letter D).



O3=63100; % maximum absolute octupolar strength in T/m^3 (from MAD-X)
e=particle_param('proton');
mom=7000*1e9;
K3=6*O3/(mom/(constants('clight'))); % K3+ (maximum normalized octupolar strength)


fid=fopen([DataDir,'octupole.dat'],'r');
letto=textscan(fid,'%s%f%f%f%f%f%f','headerlines',1); fclose(fid);
oct_length=cell2mat(letto(:,3));
betax=cell2mat(letto(:,4));
betay=cell2mat(letto(:,5));
dx=cell2mat(letto(:,6));
dy=cell2mat(letto(:,7));

indF=find(abs(betax-max(betax))<10); % focusing octupoles
indD=find(abs(betax-min(betax))<10); % defocusing octupoles

% note: additional minus sign for the defocusing octupoles because O3D=-O3F for
% the same current in foc. and defoc. octupoles

axF=sum(oct_length(indF).*betax(indF).^2)*K3/(16*pi);
axD=-sum(oct_length(indD).*betax(indD).^2)*K3/(16*pi);
ayF=sum(oct_length(indF).*betay(indF).^2)*K3/(16*pi);
ayD=-sum(oct_length(indD).*betay(indD).^2)*K3/(16*pi);
axyF=-sum(oct_length(indF).*betax(indF).*betay(indF))*K3/(8*pi);
axyD=sum(oct_length(indD).*betax(indD).*betay(indD))*K3/(8*pi);

a=[axF  axD; ayF  ayD; axyF axyD];
% disp(a)
dlmwrite([DataDir,'a.dat'],a);

qsecxF=sum(oct_length(indF).*betax(indF).*dx(indF).^2)*K3/(4*pi);
qsecxD=-sum(oct_length(indD).*betax(indD).*dx(indD).^2)*K3/(4*pi);
qsecyF=-sum(oct_length(indF).*betay(indF).*dx(indF).^2)*K3/(4*pi);
qsecyD=sum(oct_length(indD).*betay(indD).*dx(indD).^2)*K3/(4*pi);

qsec=[qsecxF qsecxD; qsecyF qsecyD];
% disp(qsec)
dlmwrite([DataDir,'qsec.dat'],qsec);    

end