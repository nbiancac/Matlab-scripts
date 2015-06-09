function [fmaxS,Q]=Q3dB(SdB, freq, fin, fstop)

ind=find(freq<fstop & freq>fin);
freq=freq(ind);
SdB=SdB(ind);
[maxS,maxind]=max(SdB);

fmaxS=freq(maxind);
fmax=interp1(SdB(maxind:end)-maxS,freq(maxind:end),-3);
fmin=interp1(SdB(1:maxind)-maxS,freq(1:maxind),-3);

% 
% ind=(SdB<(maxS-3));
% 
% indmax=find(diff(ind)==1);
% indmin=find(diff(ind)==-1);
% fmin=freq(indmin);
% fmax=freq(indmax);
Q=fmaxS/(fmax-fmin);

end