function S=S4P(dire,type)
% S=S4P(dire,type)
letto=importdata(dire);
freq=letto.data(:,1);
freq(isnan(freq))=[];
SX1=(letto.data(:,2));
SX2=((letto.data(:,4)));
SX3=((letto.data(:,6)));
SX4=((letto.data(:,8)));

S11=SX1(1:4:end);
S21=SX1(2:4:end);
S31=SX1(3:4:end);
S41=SX1(4:4:end);

S12=SX2(1:4:end);
S22=SX2(2:4:end);
S32=SX2(3:4:end);
S42=SX2(4:4:end);

S13=SX3(1:4:end);
S23=SX3(2:4:end);
S33=SX3(3:4:end);
S43=SX3(4:4:end);

S14=SX4(1:4:end);
S24=SX4(2:4:end);
S34=SX4(3:4:end);
S44=SX4(4:4:end); 

S=eval(type);

end