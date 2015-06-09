function S=S2P(dire,type)
% S=S2P(dire,type)
    letto=importdata(dire);
    freq=letto.data(:,1);
    S11=(letto.data(:,2).*exp(1i.*letto.data(:,3)*pi/180));
    S21=(letto.data(:,4).*exp(1i.*letto.data(:,5)*pi/180));
    S12=(letto.data(:,6).*exp(1i.*letto.data(:,7)*pi/180));
    S22=(letto.data(:,8).*exp(1i.*letto.data(:,9)*pi/180));
    S=eval(type);
end