function CC=sumCCmodel(CC1,CC2)
% function CC=sumCCmodel(CC1,CC2)
%
% Sums two different CC HOMs sets.


    freq=sort(unique([CC1.freq,CC2.freq]));
    Z1=interp1(CC1.freq,CC1.Z,freq);
    Z2=interp1(CC2.freq,CC2.Z,freq);
    Z=Z1+Z2;
    CC.freq=freq; CC.Z=Z;

    ind=isnan(CC.Z);
    
    if ~isempty(ind)
        disp('cleaning nans..')
        CC.Z(ind)=[];
        CC.freq(ind)=[];
    end
end