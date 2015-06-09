function [freq,Z]=sum_imp(freq1,Z1,freq2,Z2)
if size(freq1,1)>size(freq1,2)
    freq1=freq1';
    Z1=conj(Z1');
end
if size(freq2,1)>size(freq2,2)
    freq2=freq2';
    Z2=conj(Z2');
end    
freq=sort(unique([freq1,freq2]));
Z1_int=interp1(freq1,Z1,freq,'linear','extrap');
Z2_int=interp1(freq2,Z2,freq,'linear','extrap');
Z=Z1_int+Z2_int;

end