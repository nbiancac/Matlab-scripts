function Zc=Zmatch2wires(d,D,s,type)

if strcmp(type,'2wires-shielded') % odd TEM in 2 wires shielded
    Z0=constants('Z0');
    p=s/d;
    q=s/D;
    Zc=Z0/(pi)*(log(2*p*(1-q^2)/(1+q^2))-(1+4*p^2)*(1-4*q^2)/(16*p^4));
    Zc=Zc/2;

elseif strcmp(type,'2wires') % ODD TEM for 2 wires
    Z0=constants('Z0');
    Zc=Z0/pi*acosh(s/(d));
    Zc=Zc/2;
 
end

    disp(['Zc=',num2str(Zc)])

end