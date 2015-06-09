function delta=skin_dept(f,mur,sigma)
% function delta=skin_dept(f,mur,sigma)
    mu=4*pi*10^-7*mur;
    omega=2*pi*f;
    delta=sqrt(2/(omega*mu*sigma));
    disp(['delta=',num2str(delta*1e3),'mm']);
end