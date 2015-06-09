% function [zjm, zjmprime]=besselzeros(mbessel, nzeros);
% updated for use in Matlab, Octave, and Lyme by
% Robert J. McGough
% Michigan State University
% 14 Dec 2008
% Improved some of the explanations and the structures, plus the interface 
% has been slightly modified.
% zjm - this is the matrix of zeros of Jm(x)
% zjmprime - this is the matrix of zeros of Jm'(x)
% mbessel - positive integer describing how many bessel functions to evaluate 
% (always start with 'J0') - number of columns
% nzeros - how many zeros to calculate - number of rows
% 
% To obtain the Table on p. 393 of Blackstock, run: 
% [zjm, zjmprime]=besselzeros(5,5);
%
% 
% For further validation, see p. 408 and p. 410 in Abramowitz and Stegun
% 
% P(L): 1 - contains zeros of Jn'(x)
% P(L): 2 - contains zeros of Jn(x)
% 'the code' for M(L) - starts counting at the first nonzero entry.

% function mjdzo
%This program is a direct conversion of the corresponding Fortran program in
%S. Zhang & J. Jin "Computation of Special Functions" (Wiley, 1996).
%online: http://iris-lee3.ece.uiuc.edu/~jjin/routines/routines.html
%
%Converted by f2matlab open source project:
%online: https://sourceforge.net/projects/f2matlab/
% written by Ben Barrowes (barrowes@alum.mit.edu)
%

%     =============================================================
%     Purpose: This program computes the zeros of Bessel functions
%     Jn(x)and Jn'(x), and arranges them in the order
%     of their values
%     Input :  NT    --- Number of total zeros(NT ó 1200)
%     Output:  ZO(L)--- Value of the L-th zero of Jn(x)and
%     Jn'(x)
%     N(L)--- n, order of Jn(x)or Jn'(x)associated
%     with the L-th zero
%     M(L)--- m, serial number of the zeros of Jn(x)
%     or Jn'(x)associated with the L-th zero
%(L is the serial number of all the
%     zeros of Jn(x)and Jn'(x))
%     P(L)--- 1(TE)or 2(TM), a code for designating the
%     zeros of Jn(x)or Jn'(x)
%     In the waveguide applications, the zeros
%     of Jn(x)correspond to TM modes and those
%     of Jn'(x)correspond to TE modes.
%     =============================================================
% nt=[];n=[];m=[];p=[];zo=[];
%  n=zeros(1,1400);
% m=zeros(1,1400);
% zo=zeros(1,1400+1);
% fprintf(1,'%s \n','nt=?');
%     READ(*,*)NT
% nt=5;
% fprintf(1,[repmat(' ',1,1),'total number of the zeros:','%5g' ' \n'],nt);
% fprintf(1,[repmat(' ',1,15),'***  please wait.  the program is running  ***' ' \n']);
% [nt,n,m,p,zo]=jdzo(nt,n,m,p,zo);
% fprintf(1,'%0.15g \n');
% ks=nt./101+1;
% for  k0=1:ks;
% fprintf(1,'%s ',' table           zeros of bessel');fprintf(1,'%s \n', ' functions jn(x)and jn''(x)');
% fprintf(1,'%0.15g \n');
% fprintf(1,'%s ',' ----------------------------------');fprintf(1,'%s \n','----------------------------------');
% for  k=1:50;
% j1=100.*(k0-1)+k+1;
% j2=j1+50;
% if(j1 <= nt+1&j2 <= nt+1);
% fprintf(1,[repmat(' ',1,1),'%4g',repmat(' ',1,3),'%2s','%4g',' -','%2g','%14.8g',repmat(' ',1,3),'|',repmat(' ',1,2),'%4g',repmat(' ',1,3),'%2s','%4g',' -','%2g','%14.8g' ' \n'],j1-1,p(j1),n(j1),m(j1),zo(j1+1), j2-1,p(j2),n(j2),m(j2),zo(j2+1));
% elseif(j1 <= nt+1&j2 > nt+1);
% fprintf(1,[repmat(' ',1,1),'%4g',repmat(' ',1,3),'%2s','%4g',' -','%2g','%14.8g',repmat(' ',1,3),'|',repmat(' ',1,2),'%4g',repmat(' ',1,3),'%2s','%4g',' -','%2g','%14.8g' ' \n'],j1-1,p(j1),n(j1),m(j1),zo(j1+1));
% end;
% end;  k=50+1;
% fprintf(1,'%s ',' ----------------------------------');fprintf(1,'%s \n','----------------------------------');
% fprintf(1,[ '\n '  ' \n']);
% end;  k0=ks+1;
%format(1x,['total number of the zeros:'],i5);
%format(1x,i4,3x,a2,i4,2h -,i2,f14.8,3x,1h|,2x,i4, 3x,a2,i4,2h -,i2,f14.8);
%format(15x,'***  please wait.  the program is running  ***');
%format[);
% end
% function [zo, n, m, p]=jdzo(nt);
function [zjm, zjmprime] = besselzeros2(mbessel, nzeros);
%     ===========================================================
%     Purpose: Compute the zeros of Bessel functions Jn(x)and
%     Jn'(x), and arrange them in the order of their
%     magnitudes
%     Input :  NT    --- Number of total zeros(NT ó 1200)
%     Output:  ZO(L)--- Value of the L-th zero of Jn(x)
%     and Jn'(x)
%     N(L)--- n, order of Jn(x)or Jn'(x)associated
%     with the L-th zero
%     M(L)--- m, serial number of the zeros of Jn(x)
%     or Jn'(x)associated with the L-th zero
%(L is the serial number of all the
%     zeros of Jn(x)and Jn'(x))
%     P(L)--- TM or TE, a code for designating the
%     zeros of Jn(x)or Jn'(x).
%     In the waveguide applications, the zeros
%     of Jn(x)correspond to TM modes and
%     those of Jn'(x)correspond to TE modes
%     Routine called:    BJNDD for computing Jn(x), Jn'(x)and
%     Jn''(x)
%     =============================================================
i=[];x=[];bj=[];dj=[];fj=[];
%n1=zeros(1,70);
%m1=zeros(1,70);
%zoc=zeros(1,70+1);

bj=zeros(1,101);
dj=zeros(1,101);
fj=zeros(1,101);
x=0.0;
%if(nt < 600),
%    xm=-1.0+2.248485.*fix(nt).^0.5-.0159382.*fix(nt)+3.208775e-4 .*fix(nt).^1.5;
%    nm=fix(14.5+.05875.*fix(nt));
%    mm=fix(.02.*fix(nt))+6;
%else
%    xm=5.0+1.445389.*fix(nt).^.5+.01889876.*fix(nt)-2.147763e-4 .*fix(nt).^1.5;
%    nm=fix(27.8+.0327.*fix(nt));
%    mm=fix(.01088.*fix(nt))+10;
%end


% the code uses an indexing scheme that is reversed relative to Blackstock.
% so, swap them here
nm = mbessel;
mm = nzeros;

zjm = zeros(mm, nm);
zjmprime = zeros(mm, nm);

%l0=0;
for  i=1:nm;
    x1=.407658+.4795504.*(i-1).^.5+.983618.*(i-1);
    x2=1.99535+.8333883.*(i-1).^.5+.984584.*(i-1);
 %   l1=0;
    for  j=1:mm;
        ifoo=1;
        if(~(i == 1 & j == 1));
        x=x1;
            while (1);
            	[i,x,bj,dj,fj]=bjndd(i,x,bj,dj,fj);
                x0=x;
                x=x-dj(i)./fj(i);
 %               if(x1 > xm);
 %                   ifoo=0;
 %                   break;
 %               end;
                if(~(abs(x-x0)> 1.0d-10)),
                    break;
                end;
            end;
        end;
        if(ifoo == 1);
 %           l1=l1+1;
 %           n1(l1)=i-1;
 %           m1(l1)=j;
 %           if(i == 1),
 %               m1(l1)=j-1;
 %           end;
 %           p1(l1)=1;
 %           zoc(l1+1)=x;
            
            zjmprime(j, i) = x;
            
            if(i <= 15);
                x1=x+3.057+.0122.*(i-1)+(1.555+.41575.*(i-1))./(j+1).^2;
            else
                x1=x+2.918+.01924.*(i-1)+(6.26+.13205.*(i-1))./(j+1).^2;
            end;
        end;

        x=x2;
        ifoo=1;
        while (1);
            [i,x,bj,dj,fj]=bjndd(i,x,bj,dj,fj);
            x0=x;
            x=x-bj(i)./dj(i);
 %           if(x > xm);
 %               ifoo=0;
 %               break;
 %           end;
            if(~(abs(x-x0)> 1.0d-10)),
                break;
            end;
        end;
        if(ifoo == 1);
 %           l1=l1+1;
 %           n1(l1)=i-1;
 %           m1(l1)=j;
 %           p1(l1)=2;
 %           zoc(l1+1)=x;
            
            zjm(j, i) = x;
            
            if(i <= 15);
                x2=x+3.11+.0138.*(i-1)+(.04832+.2804.*(i-1))./(j+1).^2;
            else
                x2=x+3.001+.0105.*(i-1)+(11.52+.48525.*(i-1))./(j+3).^2;
            end;
        end;
    end;
    

  %  l=l0+l1;
  %  l2=l;
    
  %  while (1);
  %      if(l0 == 0),
  %          for  k=1:l;
  %             zo(k+1)=zoc(k+1);
  %              n(k)=fix(n1(k));
  %              m(k)=fix(m1(k));
  %              p(k)=p1(k);
  %          end;  
%            k=l+1;
  %          l1=0;
  %      elseif(l0 ~= 0);

%            if(zo(l0+1)>= zoc(l1+1));
%                zo(l0+l1+1)=zo(l0+1);
%                n(l0+l1)=fix(fix(n(l0)));
%                m(l0+l1)=fix(fix(m(l0)));
%                p(l0+l1)=p(l0);
%                l0=l0-1;
%            else
 %               zo(l0+l1+1)=zoc(l1+1);
 %               n(l0+l1)=fix(n1(l1));
 %               m(l0+l1)=fix(m1(l1));
 %               p(l0+l1)=p1(l1);
%                l1=l1-1;
 %           end;
 %       end;
 %       if(~(l1 ~= 0))
 %           break;
 %       end;
  %  end;
  %  l0=l2;
end;

% trim off extra leading zero s.t. the resulting arrays are all the same length
% nzo = length(zo);
% zo = zo(2:nzo);

return;
end
function [n,x,bj,dj,fj]=bjndd(n,x,bj,dj,fj,varargin);
%     =====================================================
%     Purpose: Compute Bessel functions Jn(x)and their
%     first and second derivatives(n= 0,1,úúú)
%     Input:   x ---  Argument of Jn(x,x ò 0)
%     n ---  Order of Jn(x)
%     Output:  BJ(n+1)---  Jn(x)
%     DJ(n+1)---  Jn'(x)
%     FJ(n+1)---  Jn"(x)
%     =====================================================
for  nt=1:900;
mt=fix(0.5.*log10(6.28.*nt)-nt.*log10(1.36.*abs(x)./nt));
if(mt > 20),
    break;
end;
end;
m=nt;
bs=0.0d0;
f0=0.0d0;
f1=1.0d-35;
for  k=m:-1:0;
f=2.0d0.*(k+1.0d0).*f1./x-f0;
if(k <= n)bj(k+1)=f; end;
if(k == 2.*fix(k./2))bs=bs+2.0d0.*f; end;
f0=f1;
f1=f;
end;  k=0-1;
for  k=0:n;
bj(k+1)=bj(k+1)./(bs-f);
end;  k=fix(n)+1;
dj(1)=-bj(2);
fj(1)=-1.0d0.*bj(1)-dj(1)./x;
for  k=1:n;
dj(k+1)=bj(k)-k.*bj(k+1)./x;
fj(k+1)=(k.*k./(x.*x)-1.0d0).*bj(k+1)-dj(k+1)./x;
end;  k=fix(n)+1;
return;
end

