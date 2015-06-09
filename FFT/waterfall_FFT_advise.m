%
%  waterfall_FFT_advise.m  ver 1.0  October 18, 2012
%
function[dt,df,mmm,NW,io,minf,maxf]=waterfall_FFT_advise(tim,amp)
%
tmx=max(tim);
tmi=min(tim);
%
n = max(size(amp));
dt=(tmx-tmi)/n;
%
disp(' ');
out4 = sprintf(' Time history length = %d ',n);
disp(out4)
%
NC=18;
%
ss=zeros(NC,1);
seg=zeros(NC,1);
i_seg=zeros(NC,1);
%
for i=4:NC
    ss(i) = 2^i;
    seg(i) =n/ss(i);
    i_seg(i) = fix(seg(i));
end
%
disp(' ');
out4 = sprintf(' Number of   Samples per   Time per    df');
out5 = sprintf('  Segments     Segment      Segment      ');
%
disp(out4)
disp(out5)
%
for i=4:NC
    if( i_seg(i)>0)
        str = int2str(i_seg(i));
        tseg=dt*ss(i);
        ddf=1/tseg;
        out4 = sprintf(' \t  %s  \t  %d  \t   %6.3f  \t  %6.3f',str,ss(i),tseg,ddf);
        disp(out4)
    end
end
%
disp(' ')
NW = input(' Choose the Number of Segments:  ');
disp(' ')
%
mmm = 2^fix(log(n/NW)/log(2));
%
df=1/(mmm*dt);
%
minf = input(' select minimum output frequency (Hz) ');
disp(' ');
maxf = input(' select maximum output frequency (Hz) ');
%
disp(' ');
disp(' Select Overlap ');
disp(' 1=none  2=50% ');
io=input(' ');