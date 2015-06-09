%
%  waterfall_FFT_core.m  ver 1.0  October 18, 2012
%
function[store,store_p,freq_p,max_a,max_f]=...
                             waterfall_FFT_core(NW,mmm,mk,freq,amp,minf,io)
%
    store=zeros(NW,mk);
    max_a=zeros(NW,1);
    max_f=zeros(NW,1);
    sa=zeros(mmm,1);
%
    jk=1;
    for ij=1:NW
%
        clear sa;
%
        max_a(ij)=0.;
        max_f(ij)=0.;
        if(io==1)   
            for k=1:mmm
                sa(k)=amp(jk);
                jk=jk+1;
            end
        else
            for k=1:mmm
                sa(k)=amp(jk);
                jk=jk+1;
            end
            jk=jk-mmm/2;
        end
%
        clear Y;
        Y= fft(sa,mmm);
%
        store(ij,1) = abs(Y(1))/mmm; 
%   
        j=1;
        for k=2:mk 
            store(ij,k) =2.*abs(Y(k))/mmm;   
            if(freq(k)>=minf)
                store_p(ij,j)=store(ij,k);
                freq_p(j)=freq(k);
                if(store_p(ij,j)>max_a(ij))
                    max_a(ij)=store_p(ij,j);
                    max_f(ij)=freq_p(j);
                end
                j=j+1;
            end
        end    
    end