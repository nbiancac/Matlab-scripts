%
%  waterfall_FFT_plots.m  ver 1.0  October 18, 2012
%
function[fig_num]=...
          waterfall_FFT_plots(NW,freq_p,time_a,store_p,max_f,max_a,fig_num)
%
disp(' ')
disp(' Choose color scheme ')
color_choice = input(' 1=default 2=red 3=gray ');
%
    figure(fig_num);
    fig_num=fig_num+1;
    waterfall(freq_p,time_a,store_p);  
    set(gcf,'renderer','OpenGL' );
%
    if(color_choice ==1)
        colormap(hsv(128));
    end  
    if(color_choice ==2)
        colormap(hsv(1));
    end    
    if(color_choice ==3)
        colormap(gray);
    end
%    
    xlabel(' Frequency (Hz)');
    ylabel(' Time (sec)'); 
    zlabel(' Magnitude'); 
    view([-12 15]);   
%
    disp(' ');
    disp(' Plot as surface plot? ')
    schoice = input(' 1=yes 2=no ');
%
    while(schoice==1)
        disp(' ')
        disp(' Choose Magnitude Format');
        disp(' 1=linear 2=log ');
        imag=input(' ');
        clear store_pp;
        if(imag==1)
            store_pp=store_p;
        else
            store_pp=log10(store_p);
    %
            disp(' ');
            disp(' Select Y-axis format ');
            iyf=input(' 1=default 2=set number of decades ');
            if(iyf==2)
                disp(' ');
                disp(' Select number of decades for Y-axis ');
                idec=input(' ');
                maxP=round(max(max(store_pp))+0.5);
                minP=maxP-idec+0.0001;
                sz=size(store_pp);
                for i=1:sz(1)
                    for j=1:sz(2)
                        if(store_pp(i,j)<minP)
                            store_pp(i,j)=minP;
                        end
                    end
                end
            end
    %
        end
        disp(' ')
        disp(' Selected Shading ');
        ie=input(' 1=flat 2=faceted 3=interp ');
%
        figure(fig_num);
        fig_num=fig_num+1;
        colormap(hsv(128))
        set(gcf,'renderer','openGL');
        surf(freq_p,time_a,store_pp);
%
        if(ie==1)
            shading flat;
        end
        if(ie==2)
            shading faceted;
        end
        if(ie==3)
            shading interp;
        end
% 
        xlabel(' Frequency (Hz)');
        ylabel(' Time (sec)'); 
        zlabel(' Magnitude');   
%       axis([xmin,xmax,ymin,ymax,zmin,zmax]);
        disp(' ');
        disp(' Replot surface with different options? ');
        disp(' 1=yes 2=no ');
        ire=input(' ');
        if(ire==2)
            break;
        end
    end
%
    disp(' ');
    disp(' Plot the Spectrogram? ')
    choice = input(' 1=yes 2=no ');
    if(choice==1)
        figure(fig_num);
        fig_num=fig_num+1;
        colormap(hsv(128));
        surf(freq_p,time_a,store_p,'edgecolor','none')
        colormap(jet); axis tight;
        view(0,90);
        ylabel('Time (sec)'); xlabel('Frequency (Hz)');    
    end
%
    disp(' ');
    disp(' Plot the Contour? ')
    c_choice = input(' 1=yes 2=no ');
    if(c_choice==1)
        while(1)
            disp(' ');
            disp(' Enter the number of contour levels ');
            n=input(' ');
            figure(fig_num);
            fig_num=fig_num+1;
            contour(freq_p,time_a,store_p,n);
            ylabel('Time (sec)'); xlabel('Frequency (Hz)');  
            grid on;
            disp(' ');
            disp(' Replot contour with different options? ');
            disp(' 1=yes 2=no ');
            ire=input(' ');
            if(ire==2)
                break;
            end        
        end
    end
%
disp(' ')   
disp(' Peak Values ');
disp(' Time(sec)  Freq(Hz)  Amplitude ');
for ij=1:NW
    out4 = sprintf(' \t  %6.3f  \t  %6.3f  \t  %6.3f',time_a(ij),max_f(ij),max_a(ij));
    disp(out4) 
end