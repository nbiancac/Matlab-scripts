figure
set(gca,'Fontsize',16)
gifname = ['phase_beat.gif'];

for i=1:210
plot(1:1:93,Sv(i,:),'or','Markerfacecolor','b'); 
title('Phase advance beating, Vertical plane','fontsize' ,14, 'interpreter','latex');
grid on;
axis([0 93,-1,+1])
set(gca,'linewidth',2);
set(gca,'FontSize',18);
set(gca, 'FontName', 'Times New Roman');
xlabel('BPM pair')
ylabel('\Delta \Phi (rad/2\pi) ')
grid on;
%pause(0.005)
    
    I = getframe(gcf);
    I = frame2im(I);
    [X, map] = rgb2ind(I, 128);
    if i==1
        imwrite(X, map, gifname, 'GIF', 'WriteMode', 'overwrite', 'DelayTime', 0.08, 'LoopCount', Inf);
    else
        imwrite(X, map, gifname, 'GIF', 'WriteMode', 'append', 'DelayTime', 0.08);
    end
end
%%

for i=1:238
figure(1)
plot(1:1:92,Sv(i,:),'or','Markerfacecolor','b');
title('Phase advance beating, Vertical plane','fontsize' ,14, 'interpreter','latex');
xlabel('BPV Pair','fontsize' ,14, 'interpreter','latex');
ylabel('$\Delta\Phi_{i\rightarrow j}$','fontsize' ,16,'interpreter','latex')
grid on;
axis([0 92,-1,+1])
if max(Sv(i,:))>0.707 A(k)=getframe; k=k+1; end;
end;
%%
close all;
fps=20;
%B=movie(A,1,fps);
movie2avi((A),'kick_movie')

%%


fig=figure;
aviobj = avifile('example.avi','compression', 'Cinepak','quality',100)
       
for i=1:238
plot(1:1:92,Sv(i,:),'or','Markerfacecolor','b');
title('Phase advance beating, Vertical plane','fontsize' ,14, 'interpreter','latex');
xlabel('BPV Pair','fontsize' ,14, 'interpreter','latex');
ylabel('$\Delta\Phi_{i\rightarrow j}$','fontsize' ,16,'interpreter','latex')
grid on;
axis([0 92,-1,+1])
if max(Sv(i,:))>0.707 && max(Sv(i,:))<0.9 
    F = getframe(fig);
    aviobj = addframe(aviobj,F);
end;
end;
close(fig)
aviobj = close(aviobj);