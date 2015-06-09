%%%%%%%%%%% 
% Movie
%%%%%%%%%%%


set(0,'defaultlinelinewidth',1)
set(0,'defaultaxesfontname','timesnewroman');
set(0,'defaulttextfontname','timesnewroman');
set(0,'DefaultAxesFontSize', 14) 
set(0,'DefaultTextFontSize', 14) 
set(0,'defaulttextinterpreter','TeX');


mainDir='/afs/cern.ch/user/n/nbiancac/scratch0/FiniteLength/ModeMatching/Transverse_simulations/';
cd(mainDir);

f=1757000000;
b=0.05;
L=0.20;
t=0.25;
c=b+t;
P_vec=[5];


for P=P_vec
    S=P
        
        zz_mesh=dlmread(['MMMdip_MeshGridZ_L0.2_Beta1_b0.05_t0.25_phi0.7854_Material_test_f1757000000_P5_S5.txt']);
        rr_mesh=dlmread(['MMMdip_MeshGridR_L0.2_Beta1_b0.05_t0.25_phi0.7854_Material_test_f1757000000_P5_S5.txt']);
        PEz=dlmread(['MMMdip_Ez_L0.2_Beta1_b0.05_t0.25_phi0.7854_Material_test_f1757000000_P5_S5.txt']);
        PEr=dlmread(['MMMdip_Er_L0.2_Beta1_b0.05_t0.25_phi0.7854_Material_test_f1757000000_P5_S5.txt']);
        PEphi=dlmread(['MMMdip_Ephi_L0.2_Beta1_b0.05_t0.25_phi0.7854_Material_test_f1757000000_P5_S5.txt']);
        PHz=dlmread(['MMMdip_Hz_L0.2_Beta1_b0.05_t0.25_phi0.7854_Material_test_f1757000000_P5_S5.txt']);
        PHr=dlmread(['MMMdip_Hr_L0.2_Beta1_b0.05_t0.25_phi0.7854_Material_test_f1757000000_P5_S5.txt']);
        PHphi=dlmread(['MMMdip_Hphi_L0.2_Beta1_b0.05_t0.25_phi0.7854_Material_test_f1757000000_P5_S5.txt']);
        Dzz=zz_mesh(2)-zz_mesh(1);
        Drr=rr_mesh(2)-rr_mesh(1);

        % Movie Ephi


            zz_sx=[-L-Dzz/2:Dzz:Dzz];
            zz_dx=[L-Dzz/2:Dzz:2*L];
            zz_c=[0+Dzz/4:Dzz/2:L];
            rr_sx=[b:Drr:c];

            omega=2*pi*f;
            figure()
            writerObj = VideoWriter(['MovieEphi_P',num2str(P),'_S',num2str(S),'.avi']);
            open(writerObj);
            axis tight
            set(gca,'nextplot','replacechildren');
            set(gcf,'Renderer','zbuffer');
            N=3; % # periods displayed
            T=1/f; % period from frequency
            V0=(PEphi);
            Nframe=60; % # frames per seconds.
            crange=[-max(max(real(V0))) max(max(real(V0)))];
            colormap(jet(100));
            time_vec=-N*T/2:T/(N*Nframe):N*T/2;
            for ii=1:length(time_vec)
                tt=time_vec(ii);
                V=real(V0.*exp(1i*omega*tt));
                pcolor(zz_mesh,rr_mesh,V); shading interp; hold on; contour(zz_mesh,rr_mesh,V); 
                xlabel('Length [m]'); ylabel('Radius [m]'); title('Field E_{z} [V/m]');caxis(crange);colorbar;
                pcolor(zz_sx,[b-Drr/2,b+Drr/16],ones(length(zz_sx),2)'*(-max(max(real(V0)))));
                pcolor(zz_dx,[b-Drr/2,b+Drr/16],ones(length(zz_dx),2)'*(-max(max(real(V0)))));
                pcolor([0+Dzz/4,0+Dzz/2],rr_sx,ones(2,length(rr_sx))'*(-max(max(real(V0)))));
                pcolor([L-Dzz/4,L-Dzz/2],rr_sx,ones(2,length(rr_sx))'*(-max(max(real(V0)))));
                pcolor(zz_c,[c-Drr/2,c],ones(length(zz_c),2)'*(-max(max(real(V0)))));
                hold off;
%                 pause(0.5)
                frame = getframe;
                writeVideo(writerObj,frame);
            end
            close(writerObj);

            % Movie Er

%             zz_sx=[-L-Dzz/2:Dzz:Dzz/2];
%             zz_dx=[L-Dzz/2:Dzz:2*L];
%             zz_c=[0+Dzz/4:Dzz/2:L];
%             rr_sx=[b+Drr/2:Drr:c];
% 
%             omega=2*pi*f;
%             figure()
%             T=1/f;
%             V0=(PEr);
%             crange=[-max(max(real(V0))) max(max(real(V0)))];
%             colormap(jet(100));
%             for tt=-T/2:T/10:T/2
%                 V=real(V0.*exp(1i*omega*tt));
%                 pcolor(zz_mesh,rr_mesh,V); shading interp; hold on; contour(zz_mesh,rr_mesh,V); 
%                 xlabel('Length [m]'); ylabel('Radius [m]'); title('Field E_{r} [V/m]');caxis(crange);colorbar;
%                 pcolor(zz_sx,[b-Drr/4,b+Drr/4],ones(length(zz_sx),2)'*(-max(max(real(V0)))));
%                 pcolor(zz_dx,[b-Drr/4,b+Drr/4],ones(length(zz_dx),2)'*(-max(max(real(V0)))));
%                 pcolor([0+Dzz/4,0+Dzz/2],rr_sx,ones(2,length(rr_sx))'*(-max(max(real(V0)))));
%                 pcolor([L-Dzz/4,L-Dzz/2],rr_sx,ones(2,length(rr_sx))'*(-max(max(real(V0)))));
%                 pcolor(zz_c,[c-Drr/2,c],ones(length(zz_c),2)'*(-max(max(real(V0)))));
%                 hold off;
%                 pause(0.3)
%             end
% 
            % Movie Hphi

%             zz_sx=[-L-Dzz/2:Dzz:Dzz/2];
%             zz_dx=[L-Dzz/2:Dzz:2*L];
%             zz_c=[0+Dzz/4:Dzz/2:L];
%             rr_sx=[b+Drr/2:Drr:c];
% 
%             omega=2*pi*f;
%             figure()
%             T=1/f;
%             V0=(PHphi);
%             crange=[-max(max(real(V0))) max(max(real(V0)))];
%             colormap(jet(100));
%             for tt=-T/2:T/10:T/2
%                 V=real(V0.*exp(1i*omega*tt));
%                 pcolor(zz_mesh,rr_mesh,V); shading interp; hold on; contour(zz_mesh,rr_mesh,V); 
%                 xlabel('Length [m]'); ylabel('Radius [m]'); title('Field H_{phi} [A/m]');caxis(crange);colorbar;
%                 pcolor(zz_sx,[b-Drr/4,b+Drr/4],ones(length(zz_sx),2)'*(-max(max(real(V0)))));
%                 pcolor(zz_dx,[b-Drr/4,b+Drr/4],ones(length(zz_dx),2)'*(-max(max(real(V0)))));
%                 pcolor([0+Dzz/4,0+Dzz/2],rr_sx,ones(2,length(rr_sx))'*(-max(max(real(V0)))));
%                 pcolor([L-Dzz/4,L-Dzz/2],rr_sx,ones(2,length(rr_sx))'*(-max(max(real(V0)))));
%                 pcolor(zz_c,[c-Drr/2,c],ones(length(zz_c),2)'*(-max(max(real(V0)))));
%                 hold off;
%                 pause(0.3)
%             end
    
end
