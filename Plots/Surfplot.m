   
    figure(); pcolor(zz_mesh,rr_mesh,abs(PEr)); shading interp; hold on; contour(zz_mesh,rr_mesh,abs(PEr)); hold off;
    xlabel('Length [m]'); ylabel('Radius [m]'); title('Field E_r [V/m]');
    figure(); pcolor(zz_mesh,rr_mesh,abs(PEz)); shading interp; hold on; contour(zz_mesh,rr_mesh,abs(PEz)); hold off;
    xlabel('Length [m]'); ylabel('Radius [m]'); title('Field E_z [V/m]');
    figure(); pcolor(zz_mesh,rr_mesh,abs(PHphi)); shading interp; hold on; contour(zz_mesh,rr_mesh,abs(PHphi)); hold off;
    xlabel('Length [m]'); ylabel('Radius [m]'); title('Field H_{phi} [A/m]');