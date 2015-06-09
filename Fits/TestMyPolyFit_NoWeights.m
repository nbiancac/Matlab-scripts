% Test the new version of the regression routine

DataX=-5:1:5;
DataY=DataX+0.5*randn(1,length(DataX));
uncY=0.5*randn(1,length(DataX));

[coeff,unc_coeff]=mypolyfit(DataX,DataY,1);
xx=-5.5:0.1:5.5;
yy=coeff(1).*xx+coeff(2);

% Best estimate
[estimate,uncertainty]=mypolyfitEstimate(DataX,DataY,1);

%Printing to the screen
sprintf('Pendenza:%f (%f) Intercetta: %f (%f)',...
    coeff(1),unc_coeff(1),coeff(2),unc_coeff(2))

figure;
axes1 = axes('FontSize',14,'LineWidth',2,'Parent',figure1);
grid(axes1,'on');
box(axes1,'on');
hold(axes1,'all');

plot(DataX,DataY,'ro','Linewidth',2);
plot(xx,yy,'g','Linewidth',2);
errorbar(DataX,estimate,uncertainty,'kx','Linewidth',2);
