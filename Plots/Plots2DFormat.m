%% PAPER
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperType', 'A4');

%% FONTS
set(0,'defaultlinelinewidth',1)
set(0,'defaultaxesfontname','timesnewroman');
set(0,'defaulttextfontname','arial');
set(0,'DefaultAxesFontSize', 14) 
set(0,'DefaultTextFontSize', 14) 
set(0,'defaulttextinterpreter','TeX');

%% BAR
t = colorbar('peer',gca);
set(get(t,'ylabel'),'String', 'v/c');
set(get(t,'ylabel'),'Rotation',-90);
ylabel_bar_pos=get(t,'ylabel');
set(ylabel_bar_pos,'Position',get(ylabel_bar_pos,'position')+[1.2 0.05 0]);

%% LEGEND
% top right inner
set(h_legend,'FontSize',12) 
set(h_legend,'FontName','timesnewroman') 
set(h_legend,'Position',[0.65    0.8    0.25    0.1048]); 
% bottom right inner
set(h_legend,'FontSize',12) 
set(h_legend,'FontName','timesnewroman') 
set(h_legend,'Position',[0.5307    0.1410    0.3786    0.2675]); 
% top left inner
set(h_legend,'FontSize',12) 
set(h_legend,'FontName','timesnewroman') 
set(h_legend,'Position',[0.1006    0.8284    0.3580    0.0848]); 

%% DYNAMIC LEGEND
plot(x,y,'DisplayName','pippo')
legend('-DynamicLegend')