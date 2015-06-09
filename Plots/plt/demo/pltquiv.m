% pltquiv.m ---------------------------------------------------------
%
% * The Pquiv.m function appears twice in the plt argument list to plot
%   two vector fields both with the base location specified by f and
%   lengths specified by v1 and v2 respectively. The first of these
%   Pquiv calls is somewhat similar to the  MatLab command
%   quiver(real(f),imag(f),real(v1),imag(v1));
% * Using the xy parameter to make room for long Trace ID names
% * Using tex commands (e.g. \uparrow) inside Trace ID names
% * Reassigning menu box items. In this example, the 'LinX' button is
%   replaced by a 'Filter' button. Its button down function (which is
%   executed when you click on 'Filter') searches for the 4th trace
%   (findobj) and swaps the contents of its user data and y-axis data.
% * Adding text items to the figure. Note that the text position is
%   specified using x and y axes coordinates
% * Using NaNs (not a number) to blank out portions of a trace
% * Using the TraceID callback function (TIDcback) to perform an action
%   when you click on a trace ID. For example, when you click on the
%   last trace ID (humps+rand) this will appear in the command window:
%   "A trace named humps+rand and color [1 0 0] was toggled".
%   (This TraceID callback was contrived to use all the substitutions,
%   and is not particularly useful.)

% ----- Author: ----- Paul Mennen
% ----- Email:  ----- paul@mennen.org

x  = (0:.08:5)';
x4 = (0:.01:5)';
t = x/5;
y = humps(t)/20;
f = complex(x,y);

v1 = complex(exp(-2*t).*sin(20*t), t .* cos(15*(1-t).^3));
v2 = exp(-1.4*t) .* exp(30i * t.^.5)/2;
y4  = 6 + rand(size(x4)) - humps(flipud(x4/5))/20;
p  = [1 .173 .107 .813 .874;  % plotting axis position
     -1 .006 .785 .131 .196]; % TraceID box position
h = plt(f,Pquiv(f,v1),Pquiv(f,v2),x4,y4,'xy',p,...
    'Xlim',[-.2 5.2],'Ylim',[0 6.6],'TraceID',...
    {'humps \div 20','velocity1 \uparrow','velocity2 \uparrow','humps+rand'},...
    'FigBKc',[0 0 .3],'Options','Slider-Y-M','FigName','pltquiv','TIDcback',...
    ['prin(1,"A trace named  %s and color [{%3w! }]] was toggled\n",' ...
      'get(@TID,"string"),get(@LINE,"color"));']);
set(h(1),'LineWidth',2);
y4 = filter([1 1 1]/3,1,y4); y4 = y4([3 3:end end]); % smoothed y4
set(h(4),'tag','h4','user',y4); % save smoothed y4 in trace user data
bfn = 'h=findobj("tag","h4"); set(h,"y",get(h,"user"),"user",get(h,"y"));';
LinXtag = findobj(gcf,'string','LinX');
set(LinXtag,'string','Filter','ButtonDownFcn',strrep(bfn,'"',''''));
set([text(2.8,5.8,{'Click on ''Filter''' 'in the menu box'});
     text(3.4,3.6,'NaN induced gap','fontangle','italic')],...
     'fontsize',14,'color','yellow');
x4(380:400) = NaN;  set(h(4),'x',x4);  % create a gap in trace 4

% Note: Identical results are achieved when:
% Pquiv(f,v1),Pquiv(f,v2)
%   is replaced by:
% Pquiv([f f],[v1 v2])
