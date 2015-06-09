% pub2.m ---------------------------------------------------------------------
%
% In this example, a plt figure is created in its usual data exploration mode
% showing 6 traces of randomly generated data. Each trace contains over 50
% thousand data points, although the display is zoomed to show only a small
% portion of the data. The idea is to use the cursor controls to pan and/or
% zoom the display to some area of interest and then press the "pub" button to
% generate a figure optimized for publication.
%
% What makes this more interesting is that when you pan to a new section of
% the data and again press the "pub" button, the publication figure is redrawn
% using subplots to show both selected portions. In a like manner, successive
% presses of the pub button further subdivide the plotting area with each new
% data range appearing above the previous ones.
%
% The x axis of the data exploration window is plotted in units of days past
% a time reference (1-Jan-2013 in this example), but custom date ticks are
% used on the x axis of the publication plot. To reduce clutter, only the day
% and month are shown for all vertical grid lines except the last one (which
% includes day, month, & year).
%
% The 'xView' option is used to enable the xVview slider which is particalarly
% useful in situations like this where you are viewing only a small portion
% of a long data record. (The xView slider appears above the primary figure.)
%
% The TraceID box is typically placed to the left of the plot, although for
% the publication figure in this demo the TraceID box is placed right on top
% of the plot (more like a legend). This means that sometimes the TraceID box
% will obscure some of the data, but note that you can easily use the mouse
% to drag the legend around to a spot that does not interfere with the plot.
%
% This example may seem somewhat contrived - and indeed it was conceived
% mostly to demonstrate as many unusual plt parameters and programming
% techniques as possible.

% ----- Author: ----- Paul Mennen
% ----- Email:  ----- paul@mennen.org

function pub2(In1)
  ylbl = {'mm/hour' 'liters/sec'};      % axis labels
  yl = {[0 80] [0 16]};                  % y limits
  pos = [10 55 900 400]; fname = 'pub2'; % figure position/name
  if ~nargin
    tid = [char('a'+26*rand(6,3)) char('0'+10*rand(6,4))]; % random trace labels
    spd = 144;                 % samples per day (delta T = 10 minutes)
    nd = 365;                  % number of days in the data set
    ns = nd*spd;               % number of samples in the data set
    t0 = '1-Jan-2013';         % time reference
    x = (0:ns-1)/spd;          % time index (days since time reference)
    y = cumsum(rand(1,ns)-.5); % generate random data ------------------------
    y2 = max(abs(y - filter(ones(1,1000)/1000,1,y)),7) - 7;
    [y2m yi] = max(y2);               y2 = (55/y2m) * y2;
    y1 = [y2(500:end) y2(1:499)]/10;  m = 10*find(y1(1:10:end)>0);
    y3 = filter([.15 .15],1,y2);      y4 = filter(zeros(1,5)+.35,1,y3);
    y5 = filter(zeros(1,7)+.2,1,y4);  y6 = filter(zeros(1,20)+.1,1,y5);
    plt(Pvbar(x(m),0,y1(m)),x,[y2;y3;y4;y5;y6],... % create main figure -------
       'Right',1,'Ylim',yl,'Ylabel',ylbl,'-Ydir','rev','Figname',fname,...
       'Linewidth',{2 1 1 1 1 1},'Options','xView slider',...
       'TraceC',[.2 .4 1; 0 .8 0; 1 1 0; 0 1 1; 1 0 1; 1 0 0],...
       'Xlabel',['Days since ' t0],'Pos',pos,'xy',[.105 .130 .855 .840],...
       'Xlim',x(yi)+[-5 8],'Xstring','datestr(@XU+@XVAL)','TraceID',tid);
    set(findobj(gcf,'tag','xstr'),'User',datenum(t0)); % save time reference
    uicontrol('string','pub','pos',[15 250 40 20],...  % create pub button
              'callback','pub2(0);');
    return;
  end;

  f = findobj('name',fname);          % generate pub figure ----------------------
  if length(f)==1
       xl={}; xt={}; n=1; nt=6;       % Here if pub button pushed for the 1st time
       arg1 = 'TraceMK';  arg2=.5;    % Repeat a parameter occuring earlier in the
                                      % parameter list (since we don't need Subplot)
       arg3 = arg1;  arg4=arg2;
  else ax = getappdata(f(2),'axis');  % Here if pub figure already exists
       n = length(ax)+1; nt = 6*n;    % number of axes/traces for this figure
       xl = get(ax,'xlim');           % get x limits for previous axes
       xt = get(ax,'xticklabel');     % get x ticks for previous axes
       if n==2 xl = {xl}; xt = {xt}; end;
       set(f(2),'tag',''); delete(f(2)); % delete previous pub figure
       arg1 = 'Subplot';  arg2 = ones(1,n) * 100/n;
       arg3 = 'Subtrace'; arg4 = repmat(6,1,n);
   end;

  s = datestr(get(findobj(gcf,'tag','xstr'),'User') + get(gca,'xtick'),1);
  c = cellstr(s(:,1:6));  c{end} = [c{end} s(end,7:11)]; % create date ticks
  xt = [xt; {c}];                           % add x ticks for newest axis
  xl = [xl; {get(gca,'xlim')}];             % add x limits for newest axis
  c = repmat([0 0 1; 0 .7 0; .7 .7 0; 0 .7 .7; .7 0 .7; 1 0 0],n,1); % trace colors
  lh = getappdata(gcf,'Lhandles');
  x = get(lh,'x'); y = get(lh,'y');                  % get x/y data from main figure
  tbox = findobj(gcf,'user','TraceID');              % get TracID's from TraceID box
  tid = get(flipud(findobj(tbox,'type','text')),'string');
  set(0,'Units','pixels'); sz = get(0,'screensize'); % get screen size
  pos = pos + [0 435 0 min(180*(n-1),sz(4)-960)];    % don't exceed screen height
  p = [ 0 .051 .052 .905 .922;                       % position: plotting area
       -1 .780 .080 .142 100/pos(4)];                % position: TraceID box
  l = plt(1,ones(1,nt),'Ylabel',ylbl{1},'Figname',fname,'Linewidth',2,...
     'GRIDc',[.85 .85 -.85],'Pos',pos,'Link',gcf,'ColorDef',0,...
     'TraceID',tid,'Options','-All Nocursor','xy',p,'TraceMK',.5,...
     'Xlabel','','TIDc',[1 .95 .9],'TraceC',c,arg1,arg2,arg3,arg4);
  for k = 1:nt                   % set x/y data for all lines of pub figure
    m = 1+mod((k-1),6);
    set(l(k),'x',x{m},'y',y{m});
  end;
  ax = getappdata(gcf,'axis');   % get list of axes
  hs = .05;                      % vertical separation between axes
  dh = .97/n;                    % vertical space allocated for each axis
  k = 1;  m = 1;  yp = hs;
  for a=ax                                         % iterate over axes
    p=get(a,'pos'); p(2)=yp; p(4)=dh-hs; yp=yp+dh; % compute axis position
    set(a,'pos',p,'xlim',xl{k},'ylim',yl{1},...;   % set xy position/limits
          'xticklabel',xt{k})
    plt('grid',a);                                 % update grid lines
    v = axes('pos',p,'color','none','yaxisl',...   % create right hand axis
        'right','ydir','rev','ycolor','blue',...
        'xlim',get(a,'xlim'),'ylim',[0 16],'xtick',[]);
    set(l(m),'parent',v);                          % move y1 to right axis
    set(get(a,'ylabel'),'string',ylbl{1});         % set y axis labels
    set(get(v,'ylabel'),'string',ylbl{2});
    k = k+1;  m = m+6;
  end;
  set([findobj(gcf,'user','grid'); l; ax'],'buttondown',''); % disable pan/zoom
  plt misc tidtop;  % traceID box must be on top so it can be moved via the mouse
