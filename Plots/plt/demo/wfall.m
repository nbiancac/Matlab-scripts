% wfall.m ---------------------------------------------------------------------
%
% * Shows one way to create a waterfall plot with hidden line removal.
% * Note how plt is called with a single trace color and no TraceID box.
% * Extensive use of the slider pseudo object to control the plotted data.
% * Figure user data is used to pass the handle structure (S) to the callback.
% * the 'Linesmoothing' option is selected (which surprisingly speeds
%   up the display on many systems).
% * A pseudo popup in "super-button" mode is used to start and stop the display.
% * The number of display updates per second is calculated every second with the
%   results shown in a large font below the plot.
% * Demonstrates the use of the plt 'closeReq' parameter.

% ----- Author: ----- Paul Mennen
% ----- Email:  ----- paul@mennen.org

function wfall(go)                      % wfall(0) or with no arguments, starts stopped
  if ~nargin go=0; else go = ~~go; end; % wfall(1) or wfall('run') starts running
  S.n  = 30;                            % number of traces
  S.sz = 256;                           % length of each trace
  z = repmat(zeros(S.sz,1),1,S.n);      % define S.n traces
  S.tr = plt(z,z,'TraceC',[0 1 0],'FigName','wfall','FigBKc',[1 1 1]/7,...
        'TraceID','','LabelY','','xy',[.173 .1 .815 .88],'Ylim',[0 10],...
        'Options','Ticks Linesmooth -X -Y','Pos',[9 45 900 600],...
        'CloseReq',{@clos date});
  p = .91:-.09:.2; p = num2cell([.01+0*p; p]',2); % slider xy positions
  S.ax   = gca;
  S.dly  = plt('slider',p{1},[ 0 0 999],'delay (ms)','',2);
  S.frq1 = plt('slider',p{2},[.1  0 .5],'Freq start');
  S.frq2 = plt('slider',p{3},[.45 0 .5],'Freq end');
  S.delt = plt('slider',p{4},[7  0 100],'Freq delta * 1024');
  S.dist = plt('slider',p{5},[2.7 .1 4],'Distortion (10^-x)');
  S.dx   = plt('slider',p{6},[3 0 20 0 S.sz-1],'Delta X','',2);
  S.dy   = plt('slider',p{7},[1 0  4 0  999],'Delta Y');
  S.go   = plt('pop',[.04 -.868 .1 .1],{' start ' ' stop '},'callbk',@start);
  S.ups  = text(.13,-.075,'','units','norm','color',[1 1 .5],'fontsize',14);
  set(S.tr,'z',zeros(1,S.sz)+NaN);      % initialize z data
  set(S.tr(1),'user',0);                % initialize count (# of plot updates)
  set(gcf,'user',S);
  plt('pop',S.go,'index',-go-1);        % start running if go is 1
% end function wfall;

function start(a,b) % start/stop button callback ----------------------------
  S = get(gcf,'user');
  n = S.n;  sz = S.sz;  sz2 = sz*2;              % fft length
  x = 1:sz;  y0 = 0:sz2-1;                       % for computing new trace data
  w = .5 -.5*cos((pi/sz)*(1:sz2));               % hanning window
  dxy = -1;  m = 0;  b = 0;  tic;
  f = plt('slider',S.frq1,'get');
  while ishandle(S.go) & plt('pop',S.go,'get')>1 % loop until stop button clicked
    set(S.tr(1),'user',get(S.tr(1),'user')+1);   % increment count (# of updates)
    dx     = plt('slider',S.dx,  'get');         % horizontal steps
    dy     = plt('slider',S.dy,  'get');         % vertical steps
    fstart = plt('slider',S.frq1,'get');         % advance sine wave frequency
    fend   = plt('slider',S.frq2,'get');
    fdelta = plt('slider',S.delt,'get')/1024;
    dlySec = plt('slider',S.dly, 'get')/1000;    % delay in seconds
    clip = 1 - 10^-plt('slider',S.dist,'get');   % add distortion by clipping
    d = dx + dy;
    if d ~= dxy;                                 % was a Delta X/Y slider moved?
      dxy = d;                                   % yes, modify the axis limits
      set(S.ax,'xlim',[0 sz+n*dx],'ylim',[-.3 11+n*dy]);
      for k=1:n set(S.tr(k),'x',x+k*dx); end;    % set new x data
    end;
    yc = zeros(1,sz);                     % hidden line removal comparison vector
    for k=1:n-1                           % draw all but the last trace
       z = get(S.tr(k+1),'z');            % get the data from the trace above
       ym = max(yc,z);                    % stay above the trace in front
       yc = ym(dx+1:end) - dy;            % slide comparison vector to the left
       yc = [yc zeros(1,sz-length(yc))];  % and zero pad to size sz
       set(S.tr(k),'y',ym+k*dy,'z',z);    % save unclipped data in z axis
    end;
    y = fft(min(max(sin(2*pi*f*y0),-clip),clip) .* w);  % compute new trace
    z = 5 + max(-5,log(abs(y(x))));        % new data comes in at the back
    set(S.tr(n),'y',max(yc,z)+n*dy,'z',z); % draw last trace
    if fend>fstart  f = f + fdelta;  if f>fend f=fstart; end;
    else            f = f - fdelta;  if f<fend f=fstart; end;
    end;
    t = toc;  m = m + 1;                   % count display updates
    if t>1 set(S.ups,'string',sprintf('%.2f updates/sec',m/t));
           m = 0; tic;
    end;
    pause(dlySec); drawnow;                    % delay the selected amount
  end;  % end while ishandle
% end function start

function clos(sd)        % close request function
  S = get(gcf,'user');   % displays start date and total # of display updates
  prin(1,'Start date:      %s\nDisplay updates: %d\n',sd,get(S.tr(1),'user'));
% end function clos