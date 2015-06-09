% pltsq.m -------------------------------------------------------------
%
% pltsq approximates a square wave by adding up the first five odd
% harmonics of a sine wave. The plot displays successive sums of these
% harmonics which approximates a square wave more closely as more harmonics
% are added together. The key point however (and the reason this demo was
% created) is that the amplitudes of these sine waves and sums are continually
% varied (periodically between plus and minus one) to produce a "real-time"
% moving display. plt is well suited to creating real-time displays, but there
% are a few concepts to learn and this demo is an excellent starting point.
% * Demonstrates how you can add GUI controls to the plt window - typically
%   something you will need to do when creating plt based applications.
% * 5 pseudo popup controls are added to the figure to the left of the plot
%   including one "super-button" to start and stop the plotting.
% * A text object appears below the plot which displays "updates/second" - 
%   a good measure of computational & graphics performance.
%   The color of this text object is toggled every time it is refreshed so that
%   you can tell the speed is being recomputed even if the result is the same.
% * The 'xy' argument is used to make room for pseudo popups as well as for
%   the wider than usual TraceIDs.
% * In the code, the position coordinates for the five popups are grouped in a
%   single array to make it easy to updated these coordinates with the plt move
%   function. For details on how this is done, refer to the gui1 & gui2 examples.
% * Normalized units are used here for the uicontrols. The "plt move function
%   also handles pixel units which is better if you don't
%   want the objects to change size when the figure window is resized.
% * The cursor callback parameter ('moveCB') and the plt('rename') call are
%   used to provide simultaneous cursor readouts for all 5 traces in the
%   TraceID box. (Although you can cursor the traces while the displays are
%   updating, to make it easier to see what is going on, click the stop button
%   first, then click anywhere in the plot area to see all 5 cursor values update.)
% * The 'Options' argument is used to turn off grid lines and to remove the
%   x and y-axis Log selectors from the menu box.
% * You can use the Erasemode popup to explore the effect of the erasemode line
%   property on drawing speed. (The erasemode property is no longer supported in
%   Matlab version R2014b or later, so pltsq.m checks the Matlab version and
%   disables the popup appropriately.) You can also effect the drawing speed
%   by varying the number of points per plot from a low of 25 points to a
%   high of 51200 points (32 cylces times 1600 points per cycle). 

% The demo folder also includes an alternate (older) version of this program
% (called pltsqALT.m) that performs exactly the same function except that
% uicontrols are used in place of the popup pseudo objects. This allows you
% to compare and contrast the relative merits of the two types of objects.
% Also note that the alternate version passes the S structure to the callbacks
% using a different method (thru the parameter list). Both methods have merit.
% The control structure of pltsqALT is a bit unusual in that if you change
% parameters with the uicontrols while the display is left running, a new while
% loop is started, leaving all the previous while loops still running. This
% doesn't cause any problems however since only the newest while loop has a
% chance to do anything and all the older while loops are properly terminated
% when you hit the stop button or exit the program. The newer version (pltsq.m)
% avoids this infinite nesting by separating the start/stop callback from the
% other control callbacks. This slightly complicates the communication of
% variables between the parameter and display loops and slows the update rate
% down by about 6%, but still it may be the morally correct approach. Despite
% this additional complexity, the newer version (with 51 lines of code) is still
% shorter than the older version (58 lines of code) because of the use of the
% more powerful pseudo objects. (The alternate version is not run by demoplt.m).

% ----- Author: ----- Paul Mennen
% ----- Email:  ----- paul@mennen.org

function pltsq(go)                      % pltsq(0) or with no arguments, starts stopped
  if ~nargin go=0; else go = ~~go; end; % pltsq(1) or pltsq('run') starts running
  S.tr = plt(0,[0 0 0 0 0],'FigName','pltsq','Options','Ticks/-Xlog-Ylog',...
    'xy',[0 .184 .093 .797 .89; -1 .01 .81 .13 .18],'moveCB',@curCB,...
    'Xlim',[0 4*pi],'Ylim',[-1.05,1.05],'LabelX','Cycles','LabelY','' );
  Tgo = {' start ' ' stop '};  Tera = {'normal' 'background' 'xor' 'none'};
  p = [.05 .67 .10 .1;                  % position: start/stop super button
       .01 .26 .10 .4;                  % position: speed popup
       .01 .25 .10 .3;                  % position: cycles popup
       .01 .15 .10 .3;                  % position: # of points popup
       .01 .15 .12 .2];                 % position: eraseMode popup
  er=version;  er = str2num(er(1:3)) < 8.4; % disable erase mode for R2014b
  S.go  = plt('pop',p(1,:),Tgo,        'callbk',@gCB,'swap',0);
  S.spd = plt('pop',p(2,:),2.^(0:10),  'callbk',@pCB,'labely','Speed','index',5);
  S.cyc = plt('pop',p(3,:),2.^(0:5),   'callbk',@pCB,'labely','Cycles');
  S.pts = plt('pop',p(4,:),25*2.^(0:6),'callbk',@pCB,'labely','Points/cycle','index',3);
  S.era = plt('pop',p(5,:),Tera,       'callbk',@pCB,'labely','EraseMode',...
              'index',1+2*er,'enable',er);
  S.cid = getappdata(gcf,'cid');                  % get cursor ID
  S.ups = text(.13,-.07,'','units','norm','fontsize',13);
  set(gcf,'user',S); curCB; pCB;                  % initialize plot
  plt('pop',S.go,'index',-go-1);                  % start running if go is 1
% end function pltsq

function gCB(a,b)                                 % start/stop super-button callback
  S = get(gcf,'user');  b=0; m=0; n=0; tic;
  while ishandle(S.go) & plt('pop',S.go,'get')>1
    if get(S.ups,'user') m=0; n=0; set(S.ups,'user',0); S = get(gcf,'user'); tic; end;
    m=m+1; n=n+1;  a = sin(S.Speed*pi*m/32768);   % update amplitude based on Speed
    t = toc;                                      % refresh updates/sec text on the next
    if t>1 & ~mod(S.Speed*n,32768)                % half cycle after every second
      set(S.ups,'color',[b 1 .5],'string',sprintf('%d updates/sec',round(n/t)));
      n=0;  b=1-b;  tic;
    end;
    for k=1:5 set(S.tr(k),'y',a*S.y(k,:)); end;   % update trace values
    drawnow;
  end;
% end function gCB

function curCB(a,b)                               % cursor callback
  S = get(gcf,'user');  if isempty(S) return; end;
  [xy m] = plt('cursor',S.cid,'get','pos');       % get cursor index (m)
  for k=1:5 v = get(S.tr(k),'y'); y(k)=v(m); end; % get y value for each trace
  t = 'Fund %5v ~, + 3rd %5v ~, + 5th %5v ~, + 7th %5v ~, + 9th %5v';
  plt('rename',prin(t,y));
% end function curCB

function pCB(a,b)                                 % popup callbacks
  S = get(gcf,'user');
  Ncyc    = str2num(get(S.cyc,'string'));         % Number of cycles
  Pts     = str2num(get(S.pts,'string'));         % Points per cycle
  S.Speed = str2num(get(S.spd,'string'));         % Amplitude step size
  set(S.ups,'user',1);                            % indicate parameters were updated
  x = [0: 1/Pts : Ncyc];                          % initialize x axis values
  plt('cursor',S.cid,'set','xlim',[0 Ncyc]);
  S.y = zeros(5,length(x));  v = S.y(1,:);
  m=1;  for k=1:5 v=v+sin(2*pi*m*x)/m;  m=m+2; S.y(k,:)=v; end;
  set(S.tr,'x',x,'y',0*x);
  if plt('pop',S.era,'get','enable') set(S.tr,'Erase',get(S.era,'String')); end;
  set(gcf,'user',S);
% end function pCB
