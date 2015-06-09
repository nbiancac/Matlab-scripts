% pltsqALT.m -------------------------------------------------------------
%
% pltsqALT approximates a square wave by adding up the first five odd
% harmonics of a sine wave. To make the GUI interface more challenging
% and interesting, the amplitude is continually varied (periodically
% between plus and minus 1).
% * Demonstrates how you can add additional GUI controls to the plt
%   window. Typically this is something you will want to do when
%   creating an application based on plt.
% * 10 uicontrols are added to the figure (2 buttons & 4 popups with labels).
% * A text object is added below the plot which is used to display
%   "updates/seccond" - a good measure of computational & graphics performance.
%   The color of this text object is toggled every time it is refreshed so that
%   you can tell the speed is being recomputed even if the result is the same.
% * The 'xy' argument makes room for the uicontrols added above
%   the plot area and to make room for the wider TraceIDs
% * The position coordinates for the uicontrols are listed in two arrays
%   (the 1st set for the popups and buttons and the 2nd set for the text
%   styles). This makes it easy to generate these position coordinates
%   (as well as for the 'xy' argument) using plt's move function. To see
%   how this is done, look at the gui1.m and gui2.m examples.
% * Normalized units are used here for the uicontrols, although the plt
%   move function also handles pixel units which is better if you don't
%   want the objects to change size when the figure window is resized.
% * The cursor callback ('moveCB') and the plt('rename') call are used
%   to provide simultaneous cursor readouts for all 5 traces in the
%   TraceID box. (Click the stop button first to make it easier to
%   cursor, and then click anywhere in the plot area to see all 5 cursor
%   values update.)
% * The 'Options' argument is used to turn off grid lines (initially)
%   and to remove the y-axis Log selector from the menu box.
% * You can use the Erasemode popup, you can explore the effect of the
%   erasemode property on drawing speed.

% This is an alternate (older) version of pltsq.m which performs exactly
% the same function but uses uicontrols instead of the popup pseudo objects
% used in pltsq.m.  This allows you to compare and contrast the relative
% merits of the two types of objects. Another minor difference is that the
% S structure containing the handles is passed to the callbacks using its
% argument list as opposed to passing it through the figure's user data
% as done in pltsq.m. Both methods have their merits.
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
% more powerful pseudo objects. (This alternate version is not run by demoplt.m).

% ----- Author: ----- Paul Mennen
% ----- Email:  ----- paul@mennen.org

function pltsqALT()
  S.tr = plt(0,[0 0 0 0 0],'FigName','pltsqALT','Options','Ticks/-Xlog-Ylog',...
        'xy',[0 .184 .093 .797 .89; -1 .01 .81 .13 .18],...
        'Xlim',[0 4*pi],'Ylim',[-1.05,1.05],'LabelX','Cycles','LabelY','' );
  S.go  = uicontrol('string','Start');
  S.stp = uicontrol('string','Stop','User',1);
  S.spd = uicontrol('string',prin('{%d!row}',2.^(0:10)));
  S.cyc = uicontrol('string',prin('{%d!row}',2.^(0:5)));
  S.pts = uicontrol('string',prin('{%d!row}',25*2.^(0:6)));
  S.era = uicontrol('string',{'normal' 'background' 'xor' 'none'});
  S.ups = text(.13,-.07,'','units','norm','fontsize',13);

  ui = [S.go S.stp S.spd S.cyc S.pts S.era]; % 2 buttons & 4 popups
  p = [.010 .740 .060 .040;  % position: Start button
       .080 .740 .050 .040;  % position: Stop button
       .010 .640 .120 .030;  % position: Speed popup
       .010 .530 .120 .040;  % position: Cycles popup
       .010 .440 .120 .030;  % position: Points/cycle popup
       .010 .340 .120 .030]; % position: EraseMode popup
  S.cid = getappdata(gcf,'cid'); % get cursor ID
  set(ui,'Units','Norm',{'Position'},num2cell(p,2),'CallBack',{@uiCB,S});
  set(ui(3:6),'Style','Popup','BackGround',[0 1 1],{'Value'},{5;1;3;3});
  set(S.stp,'CallBack','set(gcbo,''User'',0);');
  % create a label to identify each popup menu
  p = [.010 .670 .120 .030;  % position: Speed popup lablel
       .010 .570 .120 .030;  % position: Cycles popup label
       .010 .470 .120 .030;  % position: Points/cycle popup label
       .010 .370 .120 .030]; % position: EraseMode popup label
  set([uicontrol; uicontrol; uicontrol; uicontrol],'Style','Text',...
      {'String'},{'Speed';'Cycles';'Points/cycle';'EraseMode'},...
      'Units','Norm',{'Position'},num2cell(p,2));
  plt('cursor',S.cid,'set','moveCB',{@curCB,S});  % set cursor callback
  curCB(S);  uiCB(0,0,S);                         % start it moving
% end function pltsqALT

function curCB(S)                                  % cursor callback
   [xy m] = plt('cursor',S.cid,'get','position');  % get cursor index (m)
   for k=1:5 v = get(S.tr(k),'y'); y(k)=v(m); end; % get y value for each trace
   t = 'Fund %5v ~, + 3rd %5v ~, + 5th %5v ~, + 7th %5v ~, + 9th %5v';
   plt('rename',prin(t,y));
%end function curCB

function uiCB(h,arg2,S)               % uicontrol callbacks
  Ncyc  = 2^get(S.cyc,'Value')/2;     % Number of cycles
  Pts   = 12.5*2^get(S.pts,'Value');  % Points per cycle
  Speed = 2^get(S.spd,'Value')/2;     % Amplitude step size
  x = [0: 1/Pts : Ncyc];              % initialize x axis values
  plt('cursor',S.cid,'set','xlim',[0 Ncyc]);
  y = zeros(5,length(x));  v = y(1,:);
  m=1;  for k=1:5 v=v+sin(2*pi*m*x)/m;  m=m+2; y(k,:)=v; end;
  emode = get(S.era,'string');
  set(S.tr,'x',x,'y',0*x,'EraseMode',emode{get(S.era,'Value')});
  set(S.stp,'User',1);  b=0; m=0; n=0; tic;
  while ishandle(S.stp) & get(S.stp,'User') % loop until stop is clicked
    m=m+1; n=n+1; a=sin(Speed*pi*m/32768);  % update amplitude based on Speed
    t = toc;                                % refresh updates/sec text on the next
    if t>1 & ~mod(Speed*n,32768)            % half cycle after every second
      set(S.ups,'color',[b 1 .5],'string',sprintf('%d updates/sec',round(n/t)));
      n=0;  b=1-b;  tic;
    end;
    for k=1:5 set(S.tr(k),'y',a*y(k,:)); end;  % update trace values
    drawnow;
  end;
%end function uiCB
