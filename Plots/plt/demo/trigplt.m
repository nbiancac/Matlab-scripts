% trigplt.m ------------------------------------------------------------
%
% This example demonstrates:
% - The use of the slider pseudo object
% - showing the line characteristics in the TraceID using the TraceMK parameter
% - setting the cursor callback with the moveCB parameter
% - setting axis, TraceID box, and MenuBox positons using the 'xy' parameter
% - setting trace characteristics with the Linewidth, Styles, and Markers parameters
% - setting an initial cursor position
% - enabling the multiCursor mode
% - how to modify the colors and fonts of the Trace IDs.
%
% Text objects are used to display help information at the top of the
% plot window. This help text dissappears when any parameter is changed
% but can be re-enabled by clicking on the help button.

% The clipboard button captures the figure as a bitmap into the clipboard

% ----- Author: ----- Paul Mennen
% ----- Email:  ----- paul@mennen.org

function trigplt
  p = [ 0 .143 .093 .790 .733;  % position: plot axis
       -1 .009 .105 .100 .217;  % position: TraceID box
       -2 .017 .620 .058 .210]; % position: MenuBox
  S.tid = {'sin' 'cos' 'tan' 'csc' 'sec' 'cot'};
  S.tr = plt(0,zeros(1,6),'FigName','trigplt: y = A sin(Bx + C) + D',...
          'Right',1:2,'Pos',[10 45 780 530],'xy',p,'FigBKc',[0 .25 .3],...
          'Xlim',[0 2*pi],'Ylim',{[-8 8] [-1.5 1.5]},...
          'moveCB',@sliderCB,'DisTrace',[0 0 0 0 0 1],'TraceID',S.tid,...
          'Linewidth',{2 2 1 1 1 1},'Styles','----:-','Markers','nnnnns',...
          'Options','multiCur -X -Y','-Ycolor',[1 .3 .3],'TraceMK',.5,...
          'LabelX','radians','LabelY',{'tan / csc / sec / cot'; 'sin / cos'});
  S.tx = text(.35,1.04,'','fontsize',14,'color','yellow','units','norm');
  S.cid = get(gca,'UserData');                             % save cursor ID
  line([-99 99],[0 0],'color',[1 .3 .3],'linestyle','--'); % make x-axis more obvious
  p = [.04 .955 .2];  dp = [.24 0 0];                      % create 4 pseudo sliders
  S.sl = [plt('slider',p+0*dp,[1   0  1],'--- A ---',@sliderCB,[4 .01]);
          plt('slider',p+1*dp,[1   0 20],'--- B ---',@sliderCB,[4  .2]);
          plt('slider',p+2*dp,[0  -4  4],'--- C ---',@sliderCB,[4  .1]);
          plt('slider',p+3*dp,[0 -.5 .5],'--- D ---',@sliderCB,[4 .01],'3 4 3')];
  S.hlp = text(1.55,7.0,...
  {'Use the sliders to change the A,B,C,D parameters';
   '     - Note that sin/cos use the right axis scale';
   '     - The remaining functions use the left axis';},'color',[1 .7 .7]);
  fn = get(gcf,'number');             % get figure number
  if ischar(fn) fn = gcf; end;        % older matlab versions don't store the figure number there
  p = [.01 .50 .07 .04;               % position: help button
       .01 .43 .07 .04];              % position: clipboard button
  set([uicontrol uicontrol],'units','norm',{'pos'},num2cell(p,2),...
     {'str'},{'help'; 'clipboard'},{'callback'},...
     {{@helpBTN,S.hlp}; ['print -f' num2str(fn) ' -dbitmap -noui']});
  set(gcf,'user',S);                  % save structure for slider/cursor callback
  sliderCB;                           % initialize the plot
  plt('cursor',S.cid,'update',81);    % move cursor position
  helpBTN(0,0,S.hlp);                 % enable help text
  tid = findobj(gcf,'user','TraceID');
  rhb = findobj(tid,'xdata',[0 .95]); % find the right hand axis identifying patches
  set(rhb,'color',[.2 0 0]);          % change them to dim red
  tx = findobj(tid,'type','text');    % find all the text objects in the traceID box
  set(tx,'fontname','Courier New','fontsize',11);   % change the font and size
% end function trigplt

function helpBTN(h,arg2,t)  % toggle help visibility
   v = get(t(1),'vis');
   if v(2)=='n' v='off'; else v='on'; end;
   set(t,'vis',v);
%end function helpBTN

function sliderCB()  % Update the plot based on current slider and cursor values
  S = get(gcf,'user');
  if isempty(S) return; end;
  set(S.hlp,'vis','off');
  A = plt('slider',S.sl(1),'get');   B = plt('slider',S.sl(2),'get');
  C = plt('slider',S.sl(3),'get');   D = plt('slider',S.sl(4),'get');
  s = S.tid{plt('cursor',S.cid,'get','activeLine')}; % name of active trace
  set(S.tx,'string',prin('%3w %s(%3wx + %3w) + %3w',A,s,B,C,D));
  xx = pi*[0:.02:2];  x = B*xx+C+1e-12;
  f = A*[sin(x); cos(x); tan(x); csc(x); sec(x); cot(x)] + D;
  f(find(abs(f)>10)) = NaN;
  set(S.tr,'x',xx,{'y'},{f(1,:); f(2,:); f(3,:); f(4,:); f(5,:); f(6,:)})
%end function sliderCB
