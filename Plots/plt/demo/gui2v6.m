% gui2v6.m -----------------------------------------------------------------
%
% This is a slightly modified version gui2.m for running on Matlab
% Ver6.1 since that version doesn't support the uipanel object.
% Here the uipanel was replaced with a frame (uicontrol). demoplt.m
% will detect that version 6 is being used and run this code instead
% of gui2.m

% Unlike the previous gui building example (gui1.m) this one includes a
% plot and actually performes a useful function (displaying the frequency
% response of the most common traditional analog filters). GUI controls
% are provided to adjust the most important filter parameters (Filter
% order, Cutoff frequency, & Passband/Stopband ripple). The capabilities
% of this program where intentionally kept modest to make it easy to learn
% how to use plt to create a plot based graphical interface.
%
% Nine calls to plt are made to create eleven pseudo objects. The first plt
% call creates the first three pseudo objects:
%  1.) a plot
%  2.) a cursor
%  3.) a grid
%  4.) an edit object (filter order)
%  5.) a popup (filter type)
%  6.) a popup (decades to display)
%  7.) a popup (number of points to display)
%  8.) a slider (passband ripple)
%  9.) a slider (stopband ripple)
% 10.) a slider (cutoff frequency)
% 11.) a slider (frequency 2)
%
% Although Matlab already has objects with these names, these pseudo objects
% are different. They provide more utility and  options. The pseudo objects
% 4 thru 7 listed above are grouped inside a uipanel titled "Parameters".
%
% You can most easily absorb the point of this example (and the previous
% one called gui1.m) by watching my video which you can find here:
%   www.mennen.org\plt\video\MatlabGUIbuilding.avi
% Address any questions or comments you may have about this example,
% the video, or plt in general to me at the email address shown below.
%
% ----- Author: ----- Paul Mennen
% ----- Email:  ----- paul@mennen.org

% Although this application uses 10 traces (5 for magnitude on the left axis
% and 5 for phase on the right axis) my first idea was to use just 5 traces
% where each trace included both the magnitude and phase information. Just one axis
% was used, which required rescaling the phase information so that it would appear
% directly above the magnitude plot. The tick marks were modified so they read in
% degrees in the phase portion of the plot which was highlighted with a gray patch
% to visually separate it from the magnitude plot. Despite the required rescaling,
% the 5 trace method was still the simplest way to go. However once I added the feature
% for displaying both magnitude and phase cursors at the same time, it became clear
% that the more straighforward 10 trace method would be shorter because plt already
% includes the "DualCur" feature for cursoring two traces at once. However displaying
% multiple data sets per trace still may be a useful trick on occation, so I have
% included the older 5 trace version of this application in the demo folder
% (gui2ALT.m). The alternate version is not run by demoplt.m.

function gui2v6()
  p = {[.125 .105 .800 .760];  % plot    position
       [.100 .882 .235 .100];  % axis    position: Parameters
       [.165 .935 .040 .030];  % edit    position: filter order
       [.110 .710 .100 .200];  % popup   position: filter type
       [.310 .750 .024 .200];  % popup   position: # of decades
       [.287 .710 .054 .200];  % popup   position: # of points
       [.350 .946 .150     ];  % slider  position: Passband ripple
       [.510 .946 .150     ];  % slider  position: Stopband ripple
       [.670 .946 .150     ];  % slider  position: Cutoff frequency
       [.830 .946 .150     ];  % slider  position: frequency 2
       {-.09 .650          }}; % text    position: eliptic transition ratio
  c = [0 1 0; 1 0 1; 0 1 1; 1 0 0; .2 .6 1];  % trace colors
  S.cfg = [which(mfilename) 'at'];  % use gui2v6.mat to save configuration data
  lbl = {'dB' [blanks(70) 'Phase \circ']};
  typ = {'low pass' 'high pass' 'band pass' 'stop band'};  pts = 100*[1 2 4 8 16];
  S.tr = plt('FigName','gui2v6',0,zeros(1,10),'Right',6:10,'Options','LogX','closeReq',@cfg,...
                  'DualCur',-5,'TraceID',{'Butter' 'Bessel', 'Cheby1'  'Cheby2' 'Elliptic'},...
                  'Ylim',{[-90 60] [-1000 200]},'LabelX','radians/sec','LabelY',lbl,...
                  'TIDcback','t=plt("show"); t=t(find(t<6)); plt("show",[t t+5]);',...
                  'xy',p{1},'TraceC',[c;c],'+Ytick',-140:20:0,'-Ytick',[-180 0 180]);
  ca = gca;
  axes('units','norm','pos',p{2},'box','on','xtick',[],'ytick',[],'color','none',...
       'linewidth',3,'xcolor',[.4 .4 .4],'ycolor',[.4 .4 .4],'tag','frame');
  axes(ca);
  S.n   = plt('edit',  p{3} ,[6 1 25],'callbk',@clb,'label',{'Order:' .05});
  S.typ = plt('pop',   p{4} ,typ,'callbk',@clb,'swap');
  S.dec = plt('pop',   p{5} ,1:5,'callbk',@clb,'index',3,'label','Decades:','hide');
  S.pts = plt('pop',   p{6} ,pts,'callbk',@clb,'index',2,'label','Points:', 'hide');
  S.Rp  = plt('slider',p{7} ,[ 2  .01   9],'Passband ripple', @clb);
  S.Rs  = plt('slider',p{8} ,[ 40  10 120],'Stopband ripple', @clb);
  S.Wn  = plt('slider',p{9} ,[.02 .001  1],'Cutoff frequency',@clb,5,'%4.3f 6 2');
  S.Wm  = plt('slider',p{10},[.2  .001  1],'frequency 2',     @clb,5,'%4.3f 6 2');
  S.etr = text(p{11}{:},'','units','norm','horiz','center','color',[.2 .6 1]);
  S.cid = get(gca,'user');
  set(gcf,'user',S);
  h = getappdata(gcf,'sli'); h(5:5:end) = [];  set(h,'buttond','plt ColorPick;');
  for k = 1:length(h) setappdata(h(k),'m',{'backgr' h}); end;
  if exist(S.cfg) load(S.cfg);  % load configuration file if it exists -------------
                  plt('edit',S.n,'value', cf{1});  plt('pop',S.typ,'index',cf{2});
                  plt('pop',S.pts,'index',cf{3});  plt('pop',S.dec,'index',cf{4});
                  plt('slider',S.Rp,'set',cf{5});  plt('slider',S.Rs,'set',cf{6});
                  plt('slider',S.Wn,'set',cf{7});  plt('slider',S.Wm,'set',cf{8});
                  set(h,'background',     cf{9});  set(gcf,'position',     cf{10});
  end;
  clb;                             % initialize plot
  plt('cursor',S.cid,'update',-1); % center cursor
% end function gui2v6

function clb() % callback function for all objects
  S = get(gcf,'user');
  ty = plt('pop',S.typ,'get');                        % get filter type index
  t = {'low' 'high' 'bandpass' 'stop'};  t = t{ty};   % get filter type name
  N  = plt('edit',S.n,'get');                         % get filter order
  dec = plt('pop',S.dec,'get');                       % get number of decades to plot
  pts = str2num(get(S.pts,'string'));                 % get # of points to plot
  X  = logspace(-dec,0,pts);  W = X*1i;               % X-axis data (radians/sec)
  Wn = plt('slider',S.Wn,'get');                      % get filter freq
  Rp = plt('slider',S.Rp,'get');                      % get passband ripple
  Rs = plt('slider',S.Rs,'get'); Rs2 = max(Rp+.1,Rs); % get stopband ripple (must be > passband)
  if ty>2 Wn = [Wn plt('slider',S.Wm,'get')];         % get frequency 2
          plt('slider',S.Wm,'set','visON');           % make frequency 2 slider visible
  else    plt('slider',S.Wm,'set','visOFF');          % make frequency 2 slider invisible
  end;
  [B,A] = butter(N,Wn,t,'s');        H{1} = polyval(B,W)./polyval(A,W);
  [B,A] = besself(N,Wn(1));          H{2} = polyval(B,W)./polyval(A,W);
  [B,A] = cheby1(N,Rp,Wn,t,'s');     H{3} = polyval(B,W)./polyval(A,W);
  [B,A] = cheby2(N,Rs,Wn,t,'s');     H{4} = polyval(B,W)./polyval(A,W);
  [B,A] = ellip(N,Rp,Rs2,Wn,t,'s');  H{5} = polyval(B,W)./polyval(A,W);
  if ty~=1 H{2}=H{2}+NaN; end;       % bessel filter only applicable for low pass
  for k=1:5   % set trace data
    set(S.tr([k k+5]),'x',X,{'y'},{20*log10(abs(H{k})); angle(H{k})*180/pi});
  end;
  plt('cursor',S.cid,'set','xlim',X([1 end])); % set Xaxis limits
  plt('cursor',S.cid,'updateH');               % update cursor
  h = find(get(S.tr(5),'y') < -Rs2);           % compute Elliptic transition ratio
  if     isempty(h)                    h = 0;
  elseif (ty-2)*(ty-3)                 h = X(h(1))/Wn(1);
  else    h = find(diff([h inf])>1);   h = Wn(1)/X(h(1));
  end;
  set(S.etr,'string',prin('Elliptic ~, transition ~, ratio: ~, %5v',h));
% end function clb

function cfg() % write configuration file
  S = get(gcf,'user');  sli = findobj(gcf,'style','slider');
  cf = { plt('edit',S.n,'get');     plt('pop',S.typ,'get');
         plt('pop',S.dec,'get');    plt('pop',S.pts,'get');
         plt('slider',S.Rp,'get');  plt('slider',S.Rs,'get');
         plt('slider',S.Wn,'get');  plt('slider',S.Wm,'get');
         get(sli(1),'background');  get(gcf,'pos')            };
  save(S.cfg,'cf');
% end function cfg
