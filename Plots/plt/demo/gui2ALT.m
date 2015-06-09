% gui2ALT.m -----------------------------------------------------------------
%
% This is an older version of gui2.m. The differences between the newer and
% older versions are described in the comment section of gui2.m. Although the
% newer version is probably better, this version may still be instructive,
% especially if you ever decide to display multiple data sets on a single trace.
%
% gui2ALT.m creates the same 11 pseudo objects as in gui2.m, so refer to the header
% section of gui2.m for a description of the controls and their functions.

% ----- Author: ----- Paul Mennen
% ----- Email:  ----- paul@mennen.org

function v = gui2ALT(cid)
  if nargin % ystring function - displays phase crosshair and value
    [xy k] = plt('cursor',cid,'get','pos');    % get cursor position
    h = plt('cursor',cid,'get','obj');         % get handles of cursor objects
    hact = h(13);                              % active cursor handle
    lh = get(hact,'user'); y = get(lh{1},'y'); % get line y data
    k = k + floor(length(y)/2) + 1;            % point to phase data
    if k > length(y) v = ''; return; end;      % make sure cursor is in range
    v = prin('%6v\\circ',7.5*y(k)-240);        % scale phase onto magnitude plot
    set(hact,'x',repmat(get(hact,'x'),1,2),... % add phase marker
             'y',[get(hact,'y') y(k)]);
    return;
  end;
  yt1 = -100:20:0;  yt2 = 8:8:56;  b = blanks(27);
  S.tr = plt('FigName','gui2ALT',0,[0 0 0 0 0],'xy',[.14 .1 .84 .77],'Options','LogX',...
    'TraceID',{'Butter' 'Bessel', 'Cheby1'  'Cheby2' 'Elliptic'},'ylim',[-100 60],...
    'ystring','gui2ALT(@CID)','LabelX','radians/sec','LabelY',[b 'dB' b 'Phase (deg)'],...
    '+ytick',[yt1 yt2],'+yticklabel',num2str([yt1 7.5*yt2-240]'));
  S.cid = get(gca,'user');
  patch([1e-9 1 1 1e-9],[5 5 59 59],[1 1 1]/10);
  p = {[.100 .885 .240 .110];  % uipanel position: Parameters
       [.165 .935 .040 .030];  % edit    position: filter order
       [.110 .710 .100 .200];  % popup   position: filter type
       [.310 .750 .024 .200];  % popup   position: # of decades
       [.287 .710 .054 .200];  % popup   position: # of points
       [.350 .946 .150     ];  % slider  position: Passband ripple
       [.510 .946 .150     ];  % slider  position: Stopband ripple
       [.670 .946 .150     ];  % slider  position: Cutoff frequency
       [.830 .946 .150     ];  % slider  position: frequency 2
       {-.11 .560          }}; % text    position: eliptic transition ratio
  typ = {'low pass' 'high pass' 'band pass' 'stop band'};  pts = 100*[1 2 4 8 16];
  uipanel('units','norm','title','Parameters','backgr',get(gcf,'color'),...
                 'pos',p{1},'high',[.4 .4 .4],'foregr',[.4 .4 .4]);
  S.n   = plt('edit',  p{2},[ 6   1   25],'callbk',@clb,'label',{'Order:' .05});
  S.typ = plt('pop',   p{3},typ,'callbk',@clb,'swap');
  S.dec = plt('pop',   p{4},1:5,'callbk',@clb,'index',3,'label','Decades:','hide');
  S.pts = plt('pop',   p{5},pts,'callbk',@clb,'index',2,'label','Points:', 'hide');
  S.Rp  = plt('slider',p{6},[ 2  .01   9],'Passband ripple', @clb);
  S.Rs  = plt('slider',p{7},[ 40  10 120],'Stopband ripple', @clb);
  S.Wn  = plt('slider',p{8},[.02 .001  1],'Cutoff frequency',@clb,5,'%4.3f 6 2');
  S.Wm  = plt('slider',p{9},[.2  .001  1],'frequency 2',     @clb,5,'%4.3f 6 2');
  S.etr = text(p{10}{:},'','units','norm','horiz','center','color',[.2 .6 1]);
  set(gcf,'user',S);
  d = [S.Rp S.Rs S.Wn S.Wm];  h = [];
  for k = 1:4  h = [h plt('slider',d(k),'get','obj')];  h(end) = []; end;
  for k = 1:length(h) setappdata(h(k),'m',{'backgr' h}); end;
  set(h,'buttond','plt ColorPick;');
  clb;
  plt('cursor',S.cid,'update',-1); % center cursor
% end function gui2

function clb(a,b) % callback function for all objects
  S = get(gcf,'user');
  ty = plt('pop',S.typ,'get');                        % get filter type index
  t = {'low' 'high' 'bandpass' 'stop'};  t = t{ty};   % get filter type name
  N  = plt('edit',S.n,'get');                         % get filter order
  dec = plt('pop',S.dec,'get');                       % get number of decades to plot
  pts = str2num(get(S.pts,'string'));                 % get # of points to plot
  W  = logspace(-dec,0,pts);                          % X-axis data (radians/sec)
  Wn = plt('slider',S.Wn,'get');                      % get filter freq
  Rp = plt('slider',S.Rp,'get');                      % get passband ripple
  Rs = plt('slider',S.Rs,'get'); Rs2 = max(Rp+.1,Rs); % get stopband ripple (must be > passband)
  if ty>2 Wn = [Wn plt('slider',S.Wm,'get')];         % get frequency 2
          plt('slider',S.Wm,'set','visON');           % make frequency 2 slider visible
  else    plt('slider',S.Wm,'set','visOFF');          % make frequency 2 slider invisible
  end;
  [B,A] = butter(N,Wn,t,'s');        H1 = frq(B,A,W);         % frequency resp: Butterworth 
  [B,A] = besself(N,Wn(1));          H2 = frq(B,A,W)/(ty==1); % frequency resp: Bessel 
  [B,A] = cheby1(N,Rp,Wn,t,'s');     H3 = frq(B,A,W);         % frequency resp: Chebyshev 1
  [B,A] = cheby2(N,Rs,Wn,t,'s');     H4 = frq(B,A,W);         % frequency resp: Chebyshef 2
  [B,A] = ellip(N,Rp,Rs2,Wn,t,'s');  H5 = frq(B,A,W);         % frequency resp: Elliptical
  set(S.tr,'x',[W NaN W],{'y'},{H1; H2; H3; H4; H5});         % set trace data
  plt('cursor',S.cid,'set','xlim',W([1 end]));                % set Xaxis limits
  h = find(-H5(1:pts) > Rs);                                  % compute Elliptic transition ratio
  if isempty(h)                      h = NaN;   
  elseif (ty-2)*(ty-3)               h = W(h(1))/Wn(1);
  else   h = find(diff([h inf])>1);  h = Wn(1)/W(h(1));
  end;
  set(S.etr,'string',prin('Elliptic ~, transition ~, ratio: ~, %5v',h));
% end function clb

function H = frq(B,A,W) % computes frequency response dB magnitude and phase
  H = j*W;  H = polyval(B,H) ./ polyval(A,H);
  H = [20*log10(abs(H)) NaN 32+angle(H)*24/pi];
