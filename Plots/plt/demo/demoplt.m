% demoplt.m -----------------------------------------------------------
%
% This program makes it easy to run any of the 25 sample programs in
% the plt/demo folder by clicking on the button labled with the
% program name. After installing plt for the first time, it is
% advantageous to run demoplt and to click on the "All Demos" button
% which will cause demoplt to run all 25 sample programs in sequence.
% Simply click on the "continue" button when you are ready to look at
% the next program in the list. By viewing the plots created by all
% the demos you will quickly get a overview of the types of plots that
% are possible using plt. This often gives you ideas about how you can
% use plt for creating your own application. plt5 is first on the list
% because it is the simplest most basic example. The other demos
% appear in alphabetical order.
%
% As each demo is run, you may peruse the code for the demo program
% currently being run in the demoplt list box. Use the list box
% scroll bars to view any portion of the code of interest. If the text
% is to big or small for comfort, adjust the fontsize using the
% fontsize popup menu in the lower right corner of the demoplt figure.
% This fontsize is saved (in demoplt.mat) along with the current
% figure colors and screen location so that the figure will look the
% same the next time demoplt is started. (Delete demoplt.mat to return
% to the original default conditions.)
%
% If you are running a version of Matlab older than 7.0 then the gui1
% button is replaced by the gui1v6 button. gui1.m uses a uitable which
% aren't supported in Matlab 6, so gui1v6 replaces the uitable with a
% radio button. Similarly gui2 is replaced by gui2v6 if you are running
% a version of Matlab older than 7.0 or if you are running version
% 8.4 (R2014b). gui2 uses a uipanelwhich isn't supported in Matlab 6,
% so gui2v6 replaces the uipanel with a uicontrol frame which serves
% pretty much the same function. R2014b supports the uipanel, but the
% v6 version is run because of a bug in R2014b relating to the stacking
% order of a uipanel.
%
% In addition to its main role as a demo program launcher, demoplt
% demonstrates the use of one of plt's pseudo objects, namely the
% ColorPick window. (Note: A pseudo object is a collection of more
% primitive Matlab objects, assembled together to perform a common
% objective.) The ColorPick pseudo object is useful whenever you want
% to allow the user to have control over the color of one of the graphic
% elements. In demoplt there are 4 such elements: The text color, the
% text background color, the button color, and the figure background
% color. The ColorPick window is activated when you left or right click
% on any of the 4 text labels near the bottom of the figure, or any of
% the uicontrol objects adjacent to these labels. (A right click on
% the figure background edit box will bring up the ColorPick window).
% When the ColorPick window appears you can use the sliders or the color
% patches to change the color of the respective graphic element. For more
% details, see the "Pseudo objects" section in the help file.
%
% An optional feature of the ColorPick object is the color change
% callback function (a function to be called whenever a new color is
% selected). This feature is demonstrated here by assigning the color
% change callback to "demoplt(0)". The callback function is used to
% "animate" the figure by briefly modifying the listbox fontsize.
% This particular animation is certainly not useful (see line 111 to
% disable it), but is there merely to demonstrate how to use the callback.
%
% Although it's unrelated to plt, demoplt also demonstrates the use of
% the close request function, which in this example is assigned to
% demoplt(1) and is used to save the currently selected colors and
% screen position to a setup file.

%
% ----- Author: ----- Paul Mennen
% ----- Email:  ----- paul@mennen.org

function demoplt(in1)
  fg = findobj('name','demoplt');
  cfile = [which(mfilename) 'at'];  % use demoplt.mat to save figure colors
  if nargin
    bx = findobj(fg,'style','listbox');    pu = findobj(fg,'style','popup');
    if in1  % here for close request function
      close(findobj('name','Color Pick'));  % close any open ColorPick windows
      fbk = get(fg,'color');  lbk = get(bx,'backgr');  lfg = get(bx,'foregr');
      bbk = get(pu,'backgr'); fpos = get(fg,'pos');   fsz = get(bx,'fontsize');
      save(cfile,'fbk','lbk','lfg','bbk','fpos','fsz','-v4'); % save the .mat setup
      closereq; % close the demoplt figure
    else        % here to do the color callback function (listbox animation)
      if length(get(bx,'tag')) return; end;  % prevent callback recursion
      set(bx,'tag','0'); fs = get(bx,'fontsize');  k = fs+1;
      while k~=fs % rotate font size until we are back where we started
        k=k+1; set(bx,'fontsize',k); pause(.01); if k>28 k=1; end;
      end;
      set(bx,'tag','');  % re-enable the animation callback
    end;
    return;
  end;
  % come here if demoplt is called without arguments
  if length(fg) disp('demoplt is already running'); return; end;
  if exist(cfile) load(cfile);       % get setup infromation from .mat if it exists
  else            fbk = [ 0 .2 .3];  % Otherwise: figure background color
                  lbk = [ 0  0  0];  %            listbox background color
                  lfg = [ 1  1 .5];  %            listbox foreground (text color)
                  bbk = [.4 .8 .8];  %            button background color
                  fpos = [];         %            use default figure position
                  fsz = 9;           %            default listbox fontsize
  end;
  set(0,'units','pixels'); ssz = get(0,'screensize');  % get screen size in pixels
  w = ssz(3);      % position the demoplt figure near the right edge of the screen
  fg = figure('name','demoplt','menu','none','numberT','off','double','off',...
              'pos',[w-660 50 655 420],'color',fbk,'closereq','demoplt(1)');
  v = version;  v = str2num(v(1:3));
  if v<7 | v>= 8.4 gui2 = 'gui2v6'; else gui2 = 'gui2'; end;
  if v<7           gui1 = 'gui1v6'; else gui1 = 'gui1'; end;
  fcn = {'plt5'    'bounce'  'circles12' 'curves C' 'dice(0)'   ... 
         'editz'   'gauss'    gui1        gui2      'movbar(1)' ...
         'plt50'   'pltn'    'pltquiv'   'pltsq(1)' 'pltvar'    ...
         'pltvbar' 'pub'     'pub2'      'subplt'   'subplt8'   ...
         'tasplt'  'trigplt' 'weight'    'wfall(1)' 'winplt'       };
  n = length(fcn);  bt = zeros(1,n+1);
  x = 8;  y = 392;  % starting location for the first demo program button
  for k=1:n         % create all 25 buttons
     bt(k) = uicontrol('Pos',[x y 65 22],'string',fcn{k},'callb',{@OneDemo,fcn{k}});
     x = x+72;  if x>600 x=8; y=y-30; end;
  end;
  bt(end) = uicontrol('Pos',[x y 137 22],'user',0,'str','All Demos','callback',{@AllDemos,fcn});
  set(bt,'units','norm','fontsize',10,'backgr',bbk);
  set(bt(end),'backgr',fliplr(bbk)); % use a different color for the All Demos button
  bx = uicontrol('style','listbox','units','norm',...  % listbox for viewing demo code
          'pos',[.01 .076 .98 .7],'string',rdfunc('demoplt'),'fontsize',fsz,...
          'fontname','Lucida Console','foregr',lfg,'backgr',lbk);
  axes('pos',[1 2 6 1]/10,'units','norm'); % create an axis to contain the text objects
  cb = 'plt ColorPick demoplt(0);';        % callback for text objects, frames & edit box
  % cb = 'plt ColorPick;'; % enable this line to disable the ColorPick animation
  u = [text(.004,-1.6,'text color:') ...
       text(.373,-1.6,'text background:') ...
       text(.677,-1.6,'button color:') ...
       text(1.07,-1.6,'figure background:')];
  set(u,'color','white','horiz','right','buttondown',cb);
  p = [.104 .023 .044 .029;  % a(1) position: text color frame
       .326 .023 .044 .029;  % a(2) position: text background frame
       .509 .023 .044 .029;  % a(3) position: button color frame
       .745 .019 .101 .040;  % a(4) position: figure background editbox
       .858 .014 .130 .052]; % a(5) position: fontsize popup
  a = [uicontrol('enable','inactive') ...
       uicontrol('enable','inactive') ...
       uicontrol('enable','inactive') ...
       uicontrol('string',prin('{%3w!  }',fbk),'fontW','bold','callback',cb) ...
       uicontrol('string',prin('{fontsize: %d!row}',4:18),'value',fsz-3,'fontsize',9,...
         'callback','set(findobj(gcf,''sty'',''l''),''fontsize'',3+get(gcbo,''value''));')];
  set(a,{'style'},{'frame';'frame';'frame';'edit';'popup'},'buttondown',cb,'units','norm',...
        {'pos'},num2cell(p,2),{'backgr'},{lfg; lbk; bbk; bbk; bbk});
  q = {{'backgr',a(1),'foregr',bx,'text color'};  % define color pick actions
       {'backgr',a(2),'backgr',bx,'text background'};
       {'backgr',a(3),'backgr',[a(4:5) bt],'button color'};
       {'string',a(4),'color', fg,'figure background'}};
  for k=1:4 setappdata(a(k),'m',q{k}); setappdata(u(k),'m',q{k}); end;
  if isempty(fpos) & w>1030
    fpos = [w-700 50 690 480];  % use a bigger demoplt figure for large screen sizes
  end;
  if length(fpos) set(fg,'pos',fpos); end;
  set(fg,'user',sum(get(fg,'pos')));
%end function demoplt

function s=rdfunc(func) % returns a cell array of strings from a program file
  f = findstr(func,' ');
  if length(f) func = func(1:f-1); end;
  f = fopen(which(func));  s = {};
  while 1  ln = fgetl(f);
           if ~ischar(ln), break, end
           s = [s; {ln}];
  end
  fclose(f);
%end function rdfunc

function Fclose % close all figures created by demoplt or plt
  a = findobj('type','figure');
  for k=1:length(a)
    b = a(k);
    if ishandle(b) & (isappdata(b,'cid') | isappdata(b,'demo')) ...
                   & ~strcmp(get(b,'Name'),'demoplt')
       close(b);
    end;
  end;
%end function Fclose

function OneDemo(h,arg2,func)
  Fclose;
  set(findobj(gcf,'style','push'),'fontw','normal');
  set(findobj(gcf,'string',func),'fontw','bold');
  set(findobj(gcf,'style','listbox'),'val',1,'string',rdfunc(func));
  if strcmp(func,'pltvar') evalin('base',func); else eval([func ';']); end;
  setappdata(gcf,'demo',0);
  a = get(0,'child');                % list of all the figure windows
  b = findobj(a,'name','demoplt');   % find demoplt in the list
  if length(b) set(0,'child',[b; a(find(a~=b))]); end; % force demoplt window on-top
%end function OneDemo 

function AllDemos(h,arg2,fcn)
  nf = length(fcn);
  n = get(h,'user') + 1;  % function to start
  if n>nf
    set(h,'user',0,'string','All Demos');
    Fclose;
    set(findobj(gcf,'style','push'),'fontw','normal');
    set(findobj(gcf,'style','listbox'),'string',{'' '' '   --- ALL DEMOS COMPLETED ---'});
  else set(h,'user',n,'string','continue');
       OneDemo(0,0,fcn{n});
  end;
%end function AllDemos
