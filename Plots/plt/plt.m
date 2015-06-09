% plt.m:   An alternative to plot & plotyy (version 13Mar15)
% Author:  Paul Mennen (paul@mennen.org)
%          Copyright (c) 2015, Paul Mennen

function [Ret1,Ret2] = plt(varargin);

global Hhcpy;

Ret1 = 0;
Narg = nargin;
Mver = version;  Mver = s2d(Mver(1:3));
if Mver < 8.4  ERAS = 'eras';  ERAXOR = 'xor';  ERANOR = 'norm';
else           ERAS = 'pi';    ERAXOR = 'v';    ERANOR = 'v';
end;

if ~Narg
  eval('evalin(''base'',''pltdef'','''')');
  w = plt('select','base','who');
  n = length(w);  b = zeros(1,n);
  options = [];
  wstruct = {}; ws=0;
  TraceIDlength = 7;
  for k=1:n
    a = plt('select','base',w{k});
    if length(w{k})>5 & findstr(w{k},'pltopt') options = [options a]; end;
    if strcmp(w{k},'TraceIDlen') TraceIDlength = a;  end;
    if isnumeric(a) & length(a)>3 b(k)=1;
    elseif isstruct(a) & numel(a)==1
      f = fieldnames(a);
      for fn=1:length(f)
        g = getfield(a,f{fn});
        if isnumeric(g) & length(g)>3 ws=ws+1; wstruct{ws} = [w{k} '.' f{fn}]; end;
      end;
    end;
  end;
  w = [w(logical(b)); wstruct'];  nov = length(w);
  if nov<2 disp('Not enough vectors to plot'); return; end;
  spix = get(0,'screenp');
  fontsz = (196-spix)/10;
  spix = spix/80 - .2;
  vspace = round(20 * spix);
  hspace = round(200 * spix);
  set(0,'unit','pix');  ssz = get(0,'screensize');
  nc = fix((ssz(4)-110)/vspace);
  nov1 = nov+3;  nc = ceil(nov1/nc);
  nr = ceil(nov1/nc);
  wid = nc*hspace+20; hei = nr*vspace+55;
  fig = figure('menu','none','numberT','off','back','off','name','PLT select','color','black',...
               'Invert','on','pos',[ssz(3)-wid-5 ssz(4)-hei-76 wid hei]);
  axes('xlim',[0 nc],'ylim',[0 nr+1]-.5,'pos',[.02/nc 0 1 1],...
       'color','black','xcol','black','ycol','black');
  hb = hei-25;
  set([uicontrol('str','Plot',    'pos',[ 25 hb 40 21],'call','plt select plot;');
       uicontrol('str','CloseAll','pos',[ 80 hb 60 21],'call','plt close;');
       uicontrol('str','Help',    'pos',[155 hb 40 21],'call','plt help;');
      ],'unit','nor','fontsi',fontsz);
  text(.05,nr-2.2,{'Right click on desired x vectors';
                   'Left click on desired y vectors';
                   'Double click for right hand axis'},...
                   'fonta','ita','fontsi',9,'color',[0 1 .5]);
  txt = zeros(nov,1); sz1=txt; sz2=txt; len=txt; nvec=txt; sxy=txt; yright=txt;
  y=nr-3; x=0;
  for k=1:nov
    y = y-1;  if y<0 y=nr-1; x=x+1;  end;
    sz = size(plt('select','base',w{k}));
    sz1(k)  = sz(1);
    sz2(k)  = sz(2);
    len(k)  = max(sz);
    nvec(k) = min(sz);
    txt(k) = text(x,y,' ');
    set(txt(k),'buttond',['plt(''select'',' int2str(k) ');']);
  end;
  set(fig,'user',[txt sz1 sz2 len nvec sxy yright]);
  setappdata(fig,'w',w);
  setappdata(fig,'opt',options);
  setappdata(fig,'IDlen',TraceIDlength);
  plt select text;
  return;
end;

y1 = varargin{1}; 
if length(y1)==1 & ishandle(y1)
  y2 = varargin{2};
  if isempty(y2) | ~isnumeric(y2)
    hcback = y1;
    varargin(1:2)=[]; Narg=Narg-2;
    y1 = varargin{1};
  end;
end;
if ~ischar(y1)  y1='none'; end;
if y1(1)=='g' y1(1)='G';  end;
switch sum(y1)
case 640
  v2 = varargin{2};
  if isnumeric(v2) | v2(1) ~= 'b'
      w = getappdata(gcf,'w');  e = get(gcf,'user');
      len = e(:,4); sxy = e(:,6); nov = length(len);
      xk = []; for k=1:nov if sxy(k)<0 xk = [xk k]; end; end;
  end;
  switch sum(v2)
    case 411, Ret1 = eval(['evalin(''base'',''' varargin{3} ''');']);
    case 453,
      lnxk = length(xk);
      for k=1:nov
        s = w{k};  if length(s)>20 s=s(1:20); end;
        s = strrep(s,'_','\_');
        s = [s '  (' int2str(e(k,2)) ',' int2str(e(k,3)) ')'];
        if sxy(k) s = [s ' \leftarrow ']; s2 = int2str(abs(sxy(k)));  end;
        c = [1 1 1];
        if lnxk
          if lnxk==1 s2 = ''; end;
          if ~sxy(k)
            if ~length(find(len(xk)==len(k))) c = .5*c; end;
          elseif sxy(k)<0                 c = [.9 .3 .3]; s = [s 'x' s2];
          else
            if e(k,7)               c = [ 1 .6 .2]; s = [s 'y' s2 'R'];
                      else                c = [ 1  1  0]; s = [s 'y' s2];
                      end;
          end;
        else if e(k,5) ~= 1 c = .5*c; end;
        end;
        if c(3)==.5 | c(3)==1 bold = 'norm'; else bold = 'bold';  end;
        set(e(k,1),'color',c,'fontw',bold,'str',s);
      end;
    case 447
      if ~max(sxy) return; end;
      options       = getappdata(gcf,'opt');
      TraceIDlength = getappdata(gcf,'IDlen');
      opt = [];
      for k=1:length(options)
        t = options{k};
        if   isnumeric(t) v = num2str(t);
             if length(t)>1
               q = v;  v = '[';
               for j=1:size(q,1) v = [v q(j,:) ';']; end;
               v = [v ']'];
             end;
        else v = ['''''' t ''''''];
        end;
        opt = [opt v ','];
      end;
      if length(xk)==1 x = w{find(sxy<0)}; else x = 'X axis'; end;
      right = [];
      ny = 0;
      nvL = 0; nvR = 0;
      func = ['plt(''''LabelX'''',''''' x ''''','];
      id = '''''TraceID'''',[''''';
      for k=1:nov
        if sxy(k)>0
          nvec = e(k,5);
          x = w{find(sxy==-sxy(k))};  y = w{k};
          s2 = strrep(y,'_','\_');
          ss = strrep(y,'_','');
          nn = TraceIDlength - (nvec>1);  s = [];
          if length(ss) > nn
            sp = findstr(ss,'.');
            if length(sp)
              nn1 = min(sp-1,floor((nn-1)/2));
              nn2 = min(nn-nn1-1,length(ss)-sp);
              nn1 = nn - nn2 - 1;
              ss = ss([1:nn1-1 sp-1:sp+nn2-1 end]);
            else
              ss = [ss(1:nn-1) ss(end)];
            end;
          end;
          while length(ss) < nn ss = [ss ' ']; end;
          st = ss;
          for v = 1:nvec
            if nvec>1 ss = [st '0'+mod(v,10)]; end;
            s = [s ss ''''';'''''];
          end;
          if e(k,7) right = [right ny+(1:nvec)];  s2R=s2; nvR=nvR+1;  else  s2L=s2; nvL=nvL+1;  end;
          ny = ny + nvec;
          func = [func x ',' y ','];
          id = [id s];
        end;
      end;
      if nvL==1 func = [func '''''LabelY'''','''''  s2L ''''',']; end;
      if nvR==1 func = [func '''''LabelYr'''',''''' s2R ''''',']; end;
      if nvR>0  func = [func '''''Right'''',' prin('[{%d }]',right) ',']; end;
      if ny>1 func = [func id(1:end-3) '],']; end;
      func = [func opt];
      func = [func(1:end-1) ');'];
      plt('select','base',func);
    otherwise,
      click = sum(get(gcf,'SelectionT'));  sk = sxy(v2);
      if click==321 | ~length(xk)
        if sk>=0
             if e(v2,5)==1 e(v2,6) = min(sxy)-1; end;
        else e(v2,6) = 0;
             e(find(sxy==-sk),6)=0;
        end;
      elseif sk>=0
        if click==434 & sk sk = 0; end;
        xki = sort(-sxy(xk(find(len(xk)==len(v2)))));
        if length(xki)
          if sk m = find(xki==sk)+1;
                if m>length(xki) sk=0; e(v2,7)=0; else sk=xki(m);  end;
          else  sk = xki(1);
          end;
          e(v2,6) = sk;
          if click==434 & sk  e(v2,7)=1; end;
        end;
      end;
      set(gcf,'user',e);
      plt('select','text');
  end;

case 756
    x = varargin{2};
    m=0;  nv=x;
    if ~x Ret1=''; Ret2=1; return; end;
    while nv>900  m=m+1;  nv=.001*nv; end;
    while nv<.9   m=m-1;  nv=1000*nv; end;
    if (m==1 & x<10000) | (m==-1 & x>=.1 ) m=0; nv=x; end;
    Ret2 = nv/x;
    if     m>5  Ret1 = Pftoa('%4w',1/Ret2);
    elseif m<-6 Ret1 = Pftoa('%5w',1/Ret2);
    else        qs=reshape('Atto FemtoPico Nano MicroMilli     Kilo Mega Giga Tera Peta ',5,12)';
                Ret1 = deblank(qs(m+7,:));
    end;
    if length(Ret1) Ret1 = [Ret1 '-']; end;

case 759
  c = datevec(varargin{2});
  if isempty(c) | isnan(c(1)) Ret1 = []; return; end;
  if Narg<3 fmt=0; else fmt=varargin{3}; end;
  frac = sprintf('%4.3f',mod(c(6),1));
  c(6) = floor(c(6));
  if strcmp(frac,'1.000')
    if c(6)==59 frac = '0.999'; else c(6) = c(6)+1;  end;
  end;
  Ret1 = [datestr(datenum(c),fmt) frac(2:end)];
  if strcmp(Ret1(7:9),'-20') Ret1(8:9) = []; end;

case 425
  if Narg==1       f = '';
  elseif sum(get(gcf,'SelectionT'))==649 f = get(gcbo,'user');
  else             f = get(gcbo,'tag');
  end;
  foobar = 0;
  if exist('foobar')
    if isempty(f)
      fc = feval('which','plt.chm'); nc=length(fc);
      fh = feval('which','plt.htm'); nh=length(fh);
      if (nc+nh)==0
         disp('plt.chm and plt.htm were not found');
         return;
      end;
      if ispc if nc f=fc; else f=fh; end;
      else    if nh f=fh; else f=fc; end;
      end;
    else
      if isempty(findstr(f,filesep))
        f = feval('which',f);
      end;
      if ~exist(f) prin(1,'File %s not found \n',f);
                   return;
      end;
    end;
    if ispc & strcmpi(f(end-2:end),'chm')
         dos(['hh ' f ' &']);
    else feval('web',['file:///' f ],'-browser');
    end;
  else
    if isempty(f) f = 'plt.chm'; end;
    if isempty(findstr(f,filesep)) f = fullfile(fileparts(GetExe),f); end;
    if ~exist(f) return; end;
    if strcmpi(f(end-2:end),'chm') dos(['hh ' f ' &']); else ibrowse(f); end;
  end;

case 335

  k = 2; b = varargin;  a1 = b{2};  nvar = length(b);
  if isnumeric(a1) & ~ischar(b{3})
     a1 = 'pos'; b = [b(1) {a1} b(2) {'choices'} b(3:nvar)]; nvar = nvar+2;
  end;
  if ischar(a1)
       g  = gca;
       bk = [0 .3 .4];
       a = axes('vis','of','XtickLabel',' ','YtickLabel',' ','TickLen',[0 0]','color',bk,'xcol',bk,'ycol',bk);
       pd = {a     ''       1      ''      99    [.1 1 .9]     1      []   'none'   -1   0 };
       Ret1 = text(0,0,'','units','data','user',pd,'interp','none','color',[1 1 .4]);
       setappdata(gcf,'pop',[getappdata(gcf,'pop') Ret1]);
       setappdata(Ret1,'ty','pop');  setappdata(a,'ty','pop');
       set(Ret1,'buttond',{@plt 'pop' Ret1 'open'});
       set(a,'user',Ret1);
       axes(g);
  else n = length(a1);  Ret1 = a1;
       if n==1 k=k+1; pd=get(Ret1,'user'); a=pd{1};
       else    argn = b{3};  argv = b{4};  s = size(argv,1);
               if s for k=1:n plt('pop',Ret1(k),argn,argv(min(k,s),:)); end;
               else for k=1:n plt('pop',Ret1(k),argn,[]); end;
               end;
               return;
       end;
  end;
  nvar = length(b);  ofs = [];  cb = 0;
  while k <= nvar
    argn  = lower(b{k}); k=k+1;
    if k>nvar argv=''; else  argv = b{k}; k=k+1; end;
    switch sum(argn)
    case {338,885},
                    if argv(1)<0 argv(1)=-argv(1); pd{5}=0; end;
                    if argv(2)<0 argv(2)=-argv(2); pd{11}=1; edg(Ret1,pd{6}); end;
                    if max(argv)>1 u='pixels'; else u='normal'; end;
                    if length(argv)<2 fp=get(a,'pos');  argv=[fp(1:3) argv]; end;
                    set(a,'pos',argv,'units',u);
    case 734,  if isnumeric(argv)
                      if sum(mod(argv,1)) fm = '{%5w!row}'; else fm = '{%d!row}'; end;
                      argv = prin(fm,argv);
                    end;
                    pd{2} = argv;
    case 434,     if ~pd{7} return; end;
                    mv = get(Ret1,'buttond');
                    if iscell(mv) mv = mv{2}(1); else mv = mv(5); end;
                    mv = (mv=='m');
                    if xor(sum(get(gcf,'SelectionT'))==649,pd{11}) & ~(sum(get(gcf,'SelectionT'))==434) | mv
                      uh = findobj(pd{8},'vis','on');
                      setappdata(Ret1,'Uhide',uh);
                      set([Ret1; uh],'vis','of');  ch = pd{2};  n = length(ch);
                      set(a,'vis','on');  axes(a);
                      s = {@plt 'pop' Ret1 'index'};  w = get(gcf,'pos');
                      p = get(a,'pos');  p = p(3);  if p<1 p = p*w(3); end;
                      p = 5/p;
                      for m=1:n 
                          ht = text(p,n+.5-m,ch{m},'units','data');
                          set(ht,'interp',pd{9},'color',pd{6},'buttond',[s {-m}]);
                          if m==pd{3} set(ht,'fontw','bol'); end;
                      end;
                      cr = get(a,'color');  z = zeros(1,n-1);
                      x = [z;z+1;z];  z = [z;z;z+NaN];   y = 1:n-1;  y = [y;y;y];
                      line('z',z(:),'y',y(:),'x',x(:),'color',cr-.2+.4*(cr<.5),'user',w);
                    else
                      n = length(pd{2});  v = pd{3};
                      rpt = getappdata(gcf,'repeat');
                      if length(rpt)>1 p = rpt(2); else p = .4; end;
                      if length(rpt) rpt=rpt(1); else rpt = .03; end;
                      g = gcf;  setappdata(g,'bdown',1);
                      set(g,'WindowButtonUp','setappdata(gcf,''bdown'',0);');
                      while ishandle(g) & getappdata(g,'bdown')
                        if v==n v=1; else v=v+1; end;
                        plt('pop',Ret1,'index',-v);
                        q=p; while q>0 & ishandle(g) & getappdata(g,'bdown') pause(.01); q=q-.01; end;
                        p = rpt;
                      end;
                    end;
                    return;
    case 410,     if isempty(argv) argv = 0; end;
                    z = find(argv==0);
                    if length(z) argv(z) = findobj(gcf,'user','grid'); end;
                    pd{8} = argv;
    case 647,   pd{5} = argv;
    case 759,  pd{6} = argv;
    case 748,  set(a,'color',argv,'xcol',argv,'ycol',argv);
    case 658,   pd{9} = argv;  set(Ret1,'interp',argv);
    case 443,     pd{11} = 1;  if isempty(argv) argv = 'none'; elseif length(argv)==1 argv = pd{6}; end;
                    edg(Ret1,argv);
    case 536,    w = abs(argv);  pd{3} = w;   v = get(a,'vis');
                    if v(2)=='n'
                      set(a,'vis','of');  ch = get(a,'child'); li = findobj(ch,'type','line');
                      dp = sum(abs(get(li,'user') - get(gcf,'pos')));
                      if dp>0 & dp<20
                        pd{11} = 1 - pd{11};
                        disp('swap toggled');
                      end; 
                      ch(find(ch==Ret1)) = [];  ch(find(ch==pd{10})) = [];
                      li = getappdata(Ret1,'li');  if ishandle(li) ch(find(ch==li)) = []; end;
                      delete(ch);
                      uh = getappdata(Ret1,'Uhide');
                      set([Ret1; uh],'vis','on');
                      ch = gcf;
                      matVer = version;  mver = (matVer([1 3])-'0') * [10; 1];
                      if mver<75 & matVer(2)=='.' & matVer(4)=='.'
                        set(ch,'Share','off');
                        set(ch,'Share','on');
                      end;
                    end;
                    cb = (argv<0);
    case 617,   pd{4} = argv;  set(Ret1,'user',pd);
    case 320,      switch sum(argv)
                      case 663,  Ret1 = get(Ret1,'str');
                      case 437,      Ret1 = pd{1};
                      case 734,   Ret1 = pd{2};
                      case {536,0}, Ret1 = pd{3};
                      case 617,    Ret1 = pd{4};
                      case 647,    Ret1 = pd{5};
                      case 759,   Ret1 = pd{6};
                      case 615,    Ret1 = pd{7};
                      case 410,      Ret1 = pd{8};
                      case 658,    Ret1 = pd{9};
                      case 512,     Ret1 = pd{10};
                      case 443,      Ret1 = pd{11};
                      otherwise,       Ret1 = pd;
                    end;
                    return;
    case 615,   pd{7} = argv;
    case {512 633},
                    if ischar(argv) argv = {argv ''}; end;
                    s = argv{1};
                    q = argv{2};
                    if isempty(q)
                      if sum(argn)==512 q = -5;
                      else                 q = 16i;
                      end;
                    end;
                    q = [real(q) imag(q)];
                    if ~q(2) & q(1)<0 hz = 'right'; else hz = 'left'; end;
                    pd{10} = text(0,0,s,'units','data','horiz',hz,'color',[.6 .6 .8],'par',a,'user',q);
                    j=3;  while length(argv)>j set(t,argv{j},argv{j+1}); j=j+2; end;
    otherwise,      set(Ret1,argn,argv);
    end;
  end;
  c = pd{2};  n = length(c);  f = pd{5};  i = pd{3};  s = c{i};
  if length(f)<2 f=[.08 f]; end;
  if f(2)==99 f(2)=n; end;
  set(a,'ylim',[0 n]);
  set(Ret1,'pos',f,'str',s,'user',pd);
  edg(Ret1);
  p = get(a,'pos');
  t = pd{10};
  if ishandle(t)
    e = p(3:4);
    if p(1)<1  r = get(gcf,'pos'); e = e .* r(3:4); end;
    q = get(t,'user')./e;  q(2)=q(2)*n;
    set(t,'pos',f+q);
  end;
  p(4) = p(4)/n;
  setappdata(Ret1,'ppos',p - [([0 .3]-f) .4 .4] .* p([3 4 3 4]));
  if cb evalRep(pd{4},{'@STR',s,'@IDX',int2str(i)}); end;

case 422

  k = 2; b = varargin;  a0 = b(1); a1 = b{2};  nvar = length(b);
  if isnumeric(a1) & isnumeric(b{3})
     if     max(a1)<4      a0 = {'edit' 'units' 'norm'};
     elseif length(a1)==2  a0 = {'edit' 'units' 'pix'};
     end;
     a1 = 'pos';  b = [a0 {a1} b(2) {'value'} b(3:nvar)];
     nvar = length(b);
  end;
  if ischar(a1)
       ui = 0;
       for z = 2:2:nvar  d = b{z};
                         if lower(d(1))=='p' ui = length(b{z+1}); break; end;
       end;
       if ~ui disp('You must specify a position for the edit pseudo object'); return; end;
       ui = (ui==4);

       pd = { 1   -1e9  1e9     ''      1    1     1  '%6w'  -1   ''    ''  };

       c = .8*get(gcf,'color');  f = [1 1 .4];  if sum(c)>1.5 f = 1-f; end;
       if ui ut='uic'; Ret1 = uicontrol('style','text','fontsi',9,'ena','inact','foregr',f,'backgr',c);
       else  ut='txt'; Ret1 = text(0,0,'1','horiz','cent','color',f);
       end;
       setappdata(gcf,ut,[getappdata(gcf,ut) Ret1]);
       setappdata(Ret1,'ty','edi');
       set(Ret1,'user',pd,'buttond',{@plt 'edit' Ret1 'click' 0});
  else n = length(a1);  Ret1 = a1;
       if n==1 k = k + 1;   pd = get(Ret1,'user');
       else    argn = b{3};  argv = b{4};  s = size(argv,1);
               if s for k=1:n plt('edit',Ret1(k),argn,argv(min(k,s),:)); end;
               else for k=1:n plt('edit',Ret1(k),argn,[]); end;
               end;
               return;
       end;
  end;
  ui = get(Ret1,'type');  ui = (ui(1) == 'u');
  if ui cp = 'foregr'; gc = gcf;   nun = 'pix';
  else  cp = 'color';    gc = gca;   nun = 'data';
  end;
  while k <= nvar
    argn  = lower(b{k}); sargn = sum(argn); k=k+1;
    if k>nvar argv=''; else  argv = b{k}; k=k+1; end;
    switch sargn
    case 518,
      if argv
        c = get(gcf,'CurrentChar');
        if ~length(c) c=666;  end;
        if argv==2 c=27; end;
        t = get(Ret1,'str');   u = get(Ret1,cp);
        p = findstr(t,'_');  n = length(t);   done = 0;  eeval = 0;
        switch abs(c)
        case 666,
        case 8,    if p>1 t(p-1)=[]; end;
        case 127,  if p==n t='_';
                   else    t(p+1)=[];
                   end;
        case 127,  if p==n t='_'; else t(p+1)=[];  end;
        case 28,   if p>1 t(p)=t(p-1); t(p-1)='_'; end;
        case 29,   if p<n t(p)=t(p+1); t(p+1)='_'; end;
        case 27,   done = 1;  t = pd{10};
        case 13,   done = 1;
                   t(p) = [];
                   if pd{7}
                     tt=t;  t = 'error';
                     try, s = eval(tt);
                        if pd{7}==-1 | pd{7}==length(s)
                          if length(s)==1
                            if isreal(s)
                                 s = max(min(s,pd{3}),pd{2});
                                 t = ['  ' Pftoa(pd{8},s) '  '];
                            else pd{6} = imag(s);
                                 t = pd{10};
                            end;
                          else  t = ['[' num2str(s) ']'];
                          end;
                          if isreal(s)  pd{1} = s;  eeval = 1; end;
                        end;
                     end;
                   else pd{1}=0; eeval=1;
                   end;
        otherwise, t = strrep(t,'_',[c '_']);
        end;
        set(Ret1,'str',t);
        if done
          set(Ret1,'user',pd,cp,u([3 1 2]))
          if ~ui set(Ret1,ERAS,ERANOR,'interp',pd{11}); end;
          set(gcf,'keypress','');
          if eeval setappdata(gcf,'OBJ',Ret1); evalRep(pd{4},{'@VAL' t '@OBJ' 'plt("misc",3)'}); end;
        end;
      else
        kp = get(gcf,'keypress');  f = findstr(kp,'click');
        if length(f)==1  kp(f+7) = '2';  eval(kp); return; end;
        if ~pd{5} return; end;
        c = sum(get(gcf,'SelectionT'));
        if c ~= 321
          v = pd{1}; w = pd{6};
          if pd{7}==1 & w
            cp = get(gc,'currentp'); cp = cp(1,1);
            rp = get(Ret1,'pos');  un = get(Ret1,'unit');
            if un(1) ~= nun(1)
              set(Ret1,'units',nun);
              rp2 = rp;  rp = get(Ret1,'pos');
              set(Ret1,'units',un,'pos',rp2);
            end;
            cp = cp > rp(1) + ui*rp(3)/2;
            cp = 2*cp - 1;
            rpt = getappdata(gcf,'repeat');
            if length(rpt)>1 p = rpt(2); else p = .4; end;
            if length(rpt) rpt=rpt(1); else rpt = .03; end;
            setappdata(gcf,'bdown',1);
            set(gcf,'WindowButtonUp','setappdata(gcf,''bdown'',0);');
            m = 0;
            while getappdata(gcf,'bdown')
              if w<0 v = v * (1 - w/100)^cp;
              else   v = v + w*cp;
              end;
              if     v>pd{3} if m break; else v=pd{2}; end;
              elseif v<pd{2} if m break; else v=pd{3}; end;
              end;
              m = 1;  pd{1} = v;  s = ['  ' Pftoa(pd{8},v) '  '];
              set(Ret1,'str',s,'user',pd);
              setappdata(gcf,'OBJ',Ret1); evalRep(pd{4},{'@VAL' s '@OBJ' 'plt("misc",3)'});
              q=p; while q>0 & getappdata(gcf,'bdown') pause(.01); q=q-.01; end;
              p = rpt;
            end;
            return;
          else c=321;
          end;
        end;
        if c==321
          t = get(Ret1,'str');
          if ~strcmp(t,'error') pd{10} = t; end;
          t = deblank(t);  while(t(1)==' ') t(1)=[]; end;
          oldc = get(Ret1,cp);
          if ~ui pd{11} = get(Ret1,'interp'); set(Ret1,'interp','none',ERAS,ERAXOR); end;
          set(Ret1,'str',[t '_'],cp,oldc([2 3 1]),'user',pd);
          set(gcf,'keypress',{@plt 'edit' Ret1 'click' 1});
        end;
      end;
    case 320,    switch sum(lower(argv))
                  case {0,541}, Ret1 = pd{1};
                  case 650,    Ret1 = [pd{2:3}];
                  case 617,    Ret1 = pd{4};
                  case 615,    Ret1 = pd{5};
                  case 428,      Ret1 = pd{6};
                  case 642,    Ret1 = pd{7};
                  case 649,    Ret1 = pd{8};
                  case 512,     Ret1 = pd{9};
                  otherwise,       Ret1 = pd(1:9);
                  end;
    case {541,323,650},
                  if     length(argv)==3  pd(2:3) = {argv(2) argv(3)};
                  elseif length(argv)==2  pd(2:3) = {argv(1) argv(2)}; set(Ret1,'user',pd);  return;
                  end;
                  pd{1}  = argv(1);  s = ['  ' Pftoa(pd{8},argv(1)) '  '];
                  set(Ret1,'str',s,'user',pd);
                  if sargn==541 setappdata(gcf,'OBJ',Ret1); evalRep(pd{4},{'@VAL',s,'@OBJ','plt("misc",3)'}); end;
    case 617, pd{4} = argv;  set(Ret1,'user',pd);
    case 428,   pd{6}   = argv;  set(Ret1,'user',pd);
    case 642, pd{7}    = argv;  set(Ret1,'user',pd);
    case 615, pd{5}    = argv;  set(Ret1,'user',pd);
    case 649, if     isnumeric(argv) argv = sprintf('%%%dw',argv);
                  elseif length(argv)<2  argv = ['%' argv 'w']; end;
                  pd{8}    = argv;  set(Ret1,'user',pd);
    case 512,  if isnumeric(argv)
                       m = argv;  argv = getappdata(m,'argv');
                  else m = 0;
                  end;
                  if ischar(argv) argv = {argv ''}; end;
                  p = get(Ret1,'pos');  u = get(Ret1,'units');  hz = 'right';
                  s = argv{1};  q = argv{2};
                  if ui
                    r = round(22 + 6.5*length(s));
                    if u(1)=='n' sz=get(gcf,'pos'); r=r/sz(3); end;
                    switch length(q)
                    case 0, q = [p(1)-r p(2) r p(4)];
                    case 1, q = [p(1)-q p(2) q p(4)];
                    case 2, w = q(1);  if ~w w = r; end;
                            xo = real(q(2));  yo = imag(q(2));
                            q = [p(1)+xo-w p(2)+yo w p(4)];
                            if ~xo  q(1) = p(1)+(p(3)-w)/2;  hz = 'cent';  end;
                    case 4, hz='cent';
                    end;
                    if m set(m,'pos',q); return; end;
                    a = uicontrol('sty','text','str',s,'units',u,'pos',q,'foregr',[.6 .6 .8],'backgr',get(Ret1,'backgr'),...
                             'horiz',hz,'fontsi',9,'buttond',get(Ret1,'buttond'),'ena','inact','user',Ret1,'tag','E');
                  else
                    p = p(1:2);
                    switch length(q)
                    case 0, r = pd{8};  t = length(get(Ret1,'str')) - 4;
                            r = max([t length(Pftoa(r,pd{2})) length(Pftoa(r,pd{3}))]);
                            r = 5 + r*4;
                            if u(1)=='n' sz=get(gcf,'pos'); r=r/sz(3); end;
                            q = p - [r 0];
                    case 1, xo=real(q); yo=imag(q); q=p+[xo yo];
                            if ~xo hz='cent'; end;
                    case {2,3}, hz='cent';
                    end;
                    if m set(m,'pos',q); return; end;
                    a = text(0,0,s,'units',u,'pos',q,'horiz',hz,'color',[.6 .6 .8],...
                             'buttond',get(Ret1,'buttond'),'user',q(1:2)-p);
                  end;
                  j=3;  while length(argv)>j set(a,argv{j},argv{j+1}); j=j+2; end;
                  pd{9} = a;  set(Ret1,'user',pd);
                  setappdata(a,'argv',argv(1:2));
    otherwise,    set(Ret1,argn,argv);
                  a = pd{9};
                  if ishandle(a) & argn(1)=='p' plt('edit',Ret1,'label',a); end;
    end;
  end;

case 643
  Hin = varargin{2};
  if length(Hin)>1
     if Narg>2  In2 = varargin{3}; else In2 = 100;   end;
     if Narg>3  In3 = varargin{4}; else In3 = '';    end;
     if Narg>4  In4 = varargin{5}; else In4 = '';    end;
     In5 = 1;  In6 = '%2w %6w %3w';  In7 = [.75 .75 .75; 0 1 1];
     for k=6:8
       if Narg<k break; end;
       v = varargin{k};
       if isnumeric(v)
          if length(v)<3 In5 = v; else In7 = v; end;
       else In6 = v;
       end;
     end;
     if isempty(In2)                      In2 = 100; end;
     if length(In2)==1                    In2 = [0 In2]; end;
     if length(In2)==2 | length(In2)==4   In2 = [mean(In2(1:2)) In2]; end;
     if length(In2)==3                    In2 = [In2 -1e99 1e99]; end;
     if length(In7(:,1))<3                In7 = [In7; 0 0 0; 0 0 0]; end; 
     if ischar(In6) & length(In6(:,1))==1
        d = find(In6 == ' ');
        if length(d)<2 In6 = ['%2w ' In6 ' %3w'];  d = find(In6 == ' '); end;
        In6 = {In6(1:d(1)-1) In6(d(1)+1:d(2)-1) In6(d(2)+1:end)};
     end;
     if iscell(In6)
       for k=1:3  d=In6{k}; if length(d)<2 In6{k} = ['%' d 'w']; end; end;
       In6 = char(In6);
     end;
     hs = zeros(1,8);
     hs(2)  = uicontrol('sty','text','str',In3,'horiz','cent');  if isempty(In3) set(hs(2),'vis','off'); end;
     Ret1       = hs(2);        cb = {@plt 'slider' Ret1};
     hs(4) = uicontrol('sty','text','str','a','user',In6,'horiz','left');
     hs(5) = uicontrol('sty','text','str','a','user',In4,'horiz','right');
     hs(3)  = uicontrol('sty','slid','user',In2(1),'call',[cb {'CBslide'}]);
     hs(6)  = uicontrol('sty','edit','backg',In7(2,:),'foreg',In7(4,:),'user',In2(4:5),'horiz','cent','call',[cb {'CBedit'}]);
     set(hs(2:5),'backg',In7(1,:),'foreg',In7(3,:));
     In5 = [In5 10];  hs([7 8]) = In5(1:2);
     if hs(7)==5 & In2(4)<0 In2(4) = 1e-99; end;
     set(hs(2),'user',hs);
     plt('slider',Ret1,'set','minmax',In2(2:5),In2(1));
     plt('slider',Ret1,'set','pos',Hin);
     setappdata(gcf,'sli',[getappdata(gcf,'sli') hs(2:6)]);
     return;
  end;

  hs   = get(Hin,'user');
  fmat = get(hs(4),'user');

  Action = varargin{3};
  if Narg>3  In1 = varargin{4};  else In1 = ''; end;
  if Narg>4  In2 = varargin{5};  else In2 = []; end;
  if Narg>5  In3 = varargin{6};  else In3 = []; end;

  switch sum(Action)
  case 662
     oldval = get(hs(3),'user');    newval = get(hs(3),'Val');  sval = [];
     smax   = get(hs(3),'max');   smin   = get(hs(3),'min');
     dv = newval-oldval;  sdv = sign(dv);
     switch hs(7)
     case 4,
        v = hs(8)*round(newval/hs(8));
        if v==oldval v = v + hs(8)*sdv; end;
        newval = v;
     case 2,
        v = round(newval);
        if v==oldval v = v + sdv; end;
        newval = v;
     case 5, sval = newval; newval = exp(sval);
     case 3,
        if newval >= get(hs(3),'user')-eps  newval=2^nextpow2(newval);
        else                                  newval=2^nextpow2(newval/2);
        end;
     case 6,
        pdv = abs(dv)/(smax-smin);
        if     pdv < .02   newval = oldval + hs(8)*sdv;
        elseif pdv < .11   newval = oldval + hs(8)*sdv*10;
        end;
        newval = hs(8)*round(newval/hs(8));
     end;
     if isempty(sval)
       if newval<smin newval = smin; elseif newval>smax newval = smax; end;
       sval = newval;
     end;
     set(hs(6),'str',Pftoa(fmat(2,:),newval));
     set(hs(3),'Val',sval,'user',newval);
     evalRep(get(hs(5),'user'),{'@VAL',sprintf('%g',newval)});

  case 555
     newval = s2d(get(hs(6),'str'));    sval = [];
     minmax = get(hs(6),'user');
     if isempty(newval)
         set(hs(6),'str',Pftoa(fmat(2,:),get(hs(3),'Val')));
     else
        switch hs(7)
        case 2,  newval = round(newval);
        case 3,  newval = 2 ^ nextpow2(newval/1.414);
        case 4, newval=hs(8)*round(newval/hs(8));
        case 5,     newval = min(max(newval,minmax(1)),minmax(2));
                       sval = log(newval);
        end;
        newval = min(max(newval,minmax(1)),minmax(2));
        Ret1 = newval;
        set(hs(6),'str',Pftoa(fmat(2,:),newval));
        if isempty(sval) sval = newval; end;
        sval = max(min(get(hs(3),'max'),sval),get(hs(3),'min'));
        if hs(7)==5 newval = exp(sval); end;
        set(hs(3),'Val',sval,'user',newval);
        evalRep(get(hs(5),'user'),{'@VAL',sprintf('%g',newval)});
     end;

  case 320
     switch sum(In1)
     case 750,         Ret1 = sum(get(hs(2),'vis')) == 221;
     case 308,             Ret1 = sum(get(hs(2),'ena')) == 221;
     case 315,             Ret1 = hs(2:6);
     case {338 885}, Ret1 = get(hs(2),'pos');  Ret1 = Ret1(1:3);
     otherwise Ret1 = s2d(get(hs(6),'str'));
     end;
  case 332
     if isnumeric(In1) In2 = In1;  In1 = 'value'; end;
     switch sum(In1)
     case 541
        set(hs(6),'str',num2str(In2));
        Ret1 = plt('slider',Hin,'CBedit');
     case 323
        set(hs(6),'str',num2str(In2));
        sv = get(hs(5),'user');  set(hs(5),'user','');
        Ret1 = plt('slider',Hin,'CBedit');
        set(hs(5),'user',sv);
     case 495,  set(hs(2:6),'vis','on')
                   if isempty(get(hs(2),'str')) set(hs(2),'vis','of'); end;
     case 557, set(hs(2:6),'vis','of');
     case 527, set(hs([3 6]),'ena','of');
     case 465,  set(hs([3 6]),'ena','on');
     case {338 885},
        a = (get(0,'screenp')-96)/10;
        p = get(gcf,'pos');  w = p(3);  h = p(4);
        if length(In2<3) In2 = [In2 120]; end;
        lmax = 7 + 7*fix(s2d(fmat(3,2:(length(deblank(fmat(3,:)))-1))));
        lmin = 7 + 7*fix(s2d(fmat(1,2:(length(deblank(fmat(1,:)))-1))));
        if isnan(lmax) lmax = 28; end;
        if isnan(lmin) lmin = 14; end;
        h0=a+17;  h1=16+a/2;  h2=17-a-1;  h3=a+17;  h5=h3+2;
        if In2(1)>=1
          if In2(3)<1 In2(3) = In2(3)*w; end;
          In2(3) = max(In2(3),65);
          ya = In2(2) - h0;
          set(hs(2:6),'unit','pix',{'pos'},...
                 {[In2(1:3) h1];
                  [In2(1),In2(2)-33,In2(3),h2];
                  [In2(1),ya,lmin,h3];
                  [In2(3)+In2(1)-lmax,ya,lmax,h3];
                  [In2(1)+lmin,ya-1,In2(3)-(lmin+lmax),h5]});
        else
          if In2(3)>=1 In2(3) = In2(3)/w; end;
          In2(3) = max(In2(3),65/w);
          h1=h1/h; h2=h2/h; h3=h3/h; h5=h5/h; ya = In2(2)-h0/h;  lmin=lmin/w;  lmax=lmax/w;
          set(hs(2:6),'unit','nor',{'pos'},...
                 {[In2(1:3) h1];
                  [In2(1),In2(2)-33/h,In2(3),h2];
                  [In2(1),ya,lmin,h3];
                  [In2(3)+In2(1)-lmax,ya,lmax,h3];
                  [In2(1)+lmin,ya-1/w,In2(3)-(lmin+lmax),h5]});
         end;
     case 421
        hs(7) = In2(1);
        if length(In2)>1 hs(8) = In2(2); end;
        set(hs(2),'user',hs);
        Ret1 = plt('slider',Hin,'CBedit');
     case 650
        lg = hs(7)==5;
        if lg mn=1e-99; else mn=-1e+99; end;
        if length(In2)<4 In2 = [In2 mn 1e99]; end;
        set(hs(5),'str',[Pftoa(fmat(3,:),In2(2)) ' ']);
        set(hs(4),'str',[' ' Pftoa(fmat(1,:),In2(1))]);
        set(hs(6),'user',In2(3:4));
        if isempty(In3) v = s2d(get(hs(6),'str'));
        else            v = In3;
        end;
        v = min(max(v,In2(3)),In2(4));
        set(hs(6),'str',Pftoa(fmat(2,:),v));
        m1 = In2(1);  m2 = In2(2);
        if lg m1=log(m1); m2=log(m2); v=log(v); end;
        set(hs(3),'min',m1,'max',m2,'Val',min(max(m1,v),m2));
        Ret1 = s2d(get(hs(6),'str'));
     case 512,  set(hs(2),'str',In2);
                   if isempty(In2) set(hs(2),'vis','of'); end;
     end;
  otherwise disp(sprintf('Invalid Action in plt(''slider''), Action=%s',Action));
  end;

case 390
  ax = varargin{2};
  if Narg==2 Action='update'; else  Action = varargin{3};  end;
  if sum(Action) == 436
    if Narg<4 grdclr = [.13 .13 .13]; else grdclr = varargin{4}; end;
    if Narg<5 erMode = ERAXOR;  else erMode = varargin{5}; end;
    if Narg<6 sty = '-';        else sty    = varargin{6}; end;
    axes(ax);  f = findobj(ax,'type','line');
    if length(f) f=get(f(1),'buttond'); end;
    gridH = line('color',grdclr,'tag','SkipCur','buttond',f,'linest',sty,'user','grid',ERAS,erMode);
    setappdata(gridH,'clr',grdclr);
    Action = 'on';
    if Narg >= 7 set(gridH,varargin{7},varargin{8}); end;
    if Narg >= 9 set(gridH,varargin{9},varargin{10}); end;
    set(ax,'child',flipud(get(ax,'child')));
  else  gridH = findobj(ax,'user','grid');
        if length(gridH) ~= 1 return; end;
        z = get(gridH,'z');
        Ret1 = length(z);
        if Ret1 Ret1 = z(1) | Ret1>=2; end;
  end;

  switch sum(Action)
  case 643, if Ret1 Action='on'; end;
                fixMark;
  case 642, if Ret1 Action='off'; else Action='on'; end;
  end;

  x=0; y=0;
  switch sum(Action)
  case 221,  z = 1;  Ret1 = 1;  l = [get(ax,'xlim') get(ax,'ylim')];
             s = get(ax,'Xscale');  t = 1;
             if s(2)=='o' &  l(2)/l(1)<getappdata(gcf,'logTR')
               t = get(ax,'XtickLabel'); if iscell(t) t = [t{:}]; else t = t(:)'; end;
               t = length(findstr(t,'.'));
             end;
             if t  xt = [l(1) get(ax,'XTICK') l(2)]; n=length(xt);
                   if xt(2) <= xt(1)   xt = xt(2:n); n=n-1; end;
                   if xt(n) <= xt(n-1) n=n-1;  xt = xt(1:n); end;
             else  xt = logTicks(l(1:2));  n = length(xt);
             end;
             s = get(ax,'Yscale');  t = 1;
             if s(2)=='o' & l(4)/l(3)<getappdata(gcf,'logTR')
               t = get(ax,'YtickLabel'); if iscell(t) t = [t{:}]; else t = t(:)'; end;
               t = length(findstr(t,'.'));
             end;
             if t  yt = [l(3) get(ax,'YTICK') l(4)]; m=length(yt);
                   if yt(2) <= yt(1)   yt = yt(2:m); m=m-1; end;
                   if yt(m) <= yt(m-1) m=m-1;  yt = yt(1:m); end;
             else  yt = logTicks(l(3:4));  m = length(yt);
             end;
             t = m + n - 4;
             if t s = ones(3,1);  y = [1 m m];  x = [1 n n];  n=n-1;  m=m-1;
                  z = [0;0;NaN] * ones(1,t);           z=z(:);
                  y = [s*yt(2:m)  yt(y)'*ones(1,n-1)]; y=y(:);
                  x = [xt(x)'*ones(1,m-1) s*xt(2:n)];  x=x(:);
             end;
  case 315, z=0; Ret1=0;
  end;
  set(gridH,'z',z,'y',y,'x',x);

case 436
  Action = varargin{2};
  if Narg > 2  In1 = varargin{3};
               if Narg > 3  In2 = varargin{4}; end;
  end;
  switch sum(Action)
  case 436
     Hhcpy(12)= In1;
     if Narg ==4  renderer = In2; else renderer = '-painters';  end;
     Hhcpy(1)=figure('menu','none','numberT','off','back','off','resize','off','pos',auxLoc(350,105),'color',[0,.4,.4],'name','Hardcopy','tag',get(gcf,'tag'));
     foobar = 0;
     if exist('foobar')
          cpath = feval('which','plt.m');
     else cpath = GetExe;
     end;
     filen = [fileparts(cpath) filesep 'pltHcpy.mat'];
     if exist(filen) load(filen);
     else
       z = [0 0];  t = z+1;  w = [0 1];
       pS1=[z;w;[2 1];w;t;w;t;z;t;w;w];
       pS2 = 'temp.bmp';
     end;
     Hhcpy(3)   = uicontrol('sty','pop','str',['Win Meta|Bit Map|HPGL|LaserJet IIp|Post Script|Encaps PS|Windows'],...
                             'call','plt hcpy ModePUcb;');
     Hhcpy(5)  = uicontrol('sty','radio','str','color','call','plt hcpy rbg1 C;');
     Hhcpy(4)     = uicontrol('sty','radio','str','B&W','call','plt hcpy rbg1 BW;');
     Hhcpy(6) = uicontrol('sty','radio','str','Clip Board','call','plt hcpy rbg2 cb;');
     Hhcpy(8)    = uicontrol('sty','radio','str','Device','call','plt hcpy rbg2 dev;');
     Hhcpy(7)   = uicontrol('sty','radio','str','File','call','plt hcpy rbg2 file;');
     Hhcpy(11) = uicontrol('sty','text','str',pS2,'horiz','left');
     Hhcpy(10)  = uicontrol('sty','text','str', 'Path/File:','horiz','left');
     Hhcpy(9)  = uicontrol('str','Print','call','plt hcpy print;','user',renderer);
     Hhcpy(2)   = uicontrol('str','Select File','call','plt hcpy OpenPBcb;');
     cncl            = uicontrol('str','Cancel','call','close(gcf)');
     h1 = Hhcpy([3 5 4 6 8 7]);
     h2 = Hhcpy([11 10]);
     h3 = [Hhcpy([9 2]) cncl];
     set(h1,{'Val'},{pS1(3,1); pS1(5,1); pS1(4,1); pS1(6,1); pS1(8,1); pS1(7,1)});
     set([h1 h2],'backg',[.8,.8,.9],'foreg','black');
     set([h1 h2 h3],{'pos'},{[115 80 100 15]; [ 10 85 65 15]; [10 70  65 15]; [245 85 95 15];
                           [245 55  95 15]; [245 70 95 15]; [10 35 330 15]; [ 10 50 65 15];
                           [110  5  55 20]; [ 10  5 85 20]; [180 5  55 20]});
     for i=1+1:12-1
       if pS1(i,2)==1  set(Hhcpy(i),'ena','on'); else set(Hhcpy(i),'ena','of'); end;
     end;
     if get(Hhcpy(8),'Val') set(Hhcpy(2),'str','Select Dev');
     else                       set(Hhcpy(2),'str','Select File');
     end;
  case 364
     switch sum(In1)
     case 67,  set(Hhcpy(5),'Val',1); set(Hhcpy(4),'Val',0);
     case 153, set(Hhcpy(4),'Val',1);    set(Hhcpy(5),'Val',0);
     end;
  case 365
     switch sum(In1)
     case 319,  set(Hhcpy(6),'Val',0);  set(Hhcpy(7),'Val',0);
                 set(Hhcpy(9),'ena','on'); set(Hhcpy(11),'ena','of');
                 set(Hhcpy(8),'Val',1);      set(Hhcpy(2),'str','Select Dev','ena','on');
     case 197,   set(Hhcpy(6),'Val',1);   set(Hhcpy(7),'Val',0);
                 set(Hhcpy(9),'ena','on'); set(Hhcpy(11),'ena','of');
                 set(Hhcpy(8),'Val',0);     set(Hhcpy(2),'ena','of');
     case 416, set(Hhcpy(6),'Val',0);  set(Hhcpy(7),'Val',1);
                 set(Hhcpy(8),'Val',0);     set(Hhcpy(2),'str','Select File','ena','on');
                 fileN = get(Hhcpy(11),'str');
                 if length(fileN)<5 fileN='none'; set(Hhcpy(11),'str',fileN); end;
                 if strcmp(fileN(1:4),'none') set(Hhcpy(9),'ena','of'); end;
                 set(Hhcpy(11),'ena','on');
     end;
  case 751
     s = get(Hhcpy(3),'Val');
     en = [1 1 0 1 1; 1 1 0 1 0; 1 1 0 1 1; 1 1 0 0 1; 0 1 0 1 1; 0 1 0 1 1; 0 0 1 1 1];
     va = [0 1 0 1 0; 0 1 0 1 0; 0 1 0 0 1; 0 1 0 0 1; 0 1 0 1 0; 0 1 0 1 0; 0 0 1 1 0];
     ext= ['wmf';     'bmp';     'hgl';     'jet';     'ps ';     'eps';     'xxx']; ext=ext(s,:);
     if s==7  set(Hhcpy(2),'str','Select Dev', 'ena','on'); set(Hhcpy(11),'ena','of');
     else     set(Hhcpy(2),'str','Select File','ena','on'); set(Hhcpy(11),'ena','on');
              fileN = get(Hhcpy(11),'str');  lf = length(fileN);
              if lf < 5 fileN='none'; set(Hhcpy(11),'str',fileN); end;
              d = findstr('.',fileN);
              if length(d) set(Hhcpy(11),'str',[fileN(1:d(end)) ext]);
              else         set(Hhcpy(11),'str',[fileN '.' ext]);
              end;
     end;
     ena = {'off' 'on'};
     set(Hhcpy(6),'ena',ena{1+en(s,1)},'Val',va(s,1));
     set(Hhcpy(7),  'ena',ena{1+en(s,2)},'Val',va(s,2));
     set(Hhcpy(8),   'ena',ena{1+en(s,3)},'Val',va(s,3));
     set(Hhcpy(5), 'ena',ena{1+en(s,4)},'Val',va(s,4));
     set(Hhcpy(4),    'ena',ena{1+en(s,5)},'Val',va(s,5));
  case 745
     if strcmp('Select File',get(Hhcpy(2),'str'))
        ext=['wmf';'bmp';'hgl';'jet';'ps ';'eps';'xxx']; ext=ext(get(Hhcpy(3),'Val'),:);
        [fileN,pathN]= uiputfile(['*.',ext],'Open (new) File');
        if fileN
           fi = fileN;
           if isempty(findstr('.',fi)) fi = [fi '.' ext]; end;
           set(Hhcpy(11),'str',[pathN fi]);  set(Hhcpy(9),'ena','on');
        end;
     elseif strcmp('Select Dev',get(Hhcpy(2),'str')) print -dsetup
     else   disp('Invalid OpenPBcb string');
     end;
  case 557
     invertStr= get(Hhcpy(12),'invert');
     set(Hhcpy(12),'invert','off');
     figure(Hhcpy(12)); drawnow;
     PrintMode = get(Hhcpy(3),'Val');
     colorFlg  = get(Hhcpy(5),'Val');
     printStr = sprintf('print -f%d',Hhcpy(12));
     renderer = get(Hhcpy(9),'user');
     options={[' ' renderer ' -dmeta '];' -dbitmap ';' -dhpgl ';' -dljet2p '};
     if colorFlg options = [options; {' -dpsc ';' -depsc ';' -dwinc '}];
     else        options = [options; {' -dps ' ;' -deps  ';' -dwin  '}];
     end;
     printStr=[printStr options{PrintMode,:} ' '];
     if get(Hhcpy(7),'Val')==1
        PathFilen = get(Hhcpy(11),'str');  printStr=[printStr '''' PathFilen ''''];
     end;
     printStr = [printStr ' -noui'];  axh = [];  axh = findobj(Hhcpy(12),'type','axes');  nax = length(axh);
     if ~colorFlg
        figC=get(Hhcpy(12),'color');
        axChild = [];   axCol = zeros(nax,3);  axSpcl = zeros(nax,1);
        for i=1:nax
           x = get(axh(i),'color');
           if strcmp('none',x) axSpcl(i) = 1; else axCol(i,:)  = x; end;
           axChild = [axChild; get(axh(i),'child')];
           xCol(i,:) = get(axh(i),'xcol');  yCol(i,:) = get(axh(i),'ycol');
           tCol(i,:) = get(get(axh(i),'title'),'color');
           if xCol(i,:) == figC  axSpcl(i) = axSpcl(i) + 2; end;
        end;
        naxChild = length(axChild);
        for i=1:naxChild
           if sum(get(axChild(i),'type'))==528   kidCol(i,:) = get(axChild(i),'facecolor');
                                                    set(axChild(i),'facecolor',[.25 .25 .25]);
           else kidCol(i,:) = get(axChild(i),'color'); set(axChild(i),'color','black');
           end;
        end;
        set(Hhcpy(12),'color','white');
        for i= 1:nax
           if axSpcl(i)==1 | axSpcl(i)==3
           else set(axh(i),'color','white');
           end;
           if axSpcl(i)>=2 set(axh(i),'xcol','white','ycol','white' );
           else            set(axh(i),'xcol','black','ycol','black' );
           end;
        end;
        for i=1:nax
           if get(get(axh(i),'title'),'color')==figC set(get(axh(i),'title'),'color','white');
           else                                    set(get(axh(i),'title'),'color','black');
           end;
           if axSpcl(i) >=2
                set(get(axh(i),'xlab'),'color','white'); set(get(axh(i),'ylab'),'color','white');
           else set(get(axh(i),'xlab'),'color','black'); set(get(axh(i),'ylab'),'color','black');
           end;
        end;
     end;
     set(Hhcpy(1),'pointer','watch'); drawnow; pause(1);
     if PrintMode == 2
        set(Hhcpy(1),'vis','of');
        refresh
        eval(printStr);
        set(Hhcpy(1),'vis','on');
     else  drawnow discard;  eval(printStr);
     end;
     if ~colorFlg & (PrintMode ~=4 | PrintMode ~=7)
        set(Hhcpy(12),'color',figC);
        for i=1:nax
           if ~rem(axSpcl(i),2) set(axh(i),'color',axCol(i,:)); end;
           set(axh(i),'xcol',xCol(i,:)); set(axh(i),'ycol',yCol(i,:));
           set(get(axh(i),'title'),'color',tCol(i,:));
        end;
        for i=1:naxChild
           if sum(get(axChild(i),'type'))==528 set(axChild(i),'facecolor',kidCol(i,:));
           else                                   set(axChild(i),'color',kidCol(i,:));
           end;
        end;
        if PrintMode==2 drawnow; else drawnow discard; end;
     end;
     pS1=zeros(12-1,2);
     for i=1+1:12-1
        pS1(i,1) = get(Hhcpy(i),'Val');
        if sum(get(Hhcpy(i),'ena'))==221 pS1(i,2)=1; end;
     end;
     pS2=get(Hhcpy(11),'str');
     foobar = 0;
     if exist('foobar')
          cpath = feval('which','plt.m');
     else cpath = GetExe;
     end;
     save([fileparts(cpath) filesep 'pltHcpy.mat'],'pS1','pS2');
     close(Hhcpy(1));
     set(Hhcpy(12),'invert',invertStr);
     clear Hhcpy
  otherwise disp([Action ' invalid Action in plt(hcpy)']);
  end;

case 670

    if Narg<3 disp('Not enough arguments plt(''cursor'',...)'); return;  end;
    CurID = varargin{2};
    Action = varargin{3};
    In1 = '';  In2 = '';  In3 = '';  In4 = '';  In5 = '';
    In6 = '';  In7 = '';  In8 = '';  In9 = '';
    if Narg > 3
      In1 = varargin{4};
      if Narg > 4
        In2 = varargin{5};
        if Narg > 5
          In3 = varargin{6};  if iscell(In3) In3 = char(In3); end;
          if Narg > 6
            In4 = varargin{7};
            if Narg > 7
              In5 = varargin{8};
              if Narg > 8
                In6 = varargin{9};  if iscell(In6) In6 = char(In6); end;
                if Narg > 9
                  In7 = varargin{10};
                  if Narg > 10
                    In8 = varargin{11};
                    if Narg > 11 In9 = varargin{12}; end;
                  end;
                end;
              end;
            end;
          end;
        end;
      end;
    end;

    CurMain = getappdata(0,'CurMain');
    sact = sum(Action);
    if sact ~= 436
      if isempty(CurID) return; end;
      Hc = get(CurMain(CurID),'user');
      ax = Hc(14);  ax2 = Hc(15);
      if isempty(get(ax2,'par')) ax2 = []; end;
      axlink = 0;
      if ishandle(ax2)
           if sum(get(ax2,'vis')) == 315  ax2=[];
           else axlink = [get(get(ax2,'ylab'),'str') '  '];
                axlink = axlink(1)~=92 | axlink(2) ~= 'd';
           end;
      end;
      hix = Hc(4);  hiy = Hc(6);
      misc = get(hix,'user');
      actv = misc(4);
      iact = 15 + actv;
      ix   = misc(1);
      hact = Hc(iact);
      lh = get(hact,'user');
      lk = lh{1};
      Ret1 = hact;  Ret2 = lk;
      if ~ishandle(hact) return; end;
      xylim = [get(ax,'xlim') get(ax,'ylim')];
      xyli = xylim;
      if get(hact,'par')==ax2 xyli(3:4) = get(ax2,'ylim'); end;
      DualCur = getappdata(ax,'DualCur');
    end;
    Narg2=Narg;
    switch sact
    case 436
      if Mver < 8.4 Hc = zeros(1,15);
      else          Hc = gobjects(1,15);
      end;
      ax = CurID(1);
      DualCur = getappdata(ax,'DualCur');
      if isempty(DualCur) DualCur = 0;  setappdata(ax,'DualCur',0); end;
      Hc(14) = ax;
      if length(CurID)==2  ax2 = CurID(2);  ax2a = ax2; else ax2 = [];  ax2a = 0; end;
      Hc(15) = ax2a;
      hf = get(ax,'par');
      if isempty(In7) vis='on'; else vis=In7; end;
      Cn = length(CurMain);
      f = find(~CurMain);
      if length(f) CurID = f(1);
      else CurMain = [CurMain 0]; CurID = Cn+1;
      end;
      if CurID>200 disp('Warning: Clear actions missing for 200 cursor inits'); end;
      Ret1 = CurID;
      CurIDstr ={@plt 'cursor' CurID};
      hLines = getappdata(ax,'Lhandles');
      skp = get(hLines,'tag');
      if ~iscell(skp) skp = {skp}; end;
      hLines = hLines(cellfun('length',skp) ~= 7);
      axes(ax);
      set(gcf,'vis','of');
      if isempty(In2) In2 = [.7 .7 .7; 0 0 0; 1 1 .5; 1 0 0; 0 0 0]; end;
      trk = ~sum(In2(2,:));
      if length(hLines)
         In2n = length(In2(:,1));
         clines = length(hLines);
         cli4 = clines + 4;
         if In2n < cli4
           In2 = In2(min(1:cli4,In2n*ones(1,cli4)),:);
         end;
         In4n = length(In4);
         if In4n<clines
           In4 = In4(min(1:clines,In4n*ones(1,clines)));
         end;
         In2n     = length(In2(:,1));
         objColor = In2(5:In2n,:);
         tact = 0;
         if length(ax2) hLines2 = findobj(ax2,'type','line');
         else           hLines2 = 0;
         end;
         if Mver < 8.4  Hc = [Hc zeros(1,clines)];
         else           Hc = [Hc gobjects(1,clines)];
         end;
         for i=1:clines
            hi = hLines(i);
            xy = [get(hi,'x'); get(hi,'y')];
            if sum(objColor(i,:)) curColor = objColor(i,:);
            else                  curColor = get(hi,'color');
            end;
            Hc(15+i)=line('x',xy(1,1),'y',xy(2,1),ERAS,ERAXOR,'color',curColor,'linestyle','none',...
                            'clipping','on','vis',vis,'user',{hi trk});
            set(Hc(15+i),'marker',In4(i),'MarkerSize',In5);
            set(hi,'buttond',[CurIDstr {'lineCB' i}]);
            if ~tact & (sum(get(hi,'vis'))==221) tact = i; end;
            if length(find(hi==hLines2)) set(Hc(15+i),'par',ax2); end;
         end;
      else disp('no lines to attach cursors to')
      end;
      if max(max(In1))>2 unitt = 'Pixels'; else unitt = 'Normal'; end;
      fontsz = (196-get(0,'screenpix'))/10;
      for i=2:3
        if i==2  cbStr  = 'x'; else cbStr  = 'y'; end;
        cbStr= [CurIDstr {'scaleAxis' cbStr}];
        if isempty(In3) In3 = ['x';'y']; end;
        Hc(i) = uicontrol(hf,'Style','text','fontsi',fontsz,'vis',vis,'str',deblank(In3(i-2+1,:)),...
                'Units',unitt,'pos',In1(i-1,:),'backg',In2(1,:),'horiz','cent','buttond',cbStr,'ena','inact');
      end;
      for i=4:7
        if trk==1 objColor = [1 1 1]; else objColor = In2(2,:); end;
        cbStr  = [CurIDstr {'editCB' i-4}];
        if i==5 | i==7   objColor=get(gcf,'color'); end;
        Hc(i)  = uicontrol(hf,'Style','edit','fontsi',fontsz,'vis',vis,'str',' ','Units',unitt,'pos',In1(i-1,:),...
                           'foreg',get(ax,'color'),'backg',objColor,'horiz','left','call',cbStr);
      end;
      set(Hc([5 7]),'vis','of');
      hix = Hc(4);  hiy = Hc(6);
      ih = [8 9 10 11];  c = [173 175 79 68];  u = {-inf inf '' ''};
      cbs = {{'peakval' 0} {'peakval' 1} {'mlsCB' ''} {'markCB' ''}};
      for k = 1:4
        m = ih(k);
        Hc(m) = uicontrol(hf,'Units',unitt,'pos',In1(m-1,:),'vis',vis,'horiz','cent',...
                  'str',char(c(k)),'fontname','symbol','fontw','bol',...
                  'call',[CurIDstr cbs{k}],'user',u{k});
      end;
      set(Hc(8), 'ui',uicontextmenu('call',[CurIDstr {'peakval' 2}]));
      set(Hc(9),'ui',uicontextmenu('call',[CurIDstr {'peakval' 3}]));
      s = ['if length(findobj(gcf,''tag'',''DeltaC'',''vis'',''on'')) ' ...
           'delete([findobj(gcf,''tag'',''mark''); findobj(gcf,''tag'',''markR'')]); ' ...
           'h = findobj(gcf,''str'',''D''); eval(get(h(end),''callb'')); ' ...
           'else h=findobj(gcf,''str'',''O''); h=h(end); u=1-get(h,''user''); set(h,''user'',u); ' ...
           'set(getappdata(gcf,''Lhandles''),''LineSmooth'',char(''of''+[0 8*u])); end;'];
      set(Hc(10),'ui',uicontextmenu('call',s));
      set(Hc(11),'fontsi',12,'ui',uicontextmenu('call', ...
         ['if length(findobj(gcf,''tag'',''DeltaC'',''vis'',''on'')) ' ...
          'h = findobj(gcf,''str'',''D''); eval(get(h(end),''callb'')); plt move res;' ...
          'else plt move; end;']));
      Hc(13) = line('x',[],'y',[],ERAS,ERAXOR,'color',In2(4,:),'vis','of');
      set(Hc(13),'marker','+','MarkerSize',5*In5,'tag','DeltaC');
      if isempty(In6) In6 = ['%7w';'%7w']; end;
      set(Hc(3),'user',In6);
      if length(In1(:,1)) >= 11
        uicontrol(hf,'style','slider','Units',unitt,'pos',In1(11,:),...
             'backg',[.3 .3 .3],'Min',0,'Max',1000,'Val',500,'user',500,...
             'vis',vis,'call',[CurIDstr {'xslider'}],'tag','xslider');
      end;
      set(get(ax,'xlab'),'buttond',[CurIDstr {'xincr'}]);
      i = [CurIDstr {'axisCB'}];
      Hc(12) = line('x',[],'y',[],ERAS,ERAXOR,'vis','of','color',In2(3,:),'buttond',i);
      set([hf ax],'buttond',i);
      if isempty(In8) monoFlag = -1; else monoFlag = In8; end;
      set(Hc(5),'user',In9);
      misc = zeros(1,8);
      if length(ax2) set(ax2,'buttond',[CurIDstr {'AxisCBaux'}]); end;
      misc(5) = monoFlag;
      misc(4) = max(1,tact);
      set(hix,'user',misc);
      set(Hc(7),'user',''); setappdata(Hc(7),'CB2','');
      set(Hc(2),'user',Hc);
      CurMain(CurID) = Hc(2);
      setappdata(0,'CurMain',CurMain);

    case {643 715}
      set(Hc(15+1:end),'vis','of');
      if isempty(In1) | ~In1 In1 = ix; end;
      if sum(get(lk,'vis')) == 315
        for i=15+1:length(Hc)
          lh = get(Hc(i),'user'); lk = lh{1};
          if sum(get(lk,'vis')) == 221 break; end;
        end;
        actv = i-15; misc(4) = actv;
        iact = 15 + actv;
        hact = Hc(iact);
        if ~ishandle(hact) return; end;
      end;
      if sum(get(lk,'vis')) == 315 return; end;
      if isempty(In2)
        x = get(lk,'x');  y = get(lk,'y');  n = length(x);
        if ~n return; end;
        if In1>n | In1<=0  In1 = round(n/2); end;
        In2 = x(In1);  In3 = y(In1);
      end;
      a = get(hact,'par');
      if length(getappdata(a,'subTr')) evalQ(get(Hc(5),'user')); return; end;
      if a == ax  xyl = xylim;  else xyl = xyli; end;
      xlim = xyl(1:2);  ylim = xyl(3:4);
      xd = 1.01 * (In2-xlim);  yd = 1.01 * (In3-ylim);
      if xd(1)<=0 xlim = xlim + xd(1); elseif xd(2)>=0 xlim = xlim+xd(2); end;
      if yd(1)<=0 ylim = ylim + yd(1); elseif yd(2)>=0 ylim = ylim+yd(2); end;
      if ~isequal(xyl,[xlim ylim]) ...
           & sum(Action)==643 ...
           & isempty(getappdata(a,'hold'))
         setlim(a,'xlim',xlim); setlim(a,'ylim',ylim);
         if a==ax plt('grid',ax); end;
         evalQ(get(Hc(5),'user'));
      end;
      set(hact,'y',In3,'x',In2,'vis','on');
      multi = getappdata(gcf,'multi');
      if length(multi)
        Lha = getappdata(gcf,'Lhandles'); n = length(Lha);
        for k = 1:n  yd = get(Lha(k),'y'); yd = yd(min(In1,length(yd)));
                     set(multi(k),'pos',[In2 yd],'str',prin(get(multi(k),'tag'),yd));
                     set(multi(k+n),'x',In2,'y',yd);
        end;
        set(multi(end),'x',[In2 In2]);
      end;
      bf = {'backg' 'foreg'};
      fmt = get(Hc(3),'user');
      lc = get(lk,'color');
      fr = repmat(max(lc.*[.9 1 .5])<=.5,1,3);
      set([hix hiy],{'str'},{Pftoa(fmt(1,:),In2); Pftoa(fmt(2,:),In3)},bf,{lc fr});
      dmode = sum(get(Hc(13),'vis'))==221;
      if dmode
         dxy = {Pftoa(fmt(1,:),In2-get(Hc(13),'x')); Pftoa(fmt(2,:),In3-get(Hc(13),'y'))};
         set(Hc([5 7]),{'str'},dxy,bf,get(hix,bf),'vis','on');
      end;
      if DualCur
        if DualCur>0 dact=15+DualCur; ow=0;
        else dact = iact-DualCur;
             ow = (dact > length(Hc));
             if ow dact=iact+DualCur; end;
        end;
        Hcd = Hc(dact);
        dh = get(Hcd,'user'); dk = dh{1};
        if sum(get(dk,'vis')) == 221
          xy = [get(dk,'x'); get(dk,'y')];
          misc(6) = xy(2,In1);
          xylen = length(xy(1,:)); if ix>xylen ix=xylen; end;
          yd = xy(2,In1);
          a = get(Hcd,'par'); b = get(a,'YaxisLoc');
          if b(1)=='r' b = get(a,'ylim'); yd = max(b(1),min(b(2),yd)); end;
          set(Hcd,'x', max(xyli(1),min(xyli(2),xy(1,In1))),'y',yd,'vis','on');
          if ~dmode set(Hc(5),'vis','of');
                    s = Pftoa(fmt(2,:),xy(2,In1));
                    if ow t = s;  s = get(Hc(6),'str'); set(Hc(6),'str',t); end;
                    dc = get(dk,'color');
                    fr = repmat(max(dc.*[.9 1 .5])<=.5,1,3);
                    set(Hc(7),bf,{dc fr},'str',s,'vis','on');
          end;
        end;
      end;
      misc([1 2 3]) = [In1 In2 In3];  set(hix,'user',misc);
      s = getappdata(Hc(7),'CB2');  evalQ(s);
      sx = getappdata(ax,'xstr');  sy = getappdata(ax,'ystr');  s = [sx sy];
      rep1 = {'@XVAL','real(@XY)',...
              '@YVAL','imag(@XY)',...
              '@XY',  'plt("cursor",@CID,"get","position")',...
              '@IDX', 'plt2nd({"cursor",@CID,"get","position"})',...
              '@LNUM','plt("cursor",@CID,"get","activeLine")',...
              '@HAND','plt2nd({"cursor",@CID,"get","activeLine"})',...
              '@XU'  ,'plt("misc",1)',...
              '@YU'  ,'plt("misc",2)',...
              '@CID', sprintf('%d',CurID)};
      for k=1:length(s)
        set(s(k),'str',evalRep2(getappdata(s(k),'evl'),rep1));
      end;
      if length(getappdata(ax,'moveCBext')) | length(find(ax==findobj(gcf,'type','axes')))
        s = get(Hc(7),'user');
        if length(s) evalRep(s,rep1); end;
      end;

    case 472
      p0 = getappdata(get(ax,'xlab'),'OldCur');
      sc = get(ax,'yscale');  liny = sc(2)=='i';
      if ischar(In1)
          aPt = get(ax2,'currentp');  cy = aPt(1,2);
          y = get(ax2,'ylim');
          if liny  y = y + p0(2) - cy;
          else     y = y * p0(2) / cy;
          end;
          set(ax2,'ylim',y);
      else
        aPt = get(ax,'currentp');
        cx = aPt(1,1);  cy = aPt(1,2);
        if In1<1 x = xylim(1:2);
                 sc = get(ax,'xscale');
                 if sc(2)=='i' x = x + p0(1) - cx;
                 else          x = x * p0(1) / cx;
                 end;
                 if length(ax2) axb = [ax ax2]; else axb = ax; end;
                 set(axb,'xlim',x);
        end;
        if In1   y = xylim(3:4);
                 if liny  yi = p0(2) - cy;  y = y + yi;
                 else     yi = p0(2) / cy;  y = y * yi;
                 end;
                 set(ax,'ylim',y);
                 if axlink
                   yr = get(ax2,'ylim');
                   if liny  yr = yr + yi * diff(yr)/diff(y);
                   else     yr = yr * (yr(2)/yr(1))^(log(yi)/log(y(2)/y(1)));
                   end;
                   set(ax2,'ylim',yr);
                 end;
        end;
      end;
      plt('grid',ax);
      fixMark;
    case 606
      p0 = getappdata(get(ax,'xlab'),'OldCur');
      if ischar(In1)
          aPt = get(ax2,'currentp');
          y = get(ax2,'ylim');  y0 = y(1);  y1 = y(2);
          cy = aPt(1,2)-y0;  if ~cy cy=1e-06; end;
          set(ax2,'ylim',[y0 y0 + abs((p0(2)-y0)*(y1-y0)/cy)]);
      else
        aPt = get(ax,'currentp');
        if In1<1 x0 = xylim(1);  x1 = xylim(2);
                 cx = aPt(1,1)-x0;  if ~cx cx=1e-06; end;
                 if length(ax2) axb = [ax ax2]; else axb = ax; end;
                 set(axb,'xlim',[x0 x0 + abs((p0(1)-x0)*(x1-x0)/cx)]);
        end;
        if In1   y0 = xylim(3);  y1 = xylim(4);  dy = y1-y0;
                 cy = aPt(1,2)-y0;  if ~cy cy=1e-06; end;
                 cy = dy/cy;  p2 = p0(2)-y0;
                 set(ax,'ylim',[y0 y0 + abs(p2*cy)]);
                 if axlink
                    sc = get(ax,'yscale');
                    if sc(2)=='i'
                       y = get(ax2,'ylim');  y0 = y(1);  y1 = y(2);
                       p2 = (y1-y0) * p2 / dy;
                       set(ax2,'ylim',[y0 y0 + abs(p2*cy)]);
                    else
                       yr = get(ax2,'ylim');  yr0 = yr(1);  yr1 = yr(2);
                       aPt = get(ax2,'currentp');
                       cy = aPt(1,2)-yr0;  if ~cy cy=1e-06; end;
                       p0 = yr0 * exp(log(p0(2)/y0) * log(yr1/yr0) / log(y1/y0));
                       set(ax2,'ylim',[yr0 yr0 + abs((p0-yr0)*(yr1-yr0)/cy)]);
                     end;
                 end;
        end;
      end;
      plt('grid',ax);
      fixMark;
    case 662
      set(hact,'vis','of');
      if DualCur set(Hc(15+DualCur),'vis','of'); end;
      aPt  = get(ax,'currentp');
      xCur = max(xylim(1),min(xylim(2),aPt(1,1)));  yCur = max(xylim(3),min(xylim(4),aPt(1,2)));
      aPt   = get(Hc(12),'user');
      xOld  = aPt(1);  yOld  = aPt(2);
      if sum(get(ax,'yscale')) == 322
            dxyOk = abs(log(xylim(3)/xylim(4))) < 50 * abs(log(yOld/yCur));
      else  dxyOk = diff(xylim(3:4)) < 50 * abs(yOld-yCur);
      end;
      if dxyOk
         if sum(get(ax,'xscale')) == 322
               dxyOk = abs(log((.0001*xylim(2)+xylim(1))/xylim(2))) < 50 * abs(log(xOld/xCur));
         else  dxyOk = diff(xylim(1:2)) < 50 * abs(xOld-xCur);
         end;
      end;
      if dxyOk
          set(Hc(12),'x',[xOld xCur xCur xOld xOld],'y',[yOld yOld yCur yCur yOld]);
          fmt=get(Hc(3),'user');
          if sum(get(Hc(12),'vis')) == 315
             set(Hc(12),'vis','on');
             set(hiy,'str',Pftoa(fmt(2,:),yOld));
             set(hix,'str',Pftoa(fmt(1,:),xOld));
             set(Hc([5 7 4 6]),'ena','on');
             set(Hc(11),'ena','of');
          end;
          bk = get(hix,'backg');
          fr = get(hix,'foreg');
          set(Hc(7),'str',Pftoa(fmt(2,:),yCur),'backg',bk,'foreg',fr,'vis','on');
          set(Hc(5),'str',Pftoa(fmt(1,:),xCur),'backg',bk,'foreg',fr,'vis','on');
      end;
      MotionZoom = getappdata(gcf,'MotionZoom');
      if length(MotionZoom) feval(MotionZoom,0,[xOld xCur yOld yCur]); end;
    case {570, 557},
      b = get(hact,'buttond'); if length(b) eval(b); return; end;
      clkType = sum(get(gcf,'SelectionT'));
      if (sum(get(Hc(12),'vis'))==315) & ~misc(7) & ~misc(8)
         CurIDstr = {@plt 'cursor' CurID};
         smv = '';
         StopMotion = 'set(gcf,''WindowButtonMotionFcn'','''',''WindowButtonUpFcn'','''');';
         aPt  = get(ax,'currentp');  cx = aPt(1,1);  cy  = aPt(1,2);
         if     cx<xylim(1) | cx>xylim(2)   dragy =  1;
         elseif cy<xylim(3) | cy>xylim(4)   dragy =  0;
         else                               dragy = -1;
         end;
         xCur = max(xylim(1),min(xylim(2),cx));
         dxy = diff(xylim);
         switch clkType
         case 649
            yCur = max(xylim(3),min(xylim(4),cy));
            ovt = 0;
            if sact==570
               Narg2 = 4;
               if sum(get(lk,'vis')) ~= 221
                  for i=15+1:length(Hc)
                      lh = get(Hc(i),'user');  lk = lh{1};
                      if sum(get(lk,'vis')) == 221 break; end;
                  end;
                  actv = i-15;
               end;
               In1 = actv;
               if length(ax2)
                 ylm = get(ax2,'ylim');  dylm = diff(ylm);
                 cvert = (cy - xylim(3)) / dxy(3);
                 mdst = 1e+99;
                 idst = 0;
                 for i=15+1:length(Hc)
                   if get(Hc(i),'par') ~= ax2 continue; end;
                   lh = get(Hc(i),'user'); lk = lh{1};
                   if sum(get(lk,'vis')) == 315 continue; end;
                   x = get(lk,'x'); y = get(lk,'y');
                   if length(x)>999 | all(diff(x)>0)
                     [toss,j] = min(abs(xCur-x));
                     dvert = (y(j) - ylm(1)) / dylm;
                     acv = abs(cvert-dvert);
                   else
                     acv = min(abs((x-xCur)/dxy(1)) + abs((y-ylm(1))/dylm - cvert));
                   end;
                   if acv < mdst  idst = i;  mdst = acv; end;
                 end;
                 if mdst < .02
                   ovt = 1;
                   ii = idst-15;
                   if ii ~= DualCur In1 = ii; end;
                 end;
               end;
            end;
            lh = get(Hc(15+In1),'user'); lk = lh{1};
            x = get(lk,'x'); y = get(lk,'y');
            if Narg2==4
               set(Hc(8), 'user',-inf);
               set(Hc(9),'user',inf);
               if misc(5)==-1
                  df = diff(x);
                  monoFlag = all(df>0) | all(df<0);
               else monoFlag = misc(5);
               end;
            else monoFlag=In2;
            end;
            if monoFlag [junk,imin] = min(abs(xCur-x));
            elseif get(lk,'par')==ax2
                  ylm = get(ax2,'ylim');  dylm = diff(ylm);
                  cvert = (cy - xylim(3)) / dxy(3);
                  [toss,imin] = min(abs((x-xCur)/dxy(1)) + abs((y-ylm(1))/dylm - cvert));
            else  [toss,imin] = min(abs((x-xCur)/dxy(1)) + abs(yCur-y)/dxy(3));
            end;
            actv = In1;  misc(4) = actv;
            iact = 15 + actv;
            hact = Hc(iact);
            set(hix,'user',misc);
            plt('cursor',CurID,'update',imin,x(imin),y(imin));
            if (sact==557) | ovt
              if Narg2==4
                smv = [CurIDstr {'lineCB' In1 monoFlag 0}];
              end;
            else setappdata(get(ax,'xlab'),'OldCur',[cx cy xylim]);
                 StopMotion = [CurIDstr {'svHist'}];
                 smv = [CurIDstr {'panAX' dragy}];
            end;
         case 321
            setappdata(get(ax,'xlab'),'OldCur',[cx cy xylim]);
            StopMotion = [CurIDstr {'svHist'}];
            smv = [CurIDstr {'zoomAX' dragy}];
         otherwise,
            set(Hc(12),'user',[max(xylim(1),min(xylim(2),cx)) max(xylim(3),min(xylim(4),cy))]);
            smv = [CurIDstr {'expbox'}];
         end;
         if length(smv) set(gcf,'WindowButtonMotionFcn',smv,'WindowButtonUpFcn',StopMotion); end;
      elseif (sum(get(Hc(12),'vis'))==221) | misc(7) | misc(8)
        if clkType==649  plt('cursor',CurID,'scale','new'); end;
        plt('cursor',CurID,'restore');
        misc([7 8]) = 0;  set(hix,'user',misc);
      end;
    case 872
      CurIDstr = {@plt 'cursor' CurID};
      aPt  = get(ax2,'currentp');  setappdata(get(ax,'xlab'),'OldCur',[aPt(1,1) aPt(1,2)]);
      if sum(get(gcf,'SelectionT'))==321 smv = [CurIDstr {'zoomAX' 'R'}];
      else          smv = [CurIDstr {'panAX'  'R'}];
      end;
      set(gcf,'WindowButtonMotionFcn',smv,'WindowButtonUpFcn','set(gcf,''WindowButtonMotionFcn'','''',''WindowButtonUpFcn'','''');');
    case 555
      x   = str2double(get(Hc(4+In1),'str'));
      fmt = get(Hc(3),'user'); fmt = fmt((1+fix(In1/2)),:);
      xpx = (sum(get(Hc(12),'vis'))==221) | misc(7) | misc(8);
      if xpx  xy = Hc(12); xy = [get(xy,'x') get(xy,'y')];
              if In1 xold=xy(3*In1); else xold=xy(1); end;
      else    xy = [get(hact,'x') get(hact,'y')];
              if ~In1 xold=xy(1); elseif In1==2 xold=xy(2); else xold=0; end;
      end;
      if length(x)
         if xpx  ixy = [1 4 5; 2 3 3; 6 7 10; 8 9 9];  xy(ixy(In1+1,:)) = x;
                 set(Hc(12),'x',xy(1:5),'y',xy(6:10));
         else
            editd = length(get(hact,'buttond'));
            if ~In1
               xd = get(lk,'x');  yd = get(lk,'y');
               if editd set(hact,'x',x); plt click EDIT 5 0;
               else [v ix] = min(abs(xd-x));
                    xv = xd(ix);
                    if xv<xyli(1) | xv>xyli(2) plt('cursor',CurID,'set','xlim',xv+[-.5 .5].*diff(xyli(1:2))); end;
                    plt('cursor',CurID,'update',ix);
               end;
            elseif In1==2  set(hact,'y',x);
                           if x<xyli(3) | x>xyli(4)
                              plt('cursor',CurID,'set','ylim',x+[-.5 .5].*diff(xyli(3:4)));
                           end;
                           if editd plt click EDIT 5 0; end;
            else set(Hc(4+In1),'str',Pftoa(fmt,misc(6)));
            end;
         end;
      else set(Hc(4+In1),'str',Pftoa(fmt,xold));
      end;
    case 641
      set(gcf,'WindowButtonMotionFcn','','WindowButtonUpFcn','');
      set(Hc(12),'x',[xylim(1:2) 0 0 0],...
                      'y',[xylim(3) 0 xylim(4) 0 0]);
      p0 = getappdata(get(ax,'xlab'),'OldCur');
      set(ax,'xlim',p0(3:4),'ylim',p0(5:6));
      plt('cursor',CurID,'scale','new',0,0);
    case 520
      expHis = get(hiy,'user');
      if isempty(expHis) expHis = [xylim 1]; end;
      lExp   = length(expHis(:,1));
      curExp = find(expHis(:,5)==1);
      skip   = [];
      switch sum(In1)
      case 330
         xy = Hc(12); xy = [get(xy,'x') get(xy,'y')];
         expLim = [sort(xy(1:2)) sort(xy([6 8])) 1];
         if isequal(xylim,expLim(1:4))
            if length(curExp) expHis(curExp,5)=1; end;
         else
            if length(curExp)
               expHis(curExp,5)=0;
               if max(curExp,lExp) < 4
                  expHis = [expHis; expLim];
               else
                  if curExp==4  expHis(2:4,:)=[expHis(3:4,:); expLim];
                  elseif curExp >= 1  expHis(curExp+1,:)=expLim;
                  end;
               end;
            else
               if lExp < 4  expHis = [expHis; expLim];
               else               expHis = [expHis(lExp-1,:); expLim];
               end;
            end;
         end;
         autoScale=0;
      case 319
          if sum(get(gcf,'SelectionT'))==321 autoScale=1;
          else
            if isempty(curExp) curExp = lExp;
            else               expHis(curExp,5)=0;  curExp = curExp-1;
            end;
            if curExp  autoScale=0;  expLim=expHis(curExp,:); expHis(curExp,5)=1;
            else       autoScale=1;
            end;
          end;
      case 441
          switch sum(In2)
          case 120,    autoScale=2;
          case 121,    autoScale=3;
          case 429, autoScale=1;
          otherwise, disp([In2 ' is not a valid action in plt(cursor,CurID,scale,auto,In2)']);
          end;
      otherwise, disp([In1 ' is not a valid action in plt(cursor,CurID,scale,In1)']);
      end;
      if autoScale == 0
          setlim(ax,'xlim',sort(expLim(1:2)));
          setlim(ax,'ylim',sort(expLim(3:4)));
          wii = 0;
      else
         numLines=length(Hc)-15;  minx=+inf;  maxx=-inf;
         temp = zeros(numLines,2);  lineList = [];   wii = 1;
         for i=1:numLines
            hLine=get(Hc(15+i),'user');
            if sum(get(hLine{1},'vis')) == 221
               xdata = get(hLine{1},'x');
               minx  = min(minx,min(xdata));
               maxx  = max(maxx,max(xdata));
               if wii ydata = get(hLine{1},'y'); end;
               lineList = [i lineList];  wii = 0;
            end;
         end;

         if wii msgbox('Possible autoscale without any visible lines','Warning','warn','modal');
         else
            for i=1:numLines
              k = Hc(15+i);
              tx = get(k,'x');  ty = get(k,'y');
              temp(i,:) = [tx(1) ty(1)];
              set(k,'x',xdata(1),'y',ydata(1));
            end;
            skip = findobj(ax,'tag','SkipCur','vis','on'); set(skip,'vis','of');
            switch autoScale
            case 1,  set(ax,'YlimM','auto');
                         expLim = [minx maxx get(ax,'ylim') 1];
                         setlim(ax,'xlim',expLim(1:2));
            case 2, expLim = [minx maxx get(ax,'ylim') 1];
                         setlim(ax,'xlim',expLim(1:2));
            case 3,
               misc=get(hix,'user');
               if ismember(actv,lineList)
                  ydata = get(lk,'y');
                  ymax = max(ydata);  ymin = min(ydata);  dy = 0.25*(ymax-ymin);
                  ymin = ymin - dy;   ymax = ymax + dy;
                  if ymin ~= ymax  set(ax,'YlimM','man'); set(ax,'ylim',[ymin ymax]); drawnow;
                  else             set(ax,'YlimM','auto');  drawnow;
                  end;
               else set(ax,'YlimM','man'); drawnow; set(ax,'YlimM','auto');
               end;
               expLim = [xylim 1];
            end;
            for i=1:numLines set(Hc(15+i),'x',temp(i,1),'y',temp(i,2)); end;
         end;
      end;
      if ~wii
         if length(ax2) setlim(ax2,'xlim',expLim(1:2)); end;
         set(hiy,'user',expHis);
         if axlink & Narg<6
           yr = get(ax2,'ylim');
           setlim(ax2,'ylim',yr(1) + (get(ax,'ylim') - xylim(3)) * diff(yr) / (xylim(4)-xylim(3)));
         end;
         axes(ax); evalQ(get(Hc(5),'user'));
      end;
      set(ax,'YlimM','man');
      set(skip,'vis','on');
      plt('grid',ax);
      xView = getappdata(gcf,'xView');
      if length(xView) set(xView{1},'x',get(ax,'xlim')); end;
    case 548,
      lx = length(get(lk,'x'));
      ty = sum(get(gcf,'SelectionT'));
      if ty==434 & getappdata(gcbo,'ty')==321 ty = 321; end;
      rpt = getappdata(gcf,'repeat');
      if length(rpt)>1 p = rpt(2); else p = .4; end;
      if length(rpt) rpt=rpt(1); else rpt = .03; end;
      setappdata(gcf,'bdown',1);
      set(gcf,'WindowButtonUp','setappdata(gcf,''bdown'',0);');
      while getappdata(gcf,'bdown')
        if ty==321 ix = max(ix-1,1);
        else        ix = min(ix+1,lx);
        end;
        setappdata(gcbo,'ty',ty);
        plt('cursor',CurID,'update',ix);
        pause(p);  p = rpt;
      end;
    case 763,
      x = get(lk,'x');   y = get(lk,'y'); lx = length(x);
      xlim = xylim(1:2);   ylim = xylim(3:4);  dxlim = diff(xlim);
      v = get(gcbo,'Val');  x1 = min(x);  x2 = max(x);   xrange = x2 - x1;
      if xrange/dxlim < 2
        dx = round(v - 500);
        if abs(dx)==10 pmove = .01; else pmove = .05; end;
        ix = max(min(lx,ix + sign(dx)*round(lx*pmove)),1);
        x = x(ix);  xd = 1.01 * (x-xlim);
        if xd(1)<0 xlim = xlim + xd(1); elseif xd(2)>0 xlim = xlim+xd(2); end;
        v = 500;
      else
        dx = v - get(gcbo,'user');
        switch round(abs(dx))
          case 10,   xlim = xlim + sign(dx)*dxlim/10;
          case 100,  xlim = xlim + sign(dx)*dxlim;
          otherwise, xlim = x1 + xrange*v/1000 + dxlim*[-.5 .5];
        end;
        if     xlim(1)<x1  xlim = [x1 x1+dxlim] - dxlim/50;
        elseif xlim(2)>x2  xlim = [x2-dxlim x2] + dxlim/50;
        end;
        xlc = mean(xlim);
        v = 1000*(xlc-x1)/xrange;
        v = max(min(v,1000),0);
        [dmy ix] = min(abs(xlc-x));
      end;
      set(gcbo,'Val',v,'user',v);
      y = y(ix);  yd = 1.01 * (y-ylim);
      if yd(1)<0 ylim = ylim + yd(1); elseif yd(2)>0 ylim = ylim+yd(2); end;
      if sum(xylim - [xlim ylim])
         set(ax,'xlim',xlim,'ylim',ylim);
         if length(ax2) set(ax2,'xlim',xlim); end;
         plt('grid',ax);
         evalQ(get(Hc(5),'user'));
      end;
      plt('cursor',CurID,'update',ix);
    case 740
      switch In1
        case 2, In1=0; set(Hc(8), 'user',-inf);
        case 3, In1=1; set(Hc(9),'user',inf);
      end;
      y = get(lk,'y'); x = get(lk,'x');
      xx = find(x <= xylim(2) & x >= xylim(1));
      if isempty(xx) x = xylim(1); xx = 1:length(x); disp('You must select a line for the min/max finder'); end;
      y = y(xx);  ly = length(y);
      if In1
         ix=find(y < [y(2:ly) inf] & ...
                 y < [inf  y(1:ly-1)] & ...
                 y > get(Hc(9),'user'));
         if isempty(ix)
              [y ix] = min(y);
         else [y i]  = min(y(ix));
              ix = ix(i);
         end;
         set(Hc(9),'user',y);
      else
         ix=find(y > [y(2:ly) -inf] & ...
                 y > [-inf y(1:ly-1)] & ...
                 y < get(Hc(8),'user'));
         if isempty(ix)
              [y ix] = max(y);
         else [y i] = max(y(ix));
              ix = ix(i);
         end;
         set(Hc(8),'user',y);
      end;
      ix = ix + xx(1) - 1;
      plt('cursor',CurID,'update',ix);
    case 772
      fmt = get(Hc(3),'user');
      set(hiy,'str',Pftoa(fmt(2,:),misc(3)));
      set(hix,'str',Pftoa(fmt(1,:),misc(2)));
      set(Hc(12),'vis','of');
      if sum(get(Hc(13),'vis')) == 315
         set(Hc([5 7]),'vis','of','str','');
      end;
      set(hact,'y',max(xyli(3),min(xyli(4),misc(3))),...
               'x',max(xyli(1),min(xyli(2),misc(2))));
      set(Hc([2:3]),'ena','inact');  set(Hc(11),'ena','on');
    case 465
        l = getappdata(gcf,'Lhandles');
        MRK ='marker';  STY = 'linest';
        k = get(gcbo,'tag');
        if isempty(k) k='1';
                      setappdata(gcbo,'mrk',get(l,MRK));
                      setappdata(gcbo,'sty',get(l,STY));
        end;
        mrk = 'o';  sty = getappdata(gcbo,'sty');
        switch k
          case '1', sty = 'none';
          case '3', mrk = getappdata(gcbo,'mrk'); k = '0';
                    
        end;
        set(gcbo,'tag',char(k+1));
        if iscell(mrk) MRK = {MRK}; end;
        if iscell(sty) STY = {STY}; end;
        set(l,MRK,mrk,STY,sty);
    case 560
      set(Hc(11),'backg',1-get(Hc(11),'backg'));
      if sum(get(Hc(13),'vis')) == 315
            setappdata(ax,'DualCsv',DualCur);
            setappdata(ax,'DualCur',0);
            set(Hc(13),'vis','on','par',get(hact,'par'),'user',actv,...
                'x',get(hact,'x'),'y',get(hact,'y'));
      else  set(Hc([13 5 7]),'vis','of');
            set(lk,'vis','of'); drawnow; set(lk,'vis','on');
            DualCur = getappdata(ax,'DualCsv');
            if DualCur setappdata(ax,'DualCur',DualCur); end;
      end;
    case 493
      fg = get(hix,'par');
      a = flipud(findobj(findobj(fg,'user','TraceID'),'type','text'));
      if length(a)<actv  d = 'Yval';
      else d = deblank(get(a(actv),'str'));
      end;
      set(findobj(fg,'user','idcur'),{'str','color'},{d,get(hact,'color')});
      bx = getappdata(fg,'bx');
      if length(bx) & ishandle(bx)
        ha = get(get(bx,'par'),'pos');  ha=ha(4);
        nb = round(.03839*ha - 2.27);
        s = get(bx,'user');  t = s{2}; s = s{1};
        d = d(1:min(9,length(d))); f = findstr(d,t);
        if isempty(f) v=''; else v = [blanks(f(1)-1) repmat('v',1,length(d))]; end;
        set(bx,'str',[s(1:ix-1) {t v} s(ix:end)],'Value',ix+1);
        drawnow;
        set(bx,'ListboxTop',max(1,ix-nb));
      end;
    case 320
      switch sum(In1)
      case 315
         Ret1 = [Hc(2:13) hact];
      case {338 885}
         if Narg2==5 lineNum=In2; else lineNum=actv; end;
         k = Hc(15+lineNum);
         Ret1 = complex(get(k,'x'),get(k,'y'));
         Ret2 = ix;
      case 1028
         Ret1 = actv;
         Ret2 = lk;
         if Narg>4 & sum(get(Ret2,'vis'))==315
           Ret1 = 0;
         end;  
      case 625
         expHis = get(hiy,'user');
         if isempty(expHis)
            Ret1  = [[xylim 1]; [zeros(3,4), -ones(3,1)]];
         else aa = 4 - length(expHis(:,1));
              Ret1 = [ expHis; [ zeros(aa,4), -ones(aa,1) ] ];
         end;
         if length(ax2) Ret1 = [Ret1; [get(ax2,'xlim') get(ax2,'ylim') 2]];
         else           Ret1 = [Ret1; [zeros(1,4) -1]];
         end;
      otherwise, disp(['Bad argument (' In1 ') in plt(cursor,get,...)']);
      end;
    case 332
      switch sum(In1)
      case 557
         h = [Hc([2:11 15+1:end]) ...
              findobj(gcf,'tag','xstr') ...
              findobj(gcf,'tag','ystr') ...
              findobj(gcf,'user','idcur') ...
              findobj(gcf,'tag','xslider')];
         h = findobj(h,'vis','on');
         set(h,'vis','of');
         setappdata(Hc(2),'hid',h);
      case 495, set(getappdata(Hc(2),'hid'),'vis','on');
      case {338 885}
         for i=2:11 set(Hc(i),'pos',In2(i-1,:)); end;
      case 423,     set(Hc(2),'str',In2);
      case 424,     set(Hc(3),'str',In2);
      case 570,   set(Hc(5),'user',In2);
      case 572,   if length(In2) & ischar(In2) & In2(1)==';'
                         setappdata(ax,'moveCBext',1);  In2(1) = [];
                      end;
                      set(Hc(7),'user',In2);
      case 622,  setappdata(Hc(7),'CB2',In2);
      case {563,443,442},
         oldEh = get(hiy,'user');
         expHis = oldEh;
         if length(expHis)
             curExp = find(expHis(:,5)==1);
             if length(curExp) expHis(curExp,5)=0; end;
         end;
         switch sum(In1)
         case 563, xyl = In2;  expHis(1,:) = [xyl 1];
         case 443
            if length(oldEh)  expHis(1,:) = [xylim(1:2) In2,1];  xyl = [xylim(1:2) In2];
            else              expHis=[];                         xyl = [xylim(1:2) In2];
            end;
         case 442
            if length(oldEh)  expHis(1,:) = [In2 xylim(3:4) 1]; xyl  = [In2 xylim(3:4)];
            else              expHis=[];                        xyl = [In2 xylim(3:4)];
            end
         end;
         set(ax,'xlim',sort(xyl(1:2)),'ylim',sort(xyl(3:4)));
         set(hiy,'user',expHis);
         if length(ax2)
            if Narg2==6 xy = In3; else  xy = get(ax2,'ylim'); end;
            set(ax2,'xlim',xyl(1:2),'ylim',xy);
            set(hix,'user',misc);
         end;
         plt('grid',ax,'update');
         evalQ(get(Hc(5),'user'));
      case 540, set(Hc(8), 'user',-inf); set(Hc(9),'user',inf);
      case 625,
         if length(In2(1,:)) < 5
            plt('cursor',CurID,'set','xylim',In2);
         elseif isequal(size(In2),[5 5])
           expHis = [];
           for i=1:5
             if i<=4
               if In2(i,5) >= 0
                 expHis=[expHis; In2(i,:)];
                 if In2(i,5) == 1
                    xylim = In2(i,1:4);
                    set(ax,'xlim',In2(i,1:2),'ylim',In2(i,3:4));
                    plt('grid',ax);
                 end;
               end
             else
               if In2(i,5)==2 set(ax2,'xlim',In2(i,1:2),'ylim',In2(i,3:4)); end;
             end;
           end;
           set(hiy,'user',expHis);
         else disp('error in plt(cursor,CurId,set,expHis,xxx), xxx is wrong shape');
         end;
      case 1028
        if In2
          misc(4) = In2;
          set(hix,'user',misc);
        end;
        if isempty(In3) In3 = -1; end;
        plt('cursor',CurID,'update',In3);
      otherwise, disp(['Invalid cursor set parameter: ' In1]);
      end;
    case 519
      delete(Hc([2:13 15+1:end]));
      CurMain(CurID) = 0;
      if ~sum(CurMain) CurMain = []; end;
      setappdata(0,'CurMain',CurMain);
    case 925
      ofs = strcmp(In1,'y');
      switch sum(get(gcf,'SelectionT'))
      case 649,
         if misc(7+ofs)
            plt('cursor',CurID,'scale','new');  misc(7)=0;  misc(8)=0;  set(hix,'user',misc);
         else
            fmt = get(Hc(3),'user'); fmt = fmt(ofs+1,:);
            of2 = ofs*2;
            set(Hc(4 +of2),'str',Pftoa(fmt,xylim(of2+1)),'ena','on' );
            set(Hc(5+of2),'str',Pftoa(fmt,xylim(of2+2)),'ena','on' ,...
                                 'backg',get(hix,'backg'),'foreg',get(hix,'foreg'),'vis','on');
            misc(7+ofs)=1;
            if ~misc(7) | ~misc(8)
               set(Hc(12),'x',xylim([1 2 2 1 1]),'y',xylim([3 3 4 4 3]),'vis','on');
            end;
            set(hix,'user',misc);
         end;
      case 321, plt('cursor',CurID,'scale','auto',In1);
      end;
    otherwise disp([Action ' is not a valid action in plt(cursor)']);
    end;

case 785
  switch s2i(varargin{2})
  case 0,
     p = gcbo;  e = getappdata(p,'edt');
     m = getappdata(e,'m'); obj = m{4}(1);
     c = get(p,'str');   v = get(p,'Val');  prop = deblank(c(v,:));
     f = getappdata(e,'f');  if length(f) close(f); end;
     if strcmp(prop,'Delete') & ishandle(obj) delete(obj); end;
     if ishandle(obj)
       s = get(obj,prop);
       if isnumeric(s)
          if length(s)==3 & max(s)<=1 & min(s)>=0 s = ctriple(s);
          else                                    s = num2str(s);
          end;
       end;
       set(e,'str',s);
     else set(e,'str','Deleted');
     end;
  case 1,
     e = gcbo;  p = getappdata(e,'pop');
     m = getappdata(e,'m'); obj = m{4};
     c = get(p,'str'); v = get(p,'Val');   prop = deblank(c(v,:));
     s = get(e,'str'); f = get(e,'user');
     if sum(get(gcf,'SelectionT'))==321 if strcmp(prop,'Color') plt('ColorPick'); end;
                   return;
     end;
     switch prop
     case 'Delete',
       if strcmpi(s,'all')
         if length(c(:,1))==9 c='line'; else c='text'; end;
         delete(findobj(get(get(e,'par'),'user'),'type',c,'tag','mark'));
       end;
     case 'Color',  plt('ColorPick');
     case {'Xdata','Ydata','Zdata'},
                s = str2num(s);  n = length(s); xy = char('X' + (prop(1)=='X'));
                for k=1:length(obj)
                  if strcmp(get(obj(k),'type'),'line') & length(get(obj(k),xy))==n set(obj(k),prop,s); end;
                end;
     otherwise, if isnumeric(get(obj(1),prop)) s = str2num(s); end;
                if size(get(p,'str'),1)<11 obj = findobj(obj,'type','line'); end;
                set(obj,prop,s);
     end;
   case 2,
     p = gcbo;  e = getappdata(p,'edt');  v = get(p,'Val');  f = get(e,'user');
     h = get(p,'user');  h = h{v};  h1 = h(1);
     pr = {'color' 'color' {'xcol' 'ycol'} 'ycol' 'color' 'color' 'color' 'linest'};
     pr = pr{v};
     if iscell(pr) s = get(h1,pr{1}); else s = get(h1,pr); end;
     if ischar(s) s = [s getappdata(h1,'er')];
     else         s = ctriple(s);
     end;
     set(e,'str',s);
     f = getappdata(e,'f');  if length(f) close(f); end;
     s = get(p,'str');
     setappdata(e,'m',{'str',e,pr,h,deblank(s(v,:))});
  case 3,
     s = get(gcbo,'str');
     v = get(getappdata(gcbo,'pop'),'Val');
     if v==8
       er = s(end);
       if length(s)>1
         if er=='n' s(end)=[]; elseif er=='x'; s(end)=[]; else er='x'; end;
       end;
       gr = findobj(get(gcf,'user'),'user','grid'); gr = gr(1);
       setappdata(gr,'er',er);  set(gr,'linest',s);
     else plt ColorPick;
     end;
  end;

case 428,
  in2 = varargin{2};
  if isnumeric(in2)
    switch in2
      case 1, Ret1 = get(findobj(gcf,'tag','xstr'),'user');
      case 2, Ret1 = get(findobj(gcf,'tag','ystr'),'user');
      case 3, Ret1 = getappdata(gcf,'OBJ');
      case 4, Ret1 = getappdata(gcf,'OBJ2');
    end;
    return;
  end;
  switch sum(in2)
    case 548,
      fs = get(gcf,'pos');  fs = fs(3:4);
      aPt = get(gcf,'currentp');
      switch Narg
      case 2
        t = gcbo;  ax = getappdata(gcf,'axis');   m = 1;
        if strcmp(get(t,'user'),'grid') t=ax(1); end;
        if     t == findobj(gcf,'user','TraceID') k=-1;
        elseif t == findobj(gcf,'tag', 'MenuBox') k=-2;
        else
          if Mver<7 k=find(~(ax-t)); else k=find(ax==t); end;
          if isempty(k)
             ty = {' xy' 'uic' 'pop' 'axi'};
             id = [  0    200   100   50  ];
             for m=2:4
               tw = ty{m}; k = getappdata(gcf,tw);
               if Mver<7 k=find(~(k-t)); else k=find(k==t); end;
               k = k + id(m);
               if length(k) break; end;
             end;
          end;
        end;
        assignin('base','hhh',t);
        if sum(get(gcf,'SelectionT'))==434 inspect(t); return; end;
        if m==3 ev=getappdata(t,'CBsv'); feval(ev{:});
                t = get(t,'par');
        elseif strcmp(get(t,'tag'),'E') t=get(t,'user');
        end;
        p = get(t,'pos');   r = get(t,'type');
        if r(3)=='p' | (r(3)=='c' & strcmp(get(t,'style'),'frame')) ...
                     | (r(1)=='a' & strcmp(get(t,'tag'),  'frame'))
           ch = {};
           r = p(1:2);  r = [r -r-p(3:4)];
           m = getappdata(gcf,'sli');  w = getappdata(gcf,'pop');  np = length(w);
           m = [w getappdata(gcf,'axi') getappdata(gcf,'uic') m(1:5:end)];
           for j=1:length(m)
             if j<= np q = getappdata(m(j),'ppos'); w = get(get(m(j),'par'),'pos');
             else      q = get(m(j),'pos');           w = q;
             end;
             s = q(1:2);  s = [s -s-q(3:4)];
             if all(s>r) ch = [ch {{m(j) w}}]; end;
           end;
           setappdata(t,'inside',ch);
        end;
        u = get(t,'units');  if u(1)=='n' aPt = aPt ./ fs; end;
        setappdata(gcf,'Dxy',{p aPt t k m});
        set(gcf,'WindowButtonMotionFcn','plt misc tidmv 0;',...
                'WindowButtonUpFcn',    'plt misc tidmv 0 0;');
      case 3
        s = getappdata(gcf,'Dxy'); t = s{3};
        u = get(t,'units');  u = u(1)=='n';
        if u  aPt = aPt ./ fs; end;
        r = getappdata(gcf,'snap');  r = r + 1e4*(r==0);
        if ~u  r = r ./ fs; end;
        rr = [r r];  r = [0 0 1./r];
        aPt = aPt-s{2}; s = s{1};
        if sum(get(gcf,'SelectionT'))==649
          aPt = [aPt 0 0];
          ch = getappdata(t,'inside');  cn = length(ch);
          for k = 1:cn
            h = ch{k}{1}; p = max(r,round((ch{k}{2}+aPt).*rr) ./ rr);
            ty = getappdata(h,'ty');
            switch ty(1)
              case 's', plt('slider',h,'set','pos',p);
              case 'p', plt('pop',h,'pos',p);
              otherwise set(h,'pos',p);
            end;
          end;
        else
          aPt = [0 0 aPt];
        end;
        s = s + aPt;
        s = max(r,round(s.*rr) ./ rr);
        tn = getappdata(t,'ty');
        if tn(1)=='e' plt('edit',t,'pos',s); else set(t,'pos',s); end;
      case 4
        set(gcf,'WindowButtonMotionFcn','','WindowButtonUpFcn','');
        s = getappdata(gcf,'Dxy');  t = s{3};  u = get(t,'units');
        p = get(t,'pos');  ty = getappdata(t,'ty');
        r = get(t,'type');
        if ty(1)=='p'    h = get(t,'user'); r = get(h,'str');
                         plt('pop',h,'index',plt('pop',h,'get'));
        elseif r(3)=='c' r = get(t,'str');
                         if iscell(r) r = r{1}; end;
        end;
        if u(1)=='n' a = strrep(prin('{ %0.3f}',p),'0.','.');
        else         a = prin('{ %4d}',round(p));
        end;
        prin(1,'%s: %3d %s;  %% %s\n',ty,s{4},a,r);
      end;
    case 579,
      switch Narg
      case 2
        t = gcbo;  ax = get(t,'par');  aPt = get(ax,'currentp');  aPt = aPt(1,1:2);
        tu = get(t,'units');
        set(t,'units','data');
        k = getappdata(gcf,'txt');  k = find(k==t);
        k = k+300;
        setappdata(gcf,'Dxy',{get(t,'pos') aPt t k tu});
        assignin('base','hhh',t);
        if sum(get(gcf,'SelectionT'))==434 inspect(t); return; end;
        set(gcf,'WindowButtonMotionFcn','plt misc txtmv 0;',...
                'WindowButtonUpFcn',    'plt misc txtmv 0 0;');
      case 3
        s = getappdata(gcf,'Dxy');
        t = s{3};  tp = s{1};  ax = get(t,'par');  aPt = get(ax,'currentp');

        set(t,'pos',tp(1:2)+aPt(1,1:2)-s{2});
      case 4
        set(gcf,'WindowButtonMotionFcn','','WindowButtonUpFcn','');
        s = getappdata(gcf,'Dxy'); t = s{3};
        set(t,'units',s{5}); p = get(t,'pos');
        ty = getappdata(t,'ty');     r = get(t,'str');  if iscell(r) r=r{1}; end;
        prin(1,'%s: %3d    %6V  %6V   ;  %% %s\n',ty,s{4},p(1:2),r);
        if ty(1)=='e' plt('edit',t,'pos',p); end;
      end;
    case 555,
      fs = get(gcf,'pos');  fs = fs(3:4);
      aPt = get(gcf,'currentp');
      switch Narg
      case 2
        w = getappdata(gcf,'sli');
        t = gcbo;
        k = 1 + 5*floor((find(w==t)+399)/5);
        if isempty(k) disp('slimv error'); return; end;
        t = w(k-400);
        assignin('base','hhh',t);
        if sum(get(gcf,'SelectionT'))==434 inspect(gcbo); return; end;
        u = get(t,'units');  if u(1)=='n' aPt = aPt ./ fs; end;
        setappdata(gcf,'Dxy',{get(t,'pos') aPt t k});
        set(gcf,'WindowButtonMotionFcn','plt misc slimv 0;',...
                'WindowButtonUpFcn',    'plt misc slimv 0 0;');
      case 3
        s = getappdata(gcf,'Dxy');  t = s{3};
        u = get(t,'units');  u = u(1)=='n';
        if u  aPt = aPt ./ fs; end;
        aPt = aPt-s{2}; s = s{1}(1:3);
        if sum(get(gcf,'SelectionT'))==649 s = s + [aPt 0];
        else                  s = max([0 0 .01],s+[0 0 aPt(1)]);
        end;
        r = getappdata(gcf,'snap');  r = r + 1e4*(r==0);
        if ~u q = get(gcf,'pos');  r = r ./ q(3:4); end;
        rr = [r r(1)];  r = [0 0 1/r(1)];
        s = max(r,round(s.*rr) ./ rr);
        plt('slider',t,'set','pos',s);
      case 4
        set(gcf,'WindowButtonMotionFcn','','WindowButtonUpFcn','');
        s = getappdata(gcf,'Dxy');  t = s{3};  u = get(t,'units');
        p = get(t,'pos');  r = get(t,'str');
        if u(1)=='n' disp(strrep(prin('sli: %3d { %0.3f}     ;  %% %s',s{4},p(1:3),r),'0.','.'));
        else         disp(prin('sli: %3d { %4d}     ;  %% %s',s{4},round(p(1:3)),r));
        end;
      end;
    case 642
      t = gco;
      aPt  = get(get(t,'par'),'currentp');
      if Narg==2
        switch sum(get(gcf,'SelectionT'))
        case 649
          setappdata(t,'Dxy',get(t,'pos')-aPt(1,:));
          set(gcf,'WindowButtonMotionFcn','plt misc marker 0;','WindowButtonUpFcn','set(gcf,''WindowButtonMotionFcn'','''',''WindowButtonUpFcn'','''');');
        case 321
          callp = 'plt MarkEdit 0;';
          calle = 'plt MarkEdit 1;';
          figure('menu','none','numberT','off','back','off','resize','off','pos',auxLoc(302,85),'color',[0,.4,.4],'name','Edit Marker','tag',get(gcf,'tag'),...
                         'closereq','plt click mark 4;');
          ed1 = uicontrol('sty','edit','pos',[8 5 130 22],'fontw','bol','call',calle,'buttond',calle);
          pu1 = 'Delete|Color|LineStyle|LineWidth|Marker|MarkerSize|Xdata|Ydata|Zdata';
          pu1 = uicontrol('sty','pop','str',pu1,'pos',[8 35 130 20],'call',callp,'Val',5);
          uicontrol('sty','text','str', 'Marker properties:','pos',[8 64 130 17]);
          ed2 = uicontrol('sty','edit','pos',[145 5 150 22],'fontw','bol','call',calle,'buttond',calle);
          pu2 = uicontrol('sty','pop','str','Delete|Color|FontAngle|FontName|FontSize|FontWeight|HorizontalAlign|Position|Rotation|String|VerticalAlign',...
                      'pos',[145 35 150 20],'call',callp,'Val',10);
          uicontrol('sty','text','str', 'String properties:','pos',[145 64 150 17]);
          set([ed1 ed2 pu1 pu2],'backg',[.8,.8,.9],'foreg','black');
          l = get(t,'user');
          setappdata(ed1,'pop',pu1); setappdata(pu1,'edt',ed1); setappdata(pu1,'obj',l);
          setappdata(ed2,'pop',pu2); setappdata(pu2,'edt',ed2); setappdata(pu2,'obj',t);
          setappdata(ed1,'m',{'str',ed1,'color',l,'Marker color'});
          setappdata(ed2,'m',{'str',ed2,'color',t,'String color'});
          if ishandle(l) l=get(l,'marker'); else l='Deleted'; end;
          if ishandle(t) t=get(t,'str');      else t='Deleted'; end;
          set(ed1,'str',l);  set(ed2,'str',t);
        end;
      else  dxy = aPt(1,:) + getappdata(t,'Dxy'); set(t,'pos',dxy(1:2));
      end;
    case 856, set(gcf,'SelectionT','alt');    plt('click',varargin{3:end});
    case 741,  set(gcf,'SelectionT','normal'); plt('click',varargin{3:end});
    case 660,   aid = findobj(gcf,'user','TraceID');
                    if length(aid)
                       uistack(aid,'top');
                    end;
    case 534,
      a = getappdata(gcf,'ucreq');
      if length(a) evalQ(a); end;
      creq = getappdata(gcf,'creq');
      if length(creq) eval(creq); end;
      plt('cursor',varargin{3},'clear');
      f = findobj('type','fig','tag',get(gcf,'tag'));
      for k = 1:length(f)
        fk = f(k);
        if fk ~= gcf
          c = get(fk,'closereq');
          if iscell(c) & length(c)==4 & strcmp(c{3},'close')
            plt('cursor',c{4},'clear');
          end;
          delete(fk);
        end;
      end;
      closereq;
  end;

case 902
  m = getappdata(gcbo,'m');
  if isempty(m) m = getappdata(get(gcbo,'par'),'m'); end;
  if iscell(m) hq=m; m=gcbo;
  elseif ishandle(m) hq = getappdata(m,'m');
  end;
  if Narg>1 ccf = varargin{2}; else ccf = ''; end;
  if strcmp(ccf,'C')
    for k=2:2:length(hq)
      h = hq{k};
      for j = 1:length(h) setappdata(h(j),'f',[]); end;
    end;
    setappdata(m,'f',[]);
    closereq;
    return;
  end;
  t = 0; cb = '';  f = getappdata(m,'f');
  if isempty(f)
    hb = [];  for k=2:2:length(hq) hb = [hb hq{k}]; end;
    if length(find(m==hb))
          b=m;
    else  b=hb(1); hb=[hb m];
    end;
    ed = 0;
    switch get(b,'type')
      case 'uicontrol', c = get(b,'backg');  t = get(b,'str');  ed = strcmp(get(b,'sty'),'edit');
      case 'text',      c = get(b,'color');  t = get(b,'str');
      case 'patch',     c = get(b,'facecolor');
      otherwise,        c = get(b,'color');
    end;
    if ed             c = [str2num(t) 0 0 0];   c = min(max(c(1:3),0),1);
    elseif ischar(t)  t = str2num(t);  if length(t)==3 & min(t)>=0 & max(t)<=1 c=t; ed=1; end;
    end;
    if Narg<3 & ed & sum(get(gcf,'SelectionT'))==649 & b==m
      t = ctriple(c);
      for k=2:2:length(hq)
        p = hq{k-1};  h = hq{k};
        if ~iscell(p) p = {p}; end;
        for j=1:length(p)
          q = p{j}; if strcmp(q(1:3),'str') set(h,q,t); else set(h,q,c); end;
        end;
      end;
      if length(ccf)>1 eval(ccf); end;
      return;
    end;
    f = figure('menu','none','numberT','off','back','off','pos',auxLoc(302,205,3),'color',[0,.4,.4],'name','Color Pick',...
               'tag',get(gcf,'tag'),'double','off','closereq','plt ColorPick C;');
    setappdata(f,'ccf',ccf);
    cb = 'plt ColorPick;';  s1 = [0 100 0 100];  d = round(c*100);
    p1 = [.02 .20 110];  p2 = [0 .29 0];
    s = [plt('slider',p1+2*p2,[d(1) s1],'Red (%)'  ,cb,2) ...
         plt('slider',p1+p2,  [d(2) s1],'Green (%)',cb,2) ...
         plt('slider',p1,     [d(3) s1],'Blue (%)' ,cb,2)];
    obj = plt('slider',s(1),'get','obj');  set(obj(5),'backg',[1 1 0]);
    ax = axes('xlim',[.8 12.15],'ylim',[.8 15.4],'color',[0 0 0],'xcol',[0,.4,.4],'ycol',[0,.4,.4],'XtickLabel',' ','YtickLabel',' ','TickLen',[0 0]',...
              'unit','nor','pos',[.408 .02 .57 .96]);
    if ischar(hq{end})
      text(-3.4,14.7,hq{end},'horiz','center','color',[1 .7 .8]);
    end;
    setappdata(f,'m',m);  setappdata(ax,'m',m);
    ph = zeros(11,11);
    for row = 1:11
      for col = 1:11
        ph(row,col) = patch(col+[0 1 1 0],row+[0 0 1 1],[0 0 0],'buttond',cb);
      end;
    end;
    c = d/100;
    pat = patch([1 12 12 1],12+[.5 .5 3 3],c,'buttond',cb);  setappdata(pat,'m',m);
    setappdata(f,'h',{c s pat ph ax});
    for k=1:length(hb) setappdata(hb(k),'f',f); setappdata(hb(k),'m',hq); end;
  end;
  h = getappdata(f,'h');
  p = get(gcbo,'par');
  s = h{2};
  if length(cb) k=1;
  else
    if strcmp(get(gcbo,'type'),'patch') & isempty(cb)
      if gcbo==h{3} c = h{1}; else c = get(gcbo,'facecolor'); end;
      for k=1:3
        obj = plt('slider',s(k),'get','obj');
        bk = get(obj(5),'backg');
        if bk(1) break; end;
      end;
    else
      y = strcmp(get(gcbo,'sty'),'edit');
      if y & isequal(get(gcbo,'user'),[0 100]) y=0; end;
      if y  c = [str2num(get(gcbo,'str')) 0 0 0];  c = min(max(c(1:3),0),1);
      else  c = get(h{3},'facecolor');
      end;
      d = round(100*c);
      k = 0;
      for j=1:3
        z = plt('slider',s(j),'get');   obj = plt('slider',s(j),'get','obj');
        if z~=d(j) & ~k bk = [1 1 0];  k = j;  if ~y d(k)=z; end;
        else            bk = [0 1 1];
        end;
        set(obj(5),'backg',bk);
      end;
      if ~k set(obj(5),'backg',[1 1 0]); return; end;
      c = d/100;
    end;
  end;
  set(h{3},'facecolor',c);
  d = round(c*100);
  for j=1:3  plt('slider',s(j),'set','val',d(j)); end;
  ph = h{4};  phc = [0 0 0];  phc(k) = d(k)/10;  v = mod(k,3)+1;  w = mod(v,3)+1;
  for row = 0:10
    phc(v) = row;
    for col = 0:10
      phc(w) = col;  set(ph(row+1,col+1),'FaceColor',phc/10);
    end;
  end;
  t = ctriple(c);
  for k=2:2:length(hq)
    p = hq{k-1};  h = hq{k};
    if ~iscell(p) p = {p}; end;
    for j=1:length(p)
      q = p{j}; if strcmp(q(1:3),'str') set(h,q,t); else set(h,q,c); end;
    end;
  end;
  ccf = getappdata(f,'ccf');
  if length(ccf)>1 eval(ccf); end;

case 434
  if Narg<2 [fi pth] = uigetfile('plt.plt','Select plt figure to open');
            if isnumeric(fi) return; end;
            fi = [pth fi];
  else      fi = varargin{2};
  end;
  ydat = [];
  feval('load',fi,'-mat');
  if isempty(ydat)
    p = find(fi=='.');
    if length(p) f = [fi(1:p(end)) 'mat']; end;
    dos(['copy ' fi ' ' f ' > NUL:']);
    load(f); delete(f);
  end;
  xdat = [xdat'; ydat'];
  plt(xdat{:},params{:});

case 431
  if Narg<2 [fi pth] = uiputfile('plt.plt','Save plt figure as');
            if isnumeric(fi) return; end;
            fi = [pth fi];
  else      fi = varargin{2};
  end;
  AX    = findobj(gcf,'tag','click');
  CurID = get(AX,'user');
  AX2   = findobj(gcf,'YAxisLoc','right','user',CurID);
    lh = getappdata(gcf,'Lhandles');
    xdat = get(lh,'x');  ydat = get(lh,'y');
    xlm =  get(AX,'xlim');   ylm =  get(AX,'ylim');
    xymult = getappdata(gcf,'xymult');
    if xymult(1) ~= 1
       mult = 1/xymult(1);
       for k = 1:length(xdat) xdat{k} = mult * xdat{k}; end;
       xlm = xlm * mult;
    end;
    ym = 1;
    for k = 1:length(ydat)
      mult = xymult(k+1);
      if mult ~= 1  
         mult = 1/mult;
         if ym ylm=ylm*mult; ym=0;  end;
         ydat{k} = mult * ydat{k};
      end;
    end;
    v = get(lh,'vis');  v = [v{:}];  
    params = [getappdata(gcf,'params') { ...
             'pos'        get(gcf,'pos') ...
             'DIStrace' v(find(v=='o')+1)=='f' ...
             'Xlim'     xlm ...
             'Ylim'     ylm     }];
    if length(AX2) params = [params {'YlimR' get(AX2,'ylim')}]; end;
    ver = '13Mar15';
    save(fi,'xdat','ydat','params','ver');

case 518
  y2    = varargin{2};
  AX    = findobj(gcf,'tag','click');
  CurID = get(AX,'user');
  AX2   = findobj(gcf,'YAxisLoc','right','user',CurID);
  AXrl  = [AX AX2];
  if isempty(AX2) AX2 = -1; end;
  CurMain = getappdata(0,'CurMain');
  Hc = get(CurMain(CurID),'user');
  if ischar(y2)
    switch sum(y2)
    case 430
      h = get(AX2,'ylab');  s = get(h,'str');
      if s(1)==92 & s(2)=='d'  s = s(6:end-5);
      else                             s = ['\div ' s ' \div'];
      end;
      set(h,'str',s);

    case 733
      xl = get(AX,'xlim'); yl = get(AX,'ylim');
      d = [-.2 .2]; e = [.1 10]; dr = 'o';
      clk = sum(get(gcf,'SelectionT'));
      if clk==321 | (clk==434 & getappdata(AX,'dir')=='i')
         d = d/-1.4; e = 1./e; dr = 'i';
      end;
      xs = get(AX,'Xscale');  ys = get(AX,'Yscale');
      if xs(2)=='i' xl = xl + diff(xl)*d;
      else          xl = e.*xl;  if diff(xl)<=0 xl = xl./e; end;
      end;
      if ys(2)=='i' yl = yl + diff(yl)*d;
      else          yl = e.*yl;  if diff(yl)<=0 yl = yl./e; end;
      end;
      set(AXrl,'xlim',xl);
      set(AX  ,'ylim',yl);
      if ishandle(AX2)
        axl = get(get(AX2,'ylab'),'str');
        if axl(1)~=92 | axl(2) ~= 'd';
          yl = get(AX2,'ylim');
          if ys(2)=='i' yl = yl + diff(yl)*d;
          else          yl = e.*yl;  if diff(yl)<=0 yl = yl./e; end;
          end;
          set(AX2,'ylim',yl);
        end;
      end;
      plt('grid',AX);
      setappdata(AX,'dir',dr);
      axes(AX); evalQ(get(Hc(5),'user'));
      xView = getappdata(gcf,'xView');
      if length(xView) set(xView{1},'x',get(AX,'xlim')); end;

    case 427
      if Narg>2 in3 = s2i(varargin{3}); else in3 = -1; end;
      switch in3
        case 3,
          cFIGbk = get(gcf,'color');
          cPLTbk = get(AX,'color');
          if isstr(cPLTbk) cPLTbk = get(AX2,'color'); end;
          cXYax  = get(AX,'xcol');
          cXYlbl = get(get(AX,'xlab'),'color');
          cDELTA = findobj(gcf,'tag','DeltaC'); cDELTA = get(cDELTA(1),'color');
          cTRACE = get(getappdata(gcf,'Lhandles'),'color');
          cTRACE = reshape([cTRACE{:}],3,length(cTRACE))';
          cFile = findobj(gcf,'style','push','str','D');  cFile  = get(cFile(end),'tag');
          gr = findobj(gcf,'user','grid');  gr = gr(1);   cGRID  = getappdata(gr,'clr');
          GridSty = get(gr,'linest');                  GridEr = get(gr,ERAS);
          if isempty(cFile) [cFile pth] = uiputfile('*.mat','Select file for saving colors');
                            cFile = [pth cFile];
          end;
          if sum(cFile)
             save(cFile,'cFIGbk','cPLTbk','cXYax','cXYlbl','cDELTA','cTRACE','cGRID','GridSty','GridEr');
             msgbox(['This program will now use colors saved in file ' cFile],'modal');
          else disp('No file was selected'); end;
          return;
        case 4,
          close(findobj('name','Color Pick'));
          g = s2i(get(gcf,'tag'));  g = findobj('type','fig','Number',g);
          gr = findobj(g,'user','grid');
          if length(gr)
            gr = gr(1);
            er = getappdata(gr,'er'); if isempty(er) closereq; return; end;
            setappdata(gr,'er',[]);
            if er=='n' er=ERANOR; else er=ERAXOR; end;
            c = get(gr,'color');
            setappdata(gr,'clr',c);
            ax = getappdata(g,'axis');  axr=ax(end);  r = get(axr,'YAxisLoc');
            if r(1)=='r'
              vis = get(axr,'vis');
             if vis(2)=='n' gx = get(g,'color'); else gx = get(axr,'color'); end;
              c = bitxor(round(255*c),round(255*gx))/255;
            end;
            set(gr,'color',c,ERAS,er);
          end;
          closereq;
          return;
      end;
      g = gcf;
      mb = findobj(g,'tag','MenuBox');
      callc = 'plt MarkEdit 3;';
      if in3==2
        figure('menu','none','numberT','off','back','off','resize','off','pos',auxLoc(302,60),'color',[0,.4,.4],'name','Edit figure colors','tag',get(gcf,'tag'),...
                      'closereq','plt click mark 4;');
        ps = 'Figure background|Plot background|Axis color|Axis color (right)|Axis labels|Delta cursor|Grid color|Grid style';
        pu = uicontrol('sty','pop','str',ps,'pos',[80 35 140 20],'call','plt MarkEdit 2;');
        ed = uicontrol('sty','edit','pos', [80  5 140 22],'fontw','bol','call',callc,'buttond',callc);
        gr = findobj(g,'user','grid');  gr = gr(1);
        er = get(gr,ERAS);  setappdata(gr,'er',er(1));
        set(gr,ERAS,ERANOR);
        set(gr,'color',getappdata(gr,'clr'));
        set([ed pu],'backg',[.8,.8,.9],'foreg','black');
        setappdata(ed,'m',{'str',ed,'color',[g mb],'Figure background'});
        setappdata(ed,'pop',pu); setappdata(pu,'edt',ed); setappdata(pu,'obj',[]);
        set(ed,'str',ctriple(get(g,'color')));
        ax = getappdata(g,'axis');  a = ax(1);  ar = ax(end);
        if isstr(get(a,'color')) ap=ar; else ap=a; end;
        lb = [get(a,'xlab') get(a,'ylab') get(mb,'child')'];
        dc = findobj(g,'tag','DeltaC');
        set(pu,'user',{[g mb];
                     ap;
                     a;
                     ar;
                     lb;
                     dc(1)
                     gr
                     gr
                    });
        return;
      end;
      sb = 0;
      if Narg>3
        cid = varargin{4};
        [lm lh] = plt('cursor',cid,'get','activeLine');
        lhs = getappdata(gcf,'Lhandles');
        ln = find(lhs == lh);
        nl = length(getappdata(AX,'Lhandles'));
        if ln > nl
          Hc = get(CurMain(cid),'user');
          sb = ln-1;
        end;
      end;
      hix = Hc(4);
      misc = get(hix,'user');
      actv = misc(4);
      iact = 15 + actv;
      hact = Hc(iact);
      if sum(get(gcf,'SelectionT'))==321 | Narg>2
        tx = [];  hb = [];  tid = [];
        if sb h=[]; else h = findobj(gcf,'user','TraceID'); end;
        if length(h) h = flipud(get(h,'child'));
                     tx = findobj(h,'type','text');  tid = tx(actv);
                     h = findobj(h,'type','line');
                     bt = get(h,'button');
                     if iscell(bt)
                       hb = h(find(cellfun('length',bt)));
                       if length(hb)>=actv tid = [tid hb(actv)]; end;
                     end;
        end;
        c = hact;
        h = get(c,'user'); h = h{1};
        if in3>=0 alll = in3;
        else
              alll = get(Hc(13),'vis');
              if alll(2)=='n'
                 h = findobj(gcf,'str','D');  h = get(h(end),'callb');
                 evalQ(h); plt click mark 1;
                 return;
              end;
              alll = 0;
        end;
        if alll
          c = [];  h = [];  a = [];
          for i=15+1:length(Hc)
            lh = get(Hc(i),'user');
            c = [c Hc(i)];  h = [h lh{1}];
          end;
          tid = [tx; hb]';
          fname = 'Edit all lines';
        else
          if length(tid) tx=tid(1); a=[];
          else a=get(h,'par'); tx=get(a,'ylab'); tid=tx;
          end;
          fname = sprintf('Edit Line %d (%s)',actv+sb,deblank(get(tx,'str')));
        end;
        callp = 'plt MarkEdit 0;';
        calle = 'plt MarkEdit 1;';
        figure('menu','none','numberT','off','back','off','resize','off','pos',auxLoc(302,85),'color',[0,.4,.4],'name',fname,'tag',get(gcf,'tag'),...
               'closereq','plt click mark 4;');
        props = 'Color|LineStyle|LineWidth|Marker|MarkerSize|Xdata|Ydata|Zdata';
        ed1 = uicontrol('sty','edit','pos',[8 5 130 22],'fontw','bol','call',calle,'buttond',calle);
        pu1 = uicontrol('sty','pop','str',props,'pos',[8 35 130 20],'call',callp);
        uicontrol('sty','text','str', 'Line properties:','pos',[8 64 130 17]);
        ed2 = uicontrol('sty','edit','pos',[145 5 150 22],'fontw','bol','call',calle,'buttond',calle);
        pu2 = uicontrol('sty','pop','str',props,'pos',[145 35 150 20],'call',callp,'Val',5);
        uicontrol('sty','text','str', 'Cursor properties:','pos',[145 64 150 17]);
        set([ed1 ed2 pu1 pu2],'backg',[.8,.8,.9],'foreg','black');
        setappdata(ed1,'pop',pu1); setappdata(pu1,'edt',ed1); setappdata(pu1,'obj',h);
        setappdata(ed2,'pop',pu2); setappdata(pu2,'edt',ed2); setappdata(pu2,'obj',c);
        setappdata(ed1,'m',{'str',ed1,'color',[h tid],'ycolor',a,'Line color'});
        setappdata(ed2,'m',{'str',ed2,'color',c,'Cursor color'});
        set(ed1,'str',ctriple(get(h(1),'color')));
        set(ed2,'str',get(c(1),'MarkerSize'));
        axes(AX);
      else
        hiy = Hc(6);
        axes(AX);
        x = get(hact,'x'); y = get(hact,'y');
        p = get(hact,'par');
        rlim = [];
        if p ~= AX
          rlim = get(p,'ylim');
          ylim = get(AX,'ylim');
          y = ylim(1) + diff(ylim) * (y - rlim(1)) / diff(rlim);
        end;
        l = line(x,y,'marker','s');
        t = text(x,y,['   (' get(hix,'str') ', ' get(hiy,'str') ')'],'units','data',...
                 'fontsi',get(p,'fontsi'),'user',l,'buttond','plt misc marker;');
        set([t l],'color',get(hact,'color'),'tag','mark');
        if length(rlim)
          set(l,'tag','markR','user',{t AX p ylim rlim});
        end;
      end;
    case 674
      if sum(get(gcf,'SelectionT'))==649
        hl = findobj(gcf,'buttond','plt click TGLlogy;');
        if strcmp(get(AX,'Yscale'),'log')
             sc='linear'; st='LinY';
        else sc='log';    st='LogY';
             y = get(AX,'ylim');    if y(1)<=0 set(AX,'ylim',abs(y(2))*[.001 1]); end;
             if ishandle(AX2)
               y = get(AX2,'ylim'); if y(1)<=0 set(AX2,'ylim',abs(y(2))*[.001 1]); end;
             end;
        end;
        set(AXrl,'Yscale',sc);  set(hl,'str',st);
      else plt('hcpy','init',gcf);
      end;
    case 673
      ax = getappdata(gcf,'axis');
      if sum(get(gcf,'SelectionT'))==649
        hl = findobj(gcf,'buttond','plt click TGLlogx;');
        p = get(AX,'pos');  p = p(1);
        axe = [];  nax = length(ax);
        if ishandle(AX2) nax=nax-1; end;
        for k = 2:nax
          q = get(ax(k),'pos');
          if p == q(1) axe = [axe ax(k)]; end;
        end;
        axf = [AXrl axe];
        if strcmp(get(AX,'Xscale'),'log')
             sc='linear'; st='LinX';
        else sc='log';    st='LogX';
             x = get(AX,'xlim'); if x(1)<=0 set(axf,'xlim',abs(x(2))*[.001 1]); end;
        end;
        set(axf,'Xscale',sc); set(hl,'str',st);
        for k = 1:length(axe) plt('grid',axe(k)); end;
      else
        h = getappdata(gcf,'Lhandles');
        x = get(h,'x'); y = get(h,'y');
        if length(h)>1 set(h,{'x'},y,{'y'},x);
        else           set(h,'x',y,'y',x);
        end;
        for k = 1:length(ax)
          a = ax(k);
          set(a,'xlim',get(a,'ylim'),'ylim',get(a,'xlim'));
          x = get(a,'xlab');  y = get(a,'ylab');
          sx = get(x,'str'); set(x,'str',get(y,'str')); set(y,'str',sx);
        end;
      end;
    case 653
      if sum(get(gcf,'SelectionT'))==649
        set(AX,'TickLen',(1-plt('grid',AX,'toggle'))*[.01 .025]);
        a = getappdata(gcf,'axis');
        if ishandle(AX2) a(end) = []; end;
        for k=2:length(a)
          set(a(k),'TickLen',(1-plt('grid',a(k),'toggle'))*[.01 .025]);
        end;
      else
        ax = getappdata(gcf,'axis');  n = length(ax);
        r = get(ax(n),'YAxisLoc');
        if r(1)=='r' axr=ax(n); n=n-1; else axr=-1; end;
        for k = 1:n
          g = findobj(ax(k),'user','grid');  if isempty(g) continue; end;
          c = getappdata(g,'clr');  st = get(g,'linest'); 
          inv = c(1) > .5;
          if st(1)=='-'  c = 2*c;  if inv c=c-1; end;  c=max(min(c,1),0); cc=c; er=ERANOR; sty=':';
          else           c = c/2;  if inv c=c+.5; end; c=max(min(c,1),0); cc=c; er=ERAXOR; sty='-';
          end;
          if k==1 & ishandle(axr) & Mver<8.4
             vis = get(axr,'vis');
             if vis(2)=='n' gx = get(gcf,'color'); else gx = get(axr,'color'); end;
             c = bitxor(round(255*c),round(255*gx))/255;
          end;
          setappdata(g,'clr',cc);
          set(g,'linest',sty,'color',c,ERAS,er);
        end;
      end;
    case 668
      if sum(get(gcf,'SelectionT'))==649
        f = get(gcf,'menu');
        if   f(1)=='f' delete(findobj(gcf,'type','uimenu')); set(gcf,'menu','none');
        else set(gcf,'menu','fig');
             v = get(0,'ShowHidden');
             set(0,'ShowHidden','on');
             a = findobj(gcf,'label','&File');
             uimenu(a,'Label','p&lt  save','separator','on','call','plt save;');
             uimenu(a,'Label','pl&t  open','call','plt open;');
             c = get(a,'child');
             set(a,'child',c([4:end-3 1:3 end-2:end]));
             a = findobj(gcf,'label','&Help');
             uimenu(a,'Label','p&lt help','separator','on','call','plt help;');
             c = get(a,'child'); set(a,'child',c([2:end 1]));
             set(c(end),'separator','on');
             set(0,'ShowHidden',v);
             m = uimenu('label','&plt');
             cb = {'plt click mark 0;' ...
                   'plt click mark 1;' ...
                   'plt click mark 2;' ...
                   'plt click mark 3;' ...
                   'plt misc clkright TGLmenu;' ...
                   'plt misc clkright TGLlogx;' ...
                   'plt("hcpy","init",gcf);' ...
                   'h=findobj(gcf,"str","O"); eval(get(get(h(end),"ui"),"callb"));' ...
                   'delete([findobj(gcf,"tag","mark"); findobj(gcf,"tag","markR")]);' ...
                   'a=getappdata(gcf,"axis"); setappdata(a(1),"DualCur",plt("cursor",get(a(1),"user"),"get","activeLine",0));'... 
                   'plt move;' ...
                   'plt move res;' };
             lb = {'<1>&Edit line (<2>e<3>) <4>Rclick Mark<5>' ...
                   '<1>Edit &all lines (<2>a<3>)' ...
                   '<1>Edit &figure colors (<2>f<3>) <4>Delta+Rclick Mark<5>' ...
                   '<1>&Save figure colors (<2>s<3>)' ...
                   '<1>&Cursor Data Window (<2>c<3>)  <4>Rclick Menu<5>' ...
                   '<1>S&wap X/Y axes (<2>w<3>) <4>Rclick LinX<5>' ...
                   '<1>&Hardcopy (<2>h<3>) <4>Rclick LinY<5>' ...
                   '<1>&Toggle line smoothing (<2>t<3>)  <4>Rclick "o"<5>' ...
                   '<1>&Delete cursor annotations (<2>d<3>)  <4>Delta+Rclick "o"<5>' ...
                   '<1>Set d&ual cursor (<2>u<3>)' ...
                   '<1>Toggle &Reposition mode (<2>r<3>)  <4>Rclick Delta<5>'...
                   '<1>Reposition &Grid size (<2>g<3>) <4>Delta+Rclick Delta<5>' };
             sp = {'of' 'of' 'of' 'of' 'on' 'of' 'of' 'of' 'of' 'of' 'on' 'of'};
             rp = {'<1>' '<html>' '<2>' '<u>' '<3>' '</u>' '<4>' '<font color="blue"><i>' '<5>' ''};
             v = version;
             if v(1)=='6' rp = {'<1>' '' '<2>' '' '<3>' '' '<4>' '[' '<5>' ']'}; end;
             for k=1:length(cb)
               s = lb{k};
               for j=1:2:length(rp) s = strrep(s,rp{j},rp{j+1}); end;
               uimenu(m,'label',s,'call',strrep(cb{k},'"',''''),'sep',sp{k});
             end;
        end;
      else
        cf = gcf;
        l = getappdata(cf,'Lhandles');
        if isempty(l) disp('No data to display'); return; end;
        tn = get(flipud(findobj(findobj(gcf,'user','TraceID'),'type','text')),'str');
        a = getappdata(gcf,'axis');
        if length(a)>1
           a = get(a,'ylab');
           for k = 2:length(a) a{k} = get(a{k},'str'); end;
           tn = [tn; a(2:end)];
        end;
        t = ' index';
        x = [];
        d = [];
        n = 0;
        nc = 0;
        np = 0;
        for k = 1:length(l)
          if strcmp(get(l(k),'vis'),'off') continue; end;
          xx = get(l(k),'x');  nn = length(xx);  ix = [1 min(nn,2) max(1,nn-1) nn]; xxi = xx(ix);
          if nn~=np | ~isequalNaN(x(ix),xxi(:))
            np = nn;  x = xx;
            sz = size(x);   if sz(1)==1 x = transpose(x); end;
            ex = [];
            if     nn>n  d = [d; NaN+ones(nn-n,nc)];  n = nn;
            elseif nn<n  ex = NaN+ones(n-nn,1);
            end;
            d = [d [x; ex]];  t = [t '     X    '];  nc=nc+1;
          end;
          y = get(l(k),'y');  sz = size(y); if sz(1)==1 y = transpose(y); end;
          d = [d [y; ex]];  nc=nc+1;
          if length(tn)>=k
            s = tn{k};
            if length(s)>=9 s = [' ' s(1:9)];
            else m = 9-length(s);  mf = floor(m/2);
                 s = [blanks(1+m-mf) s blanks(mf)];
            end;
            t = [t s];
          end;
        end;
        s = {}; clr = get(cf,'color');
        w = 70*nc + 80;
        sz = get(0,'screens');  p = get(cf,'pos');
        if p(3)<1 p = p .* sz([1 2 1 2]); end;
        x1 = p(1) + p(3) + 12;  wa = sz(3) - x1;
        if p(1) > wa  x1 = 5;  wa = p(1)-5; end;
        y1 = p(2);  ha = p(4);
        if wa > w
          if x1==5 x1=wa-w-7; end; wa=w;
        else
          Y1 = p(2) + p(4) + 35;  Ha=sz(4)-Y1-80;
          if p(2) > Ha+80  Y1=59; Ha=p(2)-94; end;
          if (Ha>400) | ((Ha>200) & (sz(4)>w))
            ha = Ha; y1 = Y1; x2 = sz(4)-20; wa = min(x2,w);
            x1 = x2-wa;  x1 = x1 - min(x1,x2-p(1)-p(3));
          elseif wa<210 wa = min(w,sz(4)-10); x1=5;
          end;
        end;
        pos = [x1 y1 wa ha];
        fnt = listfonts; fnt = [fnt{:}];
        if     length(findstr(fnt,'Consolas'))       fnt = 'Lucida Sans Typewriter';
        elseif length(findstr(fnt,'Lucida Console')) fnt = 'Lucida Console';
        else                                         fnt = 'Courier';
        end;
        for k = 1:n  s = [s {[sprintf('%5d  ',k) strrep(prin('%8V  ',d(k,:)),'NaN','   ')]}]; end;
        f = figure('menu','none','numberT','off','back','off','name',get(cf,'name'),'pos',pos,'color',clr,'tag',get(gcf,'tag'));
        bx = uicontrol('style','listbox','pos',[1 1 wa ha-30],...
        'foreg',[1 1 0],'backg',[0 0 0],'user',{s t},'FontName',fnt);
        fs = uicontrol('sty','pop','str',prin('{fontsize: %d!row}',4:18),'Val',6,'pos',[wa-100 ha-19 85 16],...
          'call','set(findobj(gcf,''sty'',''l''),''fontsize'',3+get(gcbo,''value''));');
        sv = uicontrol('str','save','user',bx,'pos',[15 ha-25 50 20],'call',...;
           '[f p]=uiputfile(''pltData.txt''); s=get(get(gcbo,''user''),''user''); t=s{2}; s=s{1}; prin(-[p f],''%s\n'',t,s{:});');
        set([bx fs sv],'unit','nor');
        setappdata(cf,'bx',bx);
        figure(cf);  plt('cursor',CurID,'MVcur');
      end;

    case 242
      misc = get(Hc(4),'user');
      lH = get(Hc(15+misc(4)),'user');
      hix = Hc(4);  hiy = Hc(6);  hix2 = Hc(5);  hiy2 = Hc(7);
      x = get(lH{1},'x');  y = get(lH{1},'y');
      xlim = get(Hc(14),'xlim');
      xx = find(x >= xlim(1)  &  x <= xlim(2));
      y = y(xx); y = y(~isnan(y)); ly = length(y);
      p = findobj(gcf,'user','idcur');  ps = get(p,'str');
      switch sum(ps)
        case 286,  s = 'RMS'; r = Pftoa('%7w',sqrt(sum(y.^2)/ly));
        case 242,  s = 'y/x';
                    xr = s2d(get(hix,'str'));
                    r = getappdata(p,'idcur');
                    r = Pftoa('%7w',s2d(r{2})/xr);
                    if sum(get(hiy2,'vis'))==221 & sum(get(hix2,'vis'))==221
                       xr = str2num(get(hix2,'str'));
                       set(hiy2,'str',s2d(get(hiy2,'str'))/xr);
                    end;
        case 288,  s = '\surdx^2+y^2';
                    xr = s2d(get(hix,'str'));
                    q = getappdata(p,'idcur');
                    r = Pftoa('%7w',abs(xr + s2d(q{2})*1j));
                    if sum(get(hiy2,'vis'))==221 & sum(get(hix2,'vis'))==221
                       xr = str2num(get(hix2,'str'));
                       set(hiy2,'str',abs(xr + s2d(q{3})*1j));
                    end;
        case 1110, r = getappdata(p,'idcur');  s = r{1}; set(hiy2,'str',r{3});  r = r{2};
        otherwise,  setappdata(p,'idcur',{ps get(hiy,'str') get(hiy2,'str')});
                    s = 'Avg';
                    r = Pftoa('%7w',sum(y)/ly);
      end;
      set(p,'str',s); set(hiy,'str',r);

    case 294
      in3 = varargin{3} - '0';
      a = getappdata(gcf,'Dedit');
      if isempty(a) return; end;
      switch in3
      case 1,
        if ishandle(a{3}) & length(get(a{3},'buttond')) return; end;
        Hc = get(CurMain(a{1}),'user');  misc = get(Hc(4),'user');
        hact = Hc(15 + misc(4));
        m = get(hact,'marker'); a{4} = m(1)+0;
        ix = a{2};
        if a{3}~=hact  a{3} = hact;  a{7} = 0;
                       if ix<4 ix=ix+6; a{2}=ix; end;
        end;
        setappdata(gcf,'Dedit',a);
        b = 'plt click EDIT 3;';
        if ix<4 fc = get(hact,'color'); else fc = 'auto'; end;
        mk = 'd>^d>^d>^';  lw = '111222111';  msz = getappdata(gcf,'EditCur');
        set(hact,'marker',mk(ix),'markerface',fc,'markersize',msz,'linewidth',lw(ix)-'0','buttond',b);
      case 2,
        hd = [findobj(gcf,'style','edit','vis','on'); findobj(gcf,'user','grid')];
        set(hd,'vis','of');
        gc = findobj(hd,'buttond','plt click EDIT 2;');
        if length(gc)>1 gc = gco; end;
        fd = find(hd==gc);
        if length(fd) hd(fd) = [];  hd = [gc; hd]; end;
        hpop = getappdata(gcf,'epopup');
        setappdata(hpop,'EdHide',hd);
        pe = get(gc,'pos');  pa = get(hpop,'par');
        set(pa,'pos',[pe(1) .01 .13 .45]);
        set(gcf,'SelectionT','normal');
        evalQ(get(hpop,'buttond'));
        set(findobj(pa,'str','Cancel'),'color',[0 0 0],'fontw','bol');
      case 3, set(gcf,'WindowButtonMotionFcn','plt click EDIT 4;',...
                      'WindowButtonUpFcn',['set(gcf,''WindowButtonMotionFcn'','''',''WindowButtonUpFcn'','''');' ' plt click EDIT 5;']);
      case 4,
        i = mod(a{2},3);
        Hc = get(CurMain(a{1}),'user');
        fmt = get(Hc(3),'user');
        y = get(get(a{3},'par'),'currentp'); x=y(1,1); y=y(1,2);
        if i~=2 set(a{3},'y',y); set(Hc(6),'str',Pftoa(fmt(2,:),y)); end;
        if i~=0 set(a{3},'x',x); set(findobj(gcf,'sty','edit','str',get(Hc(4),'str')),...
                                                'str',Pftoa(fmt(1,:),x));
        end;
        m = getappdata(gcf,'MotionEdit'); if length(m) feval(m,0); end;
      case 5,
        cid = a{1};  ix = a{2};  hact = a{3};
        Hc = get(CurMain(cid),'user');  misc = get(Hc(4),'user');
        lh = get(hact,'user');  lk = lh{1};
        x = get(lk,'x');  y = get(lk,'y');   i = misc(1);
        x1 = get(hact,'x');  y1 = get(hact,'y');
        if     ix<4 ilast = a{7};
                    x0 = x(ilast);  y0 = y(ilast);
                    inc = 1 - 2*(ilast>i);  steps = max(1,abs(ilast-i));
                    dx = (x1-x0)/steps; dy = (y1-y0)/steps;
                    for k = ilast:inc:i x(k)=x0; y(k)=y0; x0=x0+dx; y0=y0+dy; end;
        elseif ix<7 ylim = get(Hc(14),'ylim');
                    if y1<ylim(1) x(i) = [];  y(i) = [];
                    else          x = [x(1:i) x1 x(i+1:end)];
                                  y = [y(1:i) y1 y(i+1:end)];  i=i+1;
                    end;
        else        x(i)=x1; y(i)=y1;
        end;
        set(lk,'x',x,'y',y); setappdata(gcf,'NewData',i);
        if (Narg==4 & mod(a{2},3)==1) | sum(get(gcf,'SelectionT'))==321
          return;
        end;
        a{7} = i;  a{8} = hact;  setappdata(gcf,'Dedit',a);
        set(a{3},'marker',char(a{4}),'markerface','auto','markersize',a{5},'linewidth',a{6},'buttond','');
        plt('cursor',cid,'update',i);
      end;

    case 511
      hpop = getappdata(gcf,'epopup');
      hd = getappdata(hpop,'EdHide');
      set(hpop,'vis','of');  set(hd,'vis','on');
      if nargin==2 ix = plt('pop',hpop,'get','index');
                   cb = get(hd(1),'call');
                   if ~iscell(cb) return; end;
                   cid = cb{3};
      else         ix = varargin{3};
                   if ischar(ix) ix = ix-'0'; end;
                   cid = getappdata(gcf,'cid');
      end;
      a = getappdata(gcf,'Dedit'); a3 = a{3};
      switch ix
      case 1,
         if sum(get(gcf,'SelectionT'))==649 plt('click','mark','-1',cid);
         else         plt click mark 2;
         end;
      case 2,
         multi = getappdata(gcf,'multi');
         if length(multi) delete(multi); setappdata(gcf,'multi',[]); return; end;
         Lha = getappdata(gcf,'Lhandles'); n = length(Lha);  cL = get(Lha,'color');
         c = (get(Lha(1),'color') + get(Lha(min(n,4)),'color'))/2;
         nn = zeros(1,n);  nnn = [nn; nn];  par = get(Lha,'par');
         tx = text(nn,nn,'');
         set(tx,{'par'},par,'unit','data',{'color'},cL,'tag','  %7w'); 
         mr = line(nnn,nnn,'color',c,'marker','o');  set(mr,{'par'},par);
         vl = line([0 0],1e9*[-1 1],'color',c,'linest',':','par',AX);
         mc = [tx; mr; vl];  set(mc,'buttond',get(AX,'buttond'));
         p = getappdata(gcf,'mcProps');
         for k = 1:2:length(p)
           prop = p{k};  prp = prop(2:end); val = p{k+1};
           switch prop(1)
             case '|',  set(vl,prp,val);
             case '+',  if iscell(val) prp  = {prp};  end; set(mr,prp,val);
             otherwise, if iscell(val) prop = {prop}; end; set(tx,prop,val);
           end;
         end;
         setappdata(gcf,'multi',mc);  plt('cursor',CurID,'update');
      case 3,
         xView = getappdata(gcf,'xView');
         if nargin<4
           if length(xView) l = xView{1};  delete([l get(l,'par')]);
                            set(AXrl,'pos',xView{2});
                            setappdata(gcf,'xView',[]);
                            return;
           end;
           numLines=length(Hc)-15;
           minx=+inf;  maxx=-inf;
           for i=1:numLines hLine = get(Hc(15+i),'user');
                            x = get(hLine{1},'x');
                            minx  = min(minx,min(x));
                            maxx  = max(maxx,max(x));
           end;
           p = get(AX,'pos');  g = get(gcf,'pos');  g = g(4);
           c1 = get(gcf,'color');  c2 = get(AXrl(end),'color');  d = c2 > 0.5;
           p2 = p - [0 0 0 20/g];  set(AXrl,'pos',p2);
           p3 = p2 + [0 p2(4)+10/g 0 12/g-p2(4)];
           v = axes('pos',p3,'xlim',[minx maxx],'ylim',[1 3],'color',c2,'xcol',c2,'ycol',c1,'XtickLabel',' ','YtickLabel',' ','TickLen',[0 0]','buttond','plt click Yedit 3 0;');
           l = line(get(AX,'xlim'),[2 2],'color',c2 + .33*(~d - d),'linewidth',16,'buttond','plt click Yedit 3 1;');
           if Mver >= 8.4 set(v,'clippingStyle','rectangle'); end;
           r = getappdata(gcf,'xvProps');
           for k = 1:2:length(r)
             prop = r{k};  val = r{k+1};
             if prop(1) == '+'  prop(1)='';  if isempty(prop) prop='pos'; val = val+p3; end;
                                set(v,prop,val);
             else               set(l,prop,val);
             end;
           end;
           setappdata(gcf,'xView',{l p 0});
           axes(AX);
         else
           l = xView{1};  v = get(l,'par');  x = get(v,'currentp');  x = x(1,1);  lx = get(l,'x');
           v4 = varargin{4}-'0';  if v4==1 & sum(get(gcf,'SelectionT'))==321 v4=0; end;
           switch v4
             case 0,
                     lx((x>mean(lx))+1) = x;  x = lx;
                     set(gcf,'WindowButtonMotionFcn','plt click Yedit 3 2;','WindowButtonUpFcn','set(gcf,''WindowButtonMotionFcn'','''',''WindowButtonUpFcn'','''');');
             case 1,
                     if sum(get(gcf,'SelectionT'))==434 x = get(v,'xlim');
                     else xView{3} = x-lx(1);  setappdata(gcf,'xView',xView); x=lx;
                          set(gcf,'WindowButtonMotionFcn','plt click Yedit 3 3;','WindowButtonUpFcn','set(gcf,''WindowButtonMotionFcn'','''',''WindowButtonUpFcn'','''');');
                     end;
             case 2, lx((x>mean(lx))+1) = x;  x = sort(lx);
             case 3, x = x-xView{3};  x = [x x+diff(lx)];
           end;
           set(l,'x',x);  set(AXrl,'xlim',x);
           axes(AX); evalQ(get(Hc(5),'user'));
         end;
      case 4,
        if ishandle(a3)
          set(a3,'marker',char(a{4}),'markerface','auto','markersize',a{5},'linewidth',a{6},'buttond','');
        end;
      otherwise,
        b = 'plt click EDIT 3;';
        if ishandle(a3) & strcmp(b,get(a3,'buttond')) plt click EDIT 5; end;
        ix = ix - 4;
        Hc = get(CurMain(cid),'user');  misc = get(Hc(4),'user');
        hact = Hc(15 + misc(4));
        if ix<4 & hact~=a{8}
          ix=ix+6;
        end;
        m = get(hact,'marker'); m = m(1)+0;
        setappdata(gcf,'Dedit',[{cid ix hact m get(hact,'markersize') get(hact,'linewidth')} a(7:8)]);
        if ix<4 fc = get(hact,'color'); else fc = 'auto'; end;
        mk = 'd>^d>^d>^';  lw = '111222111';  msz = getappdata(gcf,'EditCur');
        set(hact,'marker',mk(ix),'markerface',fc,'markersize',msz,'linewidth',lw(ix)-'0','buttond',b);
      end;
    end;
    plt('grid',AX);
  else y8 = y2;  j = y8(2);
       clkType = sum(get(gcf,'SelectionT'));
       p = {'color'; 'marker'; 'linest'; 'linewidth'}; q = {[0 .3 .3] 'none' '-' 9};
       if clkType==649 | Narg>2
          k = y8(1);  s = get(k,'vis'); s = s(2)=='n';
          if s set(k,'vis','of'); else set(k,'vis','on'); end;
          if ishandle(j) & length(get(j,'par'))
            mk = getappdata(j,'mk');
            if s   set(j,'fonta','ita','fontw','nor'); set(mk,p,q);
            else   set(j,'fonta','nor','fontw','bol'); set(mk,p,get(k,p));
            end;
          end;
       else
          for k = findobj(get(j,'par'),'type','text')'
            t = get(k,'buttond');
            if length(t{end})==2
              t = t{end}(1);
              mk = getappdata(k,'mk');
              if clkType==434 | k==j  set(t(1),'vis','on');  set(k,'fonta','nor','fontw','bol');  set(mk,p,get(t,p));
              else                      set(t(1),'vis','of');  set(k,'fonta','ita','fontw','nor');  set(mk,p,q);
              end;
            end;
          end;
       end;
       DualCur = getappdata(AX,'DualCur');
       Lh = getappdata(AX,'Lhandles');
       ls = findobj(Lh,'par',AX);
       v = 'off';
       for k=1:length(ls)
         if strcmp(get(ls(k),'vis'),'on') v = 'on';  break; end;
       end;
       set(get(AX,'ylab'),'vis',v);
       if ishandle(AX2)
         ls = findobj(Lh,'par',AX2);
         v = 'off';  leftC = get(AX2,'color');  gridXOR = leftC;
         for k=1:length(ls)
           if strcmp(get(ls(k),'vis'),'on')
             v = 'on';  leftC = 'none';  gridXOR = get(gcf,'color');
             set(AX2,'xlim',get(AX,'xlim'));
             break;
           end;
         end;
         set(AX2,'vis',v);  set(AX,'color',leftC);
         gridH = findobj(AX,'user','grid');
         if length(gridH)==1
           er = get(gridH,ERAS);
           if er(1)=='x'
             set(gridH,'color',bitxor(round(255*getappdata(gridH,'clr')),round(255*gridXOR))/255);
           end;
         end;
       end;
       TIDback = [];
       if Narg<3 & ishandle(j) TIDback = getappdata(get(j,'par'),'TIDcback'); end;
       if length(TIDback)
          setappdata(gcf,'OBJ',y2(1));  setappdata(gcf,'OBJ2',y2(2));
          evalRep(TIDback,{'@LINE' 'plt("misc",3)' '@TID' 'plt("misc",4)'});
       end;
       axes(AX);
  end;

case 708,
  ax = getappdata(gcf,'axis');  ax = ax(1);
  cid = getappdata(gcf,'cid');  nc = length(cid);
  ylbl = get(ax,'ylab');
  obj = getappdata(ylbl,'obj');
  if length(obj)
    setappdata(ylbl,'obj',[]);
    for k=1:nc plt('cursor',cid(k),'set','visON'); end;
    if obj(1) set(obj,'vis','on'); end;
  else
    for k=1:nc plt('cursor',cid(k),'set','visOFF'); end;
    a = findobj(gcf,'str','Zout');
    if length(a) a = get(a(1),'par');
                 obj = [a; get(a,'child')];
                 set(obj,'vis','of');
    else         obj = 0;
    end;
    setappdata(ylbl,'obj',obj);
  end;

case 449,
  h = getappdata(gcf,'Lhandles');  n = length(h);
  q = 3-cellfun('length',get(h,'vis'));
  if Narg<2 Ret1 = find(q'); return; end;
  v = zeros(n,1);  w=v;  e = varargin{2};
  if ischar(e) & length(e) v=v+1;  else v(e)=1; end;
  v = find(xor(v,q));
  t = findobj(gcf,'user','TraceID');
  if length(t)
    t = [flipud(findobj(t,'type','text')); w];
  end;
  for k=1:length(v) plt('click',[h(v(k)) t(v(k))],0); end;

case 439,
  if Narg==1
     if isempty(getappdata(gcf,'mv')) plt move init; end;
     if getappdata(gcf,'mv') plt move off; else plt move on; end;
     return;
  end;
  ui = getappdata(gcf,'uic');  Nui = length(ui);
  si = getappdata(gcf,'sli');  Nsi = length(si);
  ax = getappdata(gcf,'axi');  Nax = length(ax);
  tx = getappdata(gcf,'txt');  Ntx = length(tx);
  pp = getappdata(gcf,'pop');  Npp = length(pp);
  v = varargin{2};
  if ~ischar(v)
    for k=1:length(v)
      t = v(k);
      switch get(t,'type')
        case 'uicontrol', setappdata(gcf,'uic',[ui t]);
        case 'text',      setappdata(gcf,'txt',[tx t]);
        case 'axes',      setappdata(gcf,'axi',[ax t]);
      end;
    end;
    return;
  end;
  v = sum(v);
  switch v
  case 330
    f = gcf;
    if isempty(getappdata(f,'mv')) plt move init; end;
    r = getappdata(f,'snap');
    figure('menu','none','numberT','off','back','off','resize','off','pos',auxLoc(270,60),'color',[0,.4,.4],'name','SnapTo resolution','user',gcf,'tag',get(gcf,'tag'));
    plt('slider',[.035 .64],[r(1) 0 200 0 10000],'X resolution',...
      'g = get(gcf,''user''); r = getappdata(g,''snap''); setappdata(g,''snap'',[@VAL r(2)])',2);
    plt('slider',[.525 .64],[r(2) 0 200 0 10000],'Y resolution',...
      'g = get(gcf,''user''); r = getappdata(g,''snap''); setappdata(g,''snap'',[r(1) @VAL])',2);
    figure(f);
    return;
  case 436
    if isempty(getappdata(gcf,'snap')) setappdata(gcf,'snap',[100 100]); end;

    t = findobj(gcf,'type','text');
    for k=1:length(t);
      if isempty(getappdata(t(k),'ty')) tx = [tx t(k)]; setappdata(t(k),'ty','txt'); end;
    end;
    setappdata(gcf,'txt',tx)

    for k=1:Nsi setappdata(si(k),'ty','sli'); end;
    t = findobj(gcf,'type','uicontrol');
    for k=1:length(t)
      if isempty(getappdata(t(k),'ty')) ui = [ui t(k)]; setappdata(t(k),'ty','uic'); end;
    end;
    t = findobj(gcf,'type','uipanel');  r = findobj(gcf,'type','uitable');
    ui = [ui t' r'];
    for k=1:length(t) setappdata(t(k),'ty','uip'); end;
    for k=1:length(r) setappdata(r(k),'ty','uit'); end;
    setappdata(gcf,'uic',ui);

    t = findobj(gcf,'type','axes');
    for k=1:length(t)
      if isempty(getappdata(t(k),'ty')) ax = [ax t(k)]; setappdata(t(k),'ty','axi');
      end;
    end;
    setappdata(gcf,'axi',ax);
    setappdata(gcf,'mv',0);
    return;
  end;
  uip = {'buttond' 'ena'};
  if v==221
    if isempty(getappdata(gcf,'mv')) plt move init; end;
    setappdata(gcf,'mv',1);
    for k=1:Nui if sum(getappdata(ui(k),'ty'))==334
                  setappdata(ui(k),'CBsv',get(ui(k),'buttond'));
                  set(ui(k),'buttond','plt misc tidmv;');
                else
                  setappdata(ui(k),'CBsv',get(ui(k),uip));
                  set(ui(k),'ena','of','buttond','plt misc tidmv;');
                end;
    end;
    for k=1:Nsi setappdata(si(k),'CBsv',get(si(k),uip));
                set(si(k),'ena','of','buttond','plt misc slimv;');
                setappdata(si(k),'ty','sli');
    end;
    for k=1:Nax setappdata(ax(k),'CBsv',get(ax(k),'buttond'));
                set(ax(k),'buttond','plt misc tidmv;');
    end;
    for k=1:Ntx setappdata(tx(k),'CBsv',get(tx(k),'buttond'));
                set(tx(k),'buttond','plt misc txtmv;');
    end;
    for k=1:Npp setappdata(pp(k),'CBsv',get(pp(k),'buttond'));
                set(pp(k),'buttond','plt misc tidmv;');
    end;
    setappdata(gcf,'CBsv',get(gcf,'buttond'));  set(gcf,'buttond','');
    gr = findobj(gcf,'user','grid');
    if length(gr) gr=gr(end); setappdata(gr,'CBsv',get(gr,'buttond')); set(gr,'buttond','plt misc tidmv;'); end;
    set(findobj(gcf,'str','D','fontn','symbol'),'str',char(222));
  else
    set(findobj(gcf,'str',char(222),'fontn','symbol'),'str','D');
    setappdata(gcf,'mv',0);
    for k=1:Nui if sum(getappdata(ui(k),'ty'))==334
                     set(ui(k),'buttond',getappdata(ui(k),'CBsv'));
                else set(ui(k),uip,getappdata(ui(k),'CBsv'));
                end;
    end;
    for k=1:Nsi set(si(k),uip,     getappdata(si(k),'CBsv')); end;
    for k=1:Nax set(ax(k),'buttond',   getappdata(ax(k),'CBsv')); end;
    for k=1:Ntx set(tx(k),'buttond',   getappdata(tx(k),'CBsv')); end;
    for k=1:Npp set(pp(k),'buttond',   getappdata(pp(k),'CBsv')); end;
    set(gcf,'buttond',getappdata(gcf,'CBsv'));
    gr = findobj(gcf,'user','grid');
    if length(gr) gr=gr(end); set(gr,'buttond',getappdata(gr,'CBsv')); end;
  end;

case 534,   for f = findobj('type','fig')'
                   b = get(f,'buttond');
                   if iscell(b) & length(b)>3 & sum(b{4})==570
                     set(f,'closereq','closereq'); close(f);
                   end;
               end;
               setappdata(0,'CurMain',[]);
case 632,  set(flipud(findobj(findobj(gcf,'user','TraceID'),'type','text')),{'str'},varargin{2}(:));
case 774, Ret1 = '13Mar15'; Ret2 = 0;

otherwise,
cTRACE  = [0  1  0;   1  0  1;   0  1  1;   1  0  0;  .2 .6  1;
           1  1  1;   1 .6 .2;   0  0  1;   1 .2 .6;  .2  1 .6;
          .6  1 .2;  .6 .2  1;   1  1  0;   0 .6  0;  .6  0 .6;
           0 .6 .6;  .6 .6  0;  .7 .7 .7;  .6  0  0;  .2 .2 .7;
          .5 .5 .5;  .7 .2 .2;  .2 .7 .2;   0  0 .6;  .3 .3 .3;
           0 .9 .4;   0 .4 .9;  .9 .4  0;  .4 .9  0;  .9  0 .4;
          .4  0 .9;  .8 .5 .5;  .5 .8 .5;  .5 .5 .8;];
posFIG  = [];
AXISp   = [];
o4      = ones(1,4);
idPOS   = o4;
axisPOS = o4;
POScid  = [3,-.075];
cFIGbk  = [.25 .15 .15];
cPLTbk  = [ 0   0   0 ];
cTID    = [];
cXYax   = [ 1   1   1 ];
cXYlbl  = [.7  .8  .95];
CURcDEF = [ 1   1  .50];
cCURSOR = CURcDEF;
cDELTA  = [ 1   0   0 ];
Grid = 'on';
cGRIDs  = 0;
GridSt  = 0;
GridSty = '-';
cGRID = [.13 .13 .13];
cDEF    = 0;
for k=35:99 cTRACE = [cTRACE; .75 * cTRACE(k-34,:)]; end;
cDEFAULT = cTRACE;
LabelX = 'X axis';
LabelY = 'Y axis (Left)';
LabelYr = 'Y axis (Right)';
Title = '';
NewLimit = '';
FigName = 'plt';
Xlim   = 'default';
Ylim   = 'default';
YlimR  = 'default';
Xscale = 'linear'; Xsc = 'LinX';
Yscale = 'linear'; Ysc = 'LinY';
Quiv = 0;
Qhead = [.3 .3];
Mbar = 0;
Xslide = 0;
LineSmooth = 0;
FigShow = 1;
HelpFile = '';
HelpFileR = '';
cFile = '#';
ENApre = ones(1,2);
ENAcur = ones(1,99);
DIStrace = 0;
Style  = 0;
Marker = 0;
AXr = 0;
Right = [];
DualCur = 0;
TIDcback = '';
TIDcolumn = 0;
MenuBox = [1 1 1 1 0 1 1 1 1];
TRACEid = reshape(sprintf('Line%2d',1:99),6,99)';
TRACEmk = 0;
Xstring = '';
Ystring = '';
moveCB = '';
axisCB = '';
ucreq = '';
axLink = 1;
SubPlot = 0;
SubTrace = [];
NoCursor = 0;
xViewOpt = 0;
multiCurOpt = 0;
aLp = {}; aLv = {};
aRp = {}; aRv = {};
lLp = {}; lLv = {};
lRp = {}; lRv = {};
lXp = {}; lXv = {};
posAX   = [.1429  .0933  .8329  .8819];
posPEAK = [.0071  .0441  .0250  .0343];
posVALY = [.0071  .0060  .0250  .0343];
posDEL  = [.0350  .0441  .0280  .0343];
posM    = [.0350  .0060  .0280  .0343];
posSLDR = [.0071  .0080  .1250  .0200];
posCXL  = [.1386  .0095  .0200  .0410];
posC1X  = [.1643  .0076  .1000  .0450];
posC2X  = [.2686  .0076  .1000  .0450];
posCYL  = [.7609  .0095  .0200  .0410];
posC1Y  = [.7865  .0076  .1000  .0450];
posC2Y  = [.8908  .0076  .1000  .0450];
posAX2  = [.0070  .0870  .0580  .0000];

fontsz = (196-get(0,'screenpix'))/10;
if sum(lower(varargin{1})) == 310
     FIG = varargin{2};  figure(FIG);  set(FIG,'vis','of');
else FIG = figure('menu','none','numberT','off','back','off','vis','of');
end;
if Mver < 8.4  FIGt = FIG;
else           FIGt = get(FIG,'number');
end;
set(FIG,'PaperPositionMode','auto','invert','off','PaperOrient','land',...
        'PaperUnits','norm','DoubleBuf','on','Invert','on');
AX = axes('unit','nor','fontsi',fontsz,'tag','click');
Ret1 = [];  nt = 0;
kparam = [];
k = 1;
pp = 1;
while k<=Narg
  y  = varargin{k};  k=k+1;
  if ischar(y)
    if k>Narg
       disp('Error using ==> plt.  Not enough input arguments.');
       disp('For help on using plt, type "help plt"');
       eval('plt help;');
       return;
    end;
    kparam = [kparam k-1 k];
    y = lower(y);  yy = varargin{k};   k=k+1;
    pfx = zeros(1,5);
    while 1
      b = findstr(y(1),'+-<>.');
      if isempty(b) break; end;
      pfx(b) = 1;  y(1) = [];
    end;
    switch sum(y)
      case 546,     Title    = yy;
      case 442,      Xlim     = yy;
      case 443,      Ylim     = yy;
      case 557,     if y(1)=='y'  YlimR = yy;  AXr = 1; else cXYax = yy; end;
      case 632,    LabelX   = yy;
      case 633,    LabelY   = yy;
      case 747,   LabelYr  = yy;  if length(yy) AXr = 1; end;
      case 542,     Right    = yy;  if length(yy) AXr = 1; end;
      case 752,   DualCur  = yy;
      case 727,   FigName  = yy;
      case 640,    cPLTbk   = yy;
      case 614,    cFIGbk   = yy;
      case 626,    cTRACE   = yy;
      case 621,    cDELTA   = yy;
      case 654,    cXYlbl   = yy;
      case 769,   cCURSOR  = yy;
      case 420,      cTID     = yy;
      case 521,     if any(yy<0) GridEr = ERANOR; else GridEr = ERAXOR; end;
                       cGRID    = abs(yy);  cGRIDs = 1;
      case 983, GridSty = yy;        GridSt = 1;
      case 676,    if iscell(yy) Style = char(yy);   else Style = yy;   end;
      case 757,   if iscell(yy) Marker = char(yy);  else Marker = yy;  end;
      case 732,   if iscell(yy) TRACEid = char(yy); else TRACEid = yy; end;
      case 743,   if length(yy)==1 & yy  yy = [yy (yy+.9)/2 .9]; end;
                       TRACEmk = yy;
      case 638,    ENAcur   = yy;
      case 847,  DIStrace = yy;
      case {338 885},
                       if ~yy(3) yy(3) = round(yy(4)/.944); elseif ~yy(4) yy(4) = round(yy(3)*.944); end;
                       posFIG   = yy;
      case 241,        [rw cl] = size(yy);
                       if rw==1 & cl==4 AXISp = [0 yy];
                       elseif cl~=5 disp('Warning: xy parameter ignored. Must have 5 columns');
                       else AXISp = yy;
                       end;
      case 775,   axisPOS  = yy(1:4);
                       switch length(yy) case 5, idPOS(3)=yy(5); case 8, idPOS=yy(5:8); end;
      case 821,  TIDcback = yy;
      case 975, TIDcolumn = yy;
      case 783,   Xstring  = yy;
      case 784,   Ystring  = yy;
      case 873,  NewLimit = yy;
      case 635,    ENApre   = yy;
      case 668,    Quiv     = yy;  if min(yy)<2 disp('No quiver tail position '); return; end;
      case 515,     Qhead    = yy;
      case 841,  HelpFile = yy;
      case 955, HelpFileR = yy;
      case 959, cFile    = yy;
      case 846,  cPLTbk = get(0,'defaultaxescolor');
                       cFIGbk = get(0,'defaultfigurecolor');
                       cXYax  = get(0,'defaultaxesxcolor');  cXYlbl = cXYax;
                       cTRACE = get(0,'defaultaxescolororder');
                       if length(yy(1,:))==3 cTRACE = [cTRACE; yy]; end;
                       cDEF = 1;
      case 636,    moveCB = yy;
      case 634,    axisCB = yy;
      case 867,  axLink = yy;
      case 430,      if Mver < 8.4 FIGt = yy;
                       else          FIGt = get(yy,'Number');
                       end;
      case 862,  ucreq = yy;
      case 777,   SubPlot = yy;
      case 857,  SubTrace = yy;  ENApre(2) = 0;
      case 1084, setappdata(gcf,'MotionEdit',yy);
      case 1115, setappdata(gcf,'MotionZoom',yy);
      case 310,
      case 780,   kq = 0;
        while kq < length(yy)
          kq = kq + 1;
          switch yy(kq)
            case 'T', Grid = 'off';
            case 'M', Mbar = 1;
            case 'X', Xsc = 'LogX';  Xscale = 'Log';
            case 'Y', Ysc = 'LogY';  Yscale = 'Log';
            case 'S', Xslide = 1;
            case 'L', LineSmooth = 1;
            case 'N', NoCursor = 1;
            case 'H', FigShow = 0;
            case 'V', xViewOpt = 1;
            case 'C', multiCurOpt = 1;
            case '-', kq = kq + 1;  km =findstr(yy(kq),'HXYGPFMZRA');
                      if length(km) if km==10 MenuBox=0; else MenuBox(km)=0; end; end;
            case '+', kq = kq + 1;  km =findstr(yy(kq),'HXYGPFMZRA');
                      if length(km) if km==10 MenuBox=ones(1,9);
                                    else if pp & km~=5 MenuBox=0; pp=0; end;
                                         MenuBox(km)=1;
                                    end;
                      end;
          end;
        end;
      otherwise,
         if sum(pfx)
           if pfx(1) aLp = [aLp {y}]; aLv = [aLv {yy}]; end;
           if pfx(2) aRp = [aRp {y}]; aRv = [aRv {yy}]; end;
           if pfx(3) lLp = [lLp {y}]; lLv = [lLv {yy}]; end;
           if pfx(4) lRp = [lRp {y}]; lRv = [lRv {yy}]; end;
           if pfx(5) lXp = [lXp {y}]; lXv = [lXv {yy}]; end;
         else
           if iscell(yy)
             if length(yy)==length(Ret1)
               y = {y};
               yy = yy(:);
             else fprintf('Warning: For parameter %s, found %d elements but expected %d\n',y,length(yy),length(Ret1));
             end;
           end;
           set(Ret1,y,yy);
         end;
    end;
  else
     if k<=Narg           yy = varargin{k}; else yy = 'a'; end;
     if ~isreal(y)        yy = imag(y);  y = real(y);
     elseif isnumeric(yy) k = k+1;
     else                 yy=y; y=1:length(yy);
     end;
     nxt = length(Ret1) + 1;
     if length(find(Quiv==nxt))
        a = Ret1(Qtail);
        x0 = get(a,'x');  x0 = transpose(x0(:));
        y0 = get(a,'y');  y0 = transpose(y0(:));
        x1 = transpose(y(:));
        y1 = transpose(yy(:));
        rNaN = repmat(NaN,size(y1));
        x01 = x0+x1;  y01 = y0+y1;
        a = x01-Qhead(1)*(x1+Qhead(2)*y1);
        b = x01-Qhead(1)*(x1-Qhead(2)*y1);
        y  = [x0; x01; rNaN; a; x01; b; rNaN];
        a = y01-Qhead(1)*(y1-Qhead(2)*x1);
        b = y01-Qhead(1)*(y1+Qhead(2)*x1);
        yy = [y0; y01; rNaN; a; y01; b; rNaN];
        y = y(:);  yy = yy(:);
        ENAcur(nxt) = 0;
        nxt = 0;
     end;
     H = line(y,yy);  Ret1 = [Ret1; H];   nt = nt + length(H);
     if nxt Qtail = nt; end;
  end;
end;

if isempty(cTID) cTID = cPLTbk; end;
if iscell(Ylim) YlimR = Ylim{2}; Ylim = Ylim{1}; AXr = 1; end;
nSP = length(SubPlot) - 1;

if (iscell(LabelX) | iscell(Xlim)) & ~nSP & nt>1 SubPlot = [100 -50 100]; nSP = 2; end;
if iscell(LabelX) LabelXr = LabelX{2};  LabelX = LabelX{1}; else LabelXr = 'X axis'; end;
if iscell(Xlim)   Xlimr   = Xlim{2};    Xlim   = Xlim{1};   else Xlimr   = 'default'; end;

nSPl = nSP;
nSPr = 0;
SPw = 1;
if nSP
  SubPlot = SubPlot/100;
  k = find(SubPlot<0);
  if length(k) SPdx = -100*SubPlot(k);   SPw = floor(SPdx);
               SPdx = SPdx - SPw;  SPw = .01*SPw;
               if SPdx>.5 SPdx = SPdx-1; end;
               SubPlot(k) = [];      nSP = nSP-1;
               nSPr = nSP + 2 - k;   nSPl = nSP - nSPr;
  end;
end;

if isempty(posFIG)
  posFIG = [9 45 700 525];
  if nSPr           posFIG(3) = 840; end;
  if nSPl | nSPr>1  posFIG(4) = 575; end;
end;

nID = nt - nSP * isempty(SubTrace);
if nID>99 & TRACEid disp(sprintf('Max # of traceIDs = %d',99)); return; end;
if ~iscell(LabelY) LabelY = {LabelY}; end;
if length(LabelY) >= nSP+2  LabelYr = LabelY{nSP+2}; AXr = 1; end;
setappdata(FIG,'params',varargin(kparam));
k = sum(TIDcolumn);
if k  TIDcolumn = [nID-k TIDcolumn]; else TIDcolumn = nID; end;
ntid = max(TIDcolumn);
ncol = length(TIDcolumn);
if ncol>1 & all(axisPOS == o4) & all(idPOS == o4)
  idPOS(3) = ncol;
  axisPOS = [.4 + ncol/2, 1, (210-11*ncol-ncol^2)/200, 1];
end;
if all(idPOS == o4) & (nt<6 | isempty(LabelY{1})) & ~nSP
   idPOS(3) = 1.2;
end;
fsep = length(findstr(cFile,filesep));
if length(cFile) & ~fsep
  foobar = 0;
  if exist('foobar')
       m = feval('dbstack');
       if length(m)
          n = m(end).name;
          nq = findstr('(',n);
          if length(nq) n = n(1:nq(1)-2); end;
          np = findstr(filesep,n);
          if length(np) np=np(end); else np=0; end;
          if length(n)-np>30 & length(m)>1  n = m(end-1).name; end;
          m = feval('which',n);
       end;
  else m = GetExe;
  end;
  [pth, name] = fileparts(m);
  if cFile(1) == '#' cFile = [name 'Color']; end;
  cFile = fullfile(pth,[cFile '.mat']);
end;
if fsep cFile = [cFile '.mat']; end;
if exist(cFile) == 2
  load(cFile);
end;
if length(ENAcur)<nt  ENAcur = [ENAcur ones(1,nt-length(ENAcur))]; end;
if length(DIStrace)<nt  DIStrace = [DIStrace zeros(1,nt-length(DIStrace))]; end;
nC = length(cTRACE(:,1));
if length(Title)
   axisPOS = axisPOS .* [1 1 1 .96];
   ntx = findstr('[TexOff]',Title);
   if length(ntx)==1
        Title(ntx:ntx+7) = [];
        title(Title,'color',cXYlbl,'handlev','on','interp','none');
   else title(Title,'color',cXYlbl,'handlev','on');
   end;
end;
setappdata(AX,'DualCur',DualCur);
if ~cGRIDs
  if AXr GridEr = ERAXOR; else GridEr = ERANOR; end;
end;
if Mver >= 8.4 & AXr & ~cGRIDs & ~GridSt
  GridSty = ':';
  cGRID = [.26 .26 .26];
end;
if AXr & nID>1
   if isempty(Right) Right = nID; end;
   axisPOS = axisPOS .* [1 1 .937 1];
   if ischar(YlimR)
     mn = inf;  mx = -inf;
     if max(Right)>length(Ret1) disp('Error: Right trace number points to non-existing data'); return; end;
     for k=Right
       y = get(Ret1(k),'y');  mn = min(mn,min(y));  mx = max(mx,max(y));
     end;
      df=(mx-mn)/20; YlimR=[mn-df mx+df];
      if ~diff(YlimR) YlimR = [mn mn+max(1e-12,mn*1e-12)];  end;
   end;
   if length(Right)==1 yclr = cTRACE(mod(Right-1,nC)+1,:); else yclr = cXYax; end;
   AXr = axes('unit','nor','fontsi',fontsz,'YAxisLoc','right','ylim',YlimR,...
              'color',cPLTbk,'xcol',cXYax,'ycol',yclr,'xtick',[]);
   if ~axLink LabelYr = ['\div ' LabelYr ' \div']; end;
   ylabel(LabelYr,'color',yclr,'handlev','on','buttond','plt click link;');
   set(Ret1(Right),'par',AXr);
   setappdata(AXr,'Lhandles',Ret1(Right));
   axes(AX);
   set(gcf,'vis','of');
   AXrl = [AX AXr];
else AXrl = AX;  AXr = []; set(AX,'Box','On');
end;
mrk = repmat('+',1,nt);
mrk(Right) = 'o';
ceq = isequal(cCURSOR,CURcDEF);  cEXPbox = cCURSOR;
if ceq & sum(cPLTbk)>2 cEXPbox=1-cEXPbox; cCURSOR=1-cCURSOR; end;
curclr = [.7 .7 .7; 0 0 0; cEXPbox; cDELTA];
ENAcurS = sum(ENAcur(1:nt));
if ceq & ENAcurS>1 cCURSOR = [0 0 0]; end;
for k=1:nt  set(Ret1(k),'color',cTRACE(mod(k-1,nC)+1,:));
            if ENAcur(k) curclr=[curclr; cCURSOR]; else  set(Ret1(k),'tag','SkipCur'); end;
end;
if Style  if length(Style(:,1)) < nt Style=Style'; end;
          for k=1:nt
            if length(findstr(Style(k,1),'+o*.xsd^v<>ph'))
                  set(Ret1(k),'linest','none','marker',Style(k,1));
            else  set(Ret1(k),'linest',Style(k,:));  end;
          end;
end;
if Marker if length(Marker(:,1)) < nt Marker=Marker'; end;
          for k=1:nt set(Ret1(k),'marker',Marker(k,:)); end;
end;
for k=1:nt if DIStrace(k) set(Ret1(k),'vis','of'); end; end;

Ret1a = ones(size(Ret1));
axS = [];
kAX0 = []; if length(AXISp) kAX0 = find(~AXISp(:,1)); end;
if length(kAX0) p0 = AXISp(kAX0(1),2:5);  AXISp(kAX0,:) = []; else p0 = []; end;
if ~nSP
  if length(p0) posAX = p0;
  else posAX = posAX .* axisPOS;
  end;
else
  if length(p0) posAX = p0;
  else posAX = (posAX + [0 .0167 0 -.0167]) .* axisPOS;
  end;
  ySP = posAX(2);  ySPr = ySP;  vSpace = posAX(4);
  ySP  = cumsum([ySP   vSpace*SubPlot(1:nSPl+1)]);
  hSP  = diff(ySP)  - .03;  ySP  = ySP  + .012;
  hSpace = posAX(3);
  if nSPr hdx = .08 + SPdx; else hdx = 0; end;
  hSpace = hSpace - hdx;
  posAX = [posAX(1) ySP(1) hSpace*SPw hSP(1)];
  dx = posFIG(3);  dx1 = 64/dx;  dx2 = 14/dx;  dx3 = 67/dx;  dx4 = 19/dx;
  dy = posFIG(4);  dy1 = 20/dy;  y1a = .0076;  y2a = y1a + 24/dy;
  x1 = .1386;  x2 = x1+dx4;  x3 = x2+dx3;  x4 = x2+nSPl*dx3;
  posCXL = [x1 y2a dx2 dy1];
  posCYL = [x1 y1a dx2 dy1];
  posC1X = [x2 y2a dx1 dy1];
  posC2X = [x3 y2a dx1 dy1];
  posC1Y = [x4 y1a dx1 dy1];
  posC2Y = posC1Y + [dx3 0 0 0];
  if nSPr
    ySPr = cumsum([ySPr vSpace*SubPlot(nSPl+2:end)]);
    hSPr = diff(ySPr) - .03;  ySPr = ySPr + .012;
    SPx = sum(posAX([1 3])) + hdx + .05*length(AXr);
    x1 = SPx-.03;  x2 = x1+dx4;  x3 = x2+dx3;  x4 = x2+(nSPr-1)*dx3;
    posCXLr = [x1 y2a dx2 dy1];
    posCYLr = [x1 y1a dx2 dy1];
    posC1Xr = [x2 y2a dx1 dy1];
    posC2Xr = [x3 y2a dx1 dy1];
    posC1Yr = [x4 y1a dx1 dy1];
    posC2Yr = posC1Yr + [dx3 0 0 0];
    p2r = [posCXLr;posCYLr;posC1Xr;posC2Xr;posC1Yr;posC2Yr];
  end;
  for k=1:nSP
    l = Ret1(nt+k-nSP);
    Ret1a(nt+k-nSP) = 0;
    if isempty(SubTrace) c = get(l,'color'); else c = cXYax; end;
    a = axes('ycol',c,'xcol',c);
    setappdata(a,'Lhandles',l);
    if length(SubTrace) setappdata(a,'subTr',1); end;
    set(l,'par',a);
    if length(LabelY)>k ylabel(LabelY{k+1}); end;
    if k == nSPl+1 xlabel(LabelXr); end;
    axS = [axS a];
  end;
end;
set(AXrl,'pos',posAX,'Xscale',Xscale,'Yscale',Yscale);
if ischar(Ylim)
   Ylim = get(AX,'ylim');
end;
Ret1a = Ret1(logical(Ret1a));
setappdata(AX,'Lhandles',Ret1a);
setappdata(FIG,'Lhandles',Ret1);
axData = [AX axS AXr];
setappdata(FIG,'axis',axData);

pb = [posPEAK;posVALY;posM;posDEL];
if Xslide
   posAX2(2) = posAX2(2) + .028;  pb(:,2) = pb(:,2) + .028;  pb = [pb; posSLDR];
end;
CurID = plt('cursor',AXrl,'init',...
  [posCXL;posCYL;posC1X;posC2X;posC1Y;posC2Y;pb],...
  curclr,'', mrk, 0.8*fontsz,'','on',[],NewLimit);

set(findobj(gcf,'style','push','str','D'),'user',cDEFAULT,'tag',cFile);
CurIDstr = {@plt 'cursor' CurID};
set(AXrl,'user',CurID);
Left = setdiff(1:nt-nSP,Right);
nLeft = length(Left);
if (length(AXr) | nSP) & nLeft==1
  yclr = cTRACE(mod(Left-1,nC)+1,:); leftclr=yclr;
else yclr = cXYax;  leftclr = cXYlbl;
end;
set(AX,'xcol',cXYax,'ycol',yclr);
if ischar(Xlim)
  if nLeft Xlim = get(AX,'xlim'); else Xlim = get(AXr,'xlim'); end;
end;
[prefix xmult] = plt('metricp',max(abs((Xlim))));
if ENApre(1) & xmult~=1
  for k=1:nt set(Ret1(k),'x',xmult*get(Ret1(k),'x')); end;
  LabelX = [prefix LabelX];
else xmult = 1;
end;
set(AXrl,'xlim',Xlim*xmult);
[prefix mult] = plt('metricp',max(abs((Ylim))));
ymult = ones(1,nt);
if ENApre(2) & mult~=1
  for k=1:nLeft
     kk = Ret1(Left(k));
     set(kk,'y',mult*get(kk,'y'));
     ymult(Left(k)) = mult;
  end;
  LabelY{1} = [prefix LabelY{1}];
else mult = 1;
end;
set(AX,'ylim',Ylim*mult);
setappdata(FIG,'xymult',[xmult ymult]);
xlabel(LabelX,'color',cXYlbl,'handlev','on');
hYlab = ylabel(LabelY{1},'color',leftclr,'handlev','on');
plt('cursor',CurID,'set','moveCB2',[CurIDstr {'MVcur'}]);
if cDEF & ~cGRIDs & sum(cPLTbk)>1.5 & (Mver>=8.4 | isempty(AXr)) cGRID = 1-cGRID; end;
plt('grid',AX,'init',cGRID,GridEr,GridSty);
set(AX,'TickLen',(1-plt('grid',AX,Grid))*[.01 .025]);
set(gcf,'vis','of');
axes('pos',posCYL+[0 0 .2 0],'vis','of');
txt = text(-.02,.45,'','fontsi',fontsz,'horiz','right','buttond','plt click RMS;','user','idcur');
if length(Xstring)
  if ischar(Xstring) & Xstring(1) == '?'
        Xstring(1)=[];
        a = uicontrol('sty','edit','unit','nor','pos',posC2X.*[1 1 1.7 1],'horiz','cent',...
                 'backg',[.2 .2 .2],'foreg',[1 1 .3]);
  else  a = text(-2.22,.45,'','color',cXYlbl);
  end;
  set(a,'fontsi',fontsz,'tag','xstr');
  setappdata(a,'evl',Xstring);
  uistack(a,'bottom');
end;
if length(Ystring)
  if ischar(Ystring) & Ystring(1) == '?'
        Ystring(1)=[];
        a = uicontrol('sty','edit','unit','nor','pos',posC2Y,'horiz','cent',...
                 'backg',[.2 .2 .2],'foreg',[1 1 .3]);
  else  a = text(.6,.45,'','color',cXYlbl);
  end;
  set(a,'fontsi',fontsz,'tag','ystr');
  setappdata(a,'evl',Ystring);
  uistack(a,'bottom');
end;
nMenu = sum(MenuBox);
aid = -1;
ahi = .035*nMenu;
if nID>1 & TRACEid
  h = 19*ntid;
  ahip = ahi * 525;
  hr = 440 / (h+ahip);
  if hr<1 ahi = ahi*hr;  h = h*hr;  end;
  aidp = idPOS.*[3 521-h 50 h]./[700 525 700 525];
  aid = axes('xlim',[0 ncol],'ylim',[-ntid 0]-.5,'color',cTID,'xcol',cPLTbk,'ycol',cPLTbk,'XtickLabel',' ','YtickLabel',' ','TickLen',[0 0]','user','TraceID',...
        'unit','nor','pos',aidp,'buttond','plt misc tidmv;');
  setappdata(aid,'TIDcback',TIDcback);
  setappdata(aid,'ty',' xy');
  cRid = cPLTbk + .16*(2*(cPLTbk<.5)-1);
  row = 1;  col = 1;
  bln = 0;
  lpl = {'color'; 'marker'; 'linest'; 'linewidth'};
  [tn tw] = size(TRACEid);
  if nID>tn TRACEid = [TRACEid; repmat(' ',nID-tn,tw)]; end;
  enap = 1;
  wmax = 0;
  for k=1:nID
    s = TRACEid(k,:);
    wmax = max(wmax,length(s));
    if all(s==' ') bln=bln+1; continue; end;
    if s(1)==']' s = s(2:end); enap = 0; end;
    isR = length(find(k==Right));
    if isR & enap cR=cRid; line(col-[1 .05],[0 0]-row,'color',cR,'LineWidth',9);
    else          cR=cPLTbk;
    end;
    d = text(col-.93,-row,s);
    ms = {@plt 'click' [Ret1(k) d]};
    set(d,'fontsi',fontsz,'fontw','bol','color',cTRACE(mod(k-1,nC)+1,:),'buttond',ms);
    setappdata(d,'ty','-');
    if TRACEmk
      mk = line(TRACEmk,TRACEmk*0-row,lpl,get(Ret1(k),lpl),'buttond',ms);
      setappdata(d,'mk',mk);
      if TRACEmk(1)<.25 set(d,'color',cR); end;
    else mk = [];
    end;
    if DIStrace(k) set(d,'fonta','ita','fontw','nor'); set(mk,lpl,{[0 .3 .3] 'none' '-' 9}); end;
    row = row+1;
    if row>TIDcolumn(col) col=col+1; row=1; end;
  end;
  if bln & ncol==1
     dy = aidp(4) * (1 - (nID-bln)/nID);
     aidp = aidp + [0 dy 0 -dy];
     set(aid,'ylim',[bln-nID 0]-.5,'pos',aidp);
  end;
  if wmax>6 & ncol==1
    if aidp(3)*posFIG(3) < 55
       aidp(3) = aidp(3)*1.13; set(aid,'pos',aidp);
    end;
  end;
end;
if nMenu
  posAX2(4) = ahi;  c = .2*cXYax + .8*cFIGbk;
  cb = sum(cXYlbl)<1.5;  cb = .75*cFIGbk + .25*[cb cb cb];
  amb = axes('unit','nor','pos',posAX2,'ylim',[-nMenu 0]-.45,'Box','On','XaxisLoc','top',...
           'color',cb,'xcol',c,'ycol',c,'XtickLabel',' ','YtickLabel',' ','TickLen',[0 0]','tag','MenuBox','buttond','plt misc tidmv;');
  setappdata(amb,'ty',' xy');
  b=0;  t=[];
  btn = [CurIDstr {'scale' 'old'}];
  txx = { 'Help',    Xsc,       Ysc,       'Grid',    'Print',   'Menu',    'Mark',    'Zout',   'XY\leftrightarrow'};
  btn = { 'plt help 1;'; 'plt click TGLlogx;'; 'plt click TGLlogy;'; 'plt click TGLgrid;'; 'plt(''hcpy'',''init'',gcf);'; 'plt click TGLmenu;'; 'plt click mark;'; 'plt click ZoomOut;'; btn };
  for k=1:length(MenuBox)
    if MenuBox(k)
       b=b-1;  te = text(.5,b,txx{k},'interp','tex');
       t = [t te];  set(te,'buttond',btn{k});  setappdata(te,'ty','-');
       if k==1 set(te,'user',HelpFile,'tag',HelpFileR); end;
    end;
  end;
  set(t,'fontsi',fontsz,'color',cXYlbl,'horiz','cent');
else amb = [];
end;
CurMain = getappdata(0,'CurMain');
Hc = get(CurMain(CurID),'user');
set(Hc(4),'buttond','plt click EDIT 1;');
set(Hc(6),'buttond','plt click EDIT 2;');
if posFIG(1)<0 posFIG = abs(posFIG);
else for k=flipud(findobj('type','fig'))'
      if get(k,'pos')==posFIG  posFIG = posFIG + [30 25 0 0]; end;
     end;
end;
set(FIG,'pos',posFIG,'Name',FigName,'color',cFIGbk,'tag',sprintf('%d',FIGt),'CloseReq',{@plt 'misc' 'close' CurID});
setappdata(FIG,'ucreq',ucreq);
if Mbar plt click TGLmenu; end;
set(AX,aLp,aLv);
set(get(AX,'ylab'),lLp,lLv);
set(get(AX,'xlab'),lXp,lXv);
if length(AXr) set(AXr,aRp,aRv);
       set(get(AXr,'ylab'),lRp,lRv);
end;
setappdata(AX,'xstr',findobj(gcf,'tag','xstr'));
setappdata(AX,'ystr',findobj(gcf,'tag','ystr'));
v = 'off';
for k=1:nLeft
  if strcmp(get(Ret1(Left(k)),'vis'),'on') v = 'on'; break; end;
end;
set(hYlab,'vis',v,'ui',uicontextmenu('call','plt hideCur;'));
leftC = cPLTbk;
if length(AXr)
  ls = Ret1(Right)';
  icr = 15+Right;  icg = find(icr>length(Hc));
  if length(icg)
    if length(SubTrace) icr(icg) = [];
    else disp('Error: Subplot data must follow all main plot data in the argument list'); return;
    end;
  end;
  set(Hc(icr),'par',AXr,ERAS,ERAXOR);
  v = 'off'; gridXOR = cPLTbk;
  for k=1:length(ls)
    if strcmp(get(ls(k),'vis'),'on')
       v = 'on';
       leftC = 'none';
       gridXOR = get(gcf,'color');
       break;
    end;
  end;
  set(AXr,'vis',v);
  if GridEr(1)=='x'
    set(findobj(AX,'user','grid'),'color',bitxor(round(255*cGRID),round(255*gridXOR))/255);
  end;
end;
plt('cursor',CurID,'MVcur');
if length(moveCB) plt('cursor',CurID,'set','moveCB',moveCB); end;
if length(axisCB) plt('cursor',CurID,'set','axisCB',axisCB); end;
set(AX,'color',leftC);
cidS = CurID;
if nSP
  if ishandle(aid) & all(idPOS([1 2 4]) == [1 1 1])
    set(aid,'unit','nor');
    p = get(aid,'pos');
    p(2) = ySP(2)-p(4)-.025;
    if nMenu & (p(2) < sum(posAX2([2 4]))+.015)
       p(2) = posAX2(2) + .025;  posAX2(2) = p(2) + p(4) + .02;
       set(amb,'pos',posAX2);
    end;
    set(aid,'pos',p);
  end;
  h = plt('cursor',CurID,'get','obj');
  q = [-1 .1 .1 .1];  u = [q;q;q;q];
  py1 = get(h(5),'pos');   p1 = [u;py1;u;q];
  p = {'pos' 'xlim' 'color' 'Xscale' 'FontSize' 'TickLength' 'Box'};
  v = get(AX,p);
  if ischar(v{3}) v{3} = get(AXr,'color'); end;
  for k=1:nSPl
    v{1}([2 4]) = [ySP(k+1) hSP(k+1)];
    set(axS(k),p,v);
    p1(5,1) = p1(5,1) - dx3;
    cidS = [cidS plt('cursor',axS(k),'init',p1,'','','+',8)];
    set(axS(k),'user',cidS(end));
    CurMain = getappdata(0,'CurMain');  Hc = get(CurMain(cidS(end)),'user');
    set(Hc(6),'buttond','plt click EDIT 2;');
  end;
  if nSPr
    v{1}([1 3]) = [SPx hSpace*(1-SPw)];  p2r = [p2r;u];
    if ischar(Xlimr) k = get(Ret1(end),'x');  v{2} = [min(k) max(k)];
    else             v{2} = Xlimr;
    end;
    for k=1:nSPr
      kr = k+nSPl;
      v{1}([2 4]) = [ySPr(k) hSPr(k)];
      set(axS(kr),p,v);
      cidS = [cidS plt('cursor',axS(kr),'init',p2r,'','','+',8)];
      set(axS(kr),'user',cidS(end));
      CurMain = getappdata(0,'CurMain');  Hc = get(CurMain(cidS(end)),'user');
      if k==1 p2r([1 2 3 4 6],:) = [u;q];
              set(Hc(4),'buttond','plt click EDIT 1;');
      end;
      p2r(5,1) = p2r(5,1) - dx3;
      set(Hc(6),'buttond','plt click EDIT 2;');
     end;
  end;
  creq = [];
  s = 'plt(''cursor'',';
  cidSS = {[s int2str(cidS(1))]};
  for k=1:nSP
    SS = [s int2str(cidS(k+1))];
    cidSS = [cidSS {SS}];
    creq = [creq SS ',''clear'');'];
    plt('grid',axS(k),'init',cGRID,GridEr,GridSty);
    if length(Grid)==3 plt('grid',axS(k),'off'); end;
  end;
  set(gcf,'vis','of');
  n = length(SubTrace);   m = length(Ret1);  a = length(axData);
  if n == a
    h = 1;
    for k = 1:n
       h2 = h + SubTrace(k) - 1;
       if h2 > m break; end;
       set(Ret1(h:h2),'parent',axData(k));
       h = h2 + 1;
    end;
  elseif n == m
    for k=1:n
      if SubTrace(k) <= a  set(Ret1(k),'parent',axData(SubTrace(k))); end;
    end;
  end;

  setappdata(FIG,'creq',creq);
  set(findobj(FIG,'user','idcur'),'vis','of');
  setappdata(FIG,'c',0);
  s1a = 'c=getappdata(gcf,''c''); if c==%d setappdata(gcf,''c'',0); else setappdata(gcf,''c'',c+1);';
  s1 = sprintf(s1a,nSPl);
  s2 = ',''set'',''activeLine'',0,plt2nd({''cursor'',%d,''get'',''pos''})); end';
  s3 = ',"set","xlim",get(gca,"xlim")); end';
  for k=0:nSPl
    if k==nSPl j=1; else j=k+2; end;
    CI = cidS(k+1);
    plt('cursor',CI,'set','moveCB2',[s1 cidSS{j} sprintf(s2,CI)]);
    plt('cursor',CI,'set','axisCB',[s1 cidSS{j} s3]);
  end;
  if nSPr
    s1 = sprintf(s1a,nSPr);
    for k=1:nSPr
      kr = k+nSPl+1;
      if k==nSPr j=nSPl+2; else j=kr+1; end;
      CI = cidS(kr);
      plt('cursor',CI,'set','moveCB2',[s1 cidSS{j} sprintf(s2,CI)]);
      plt('cursor',CI,'set','axisCB',[s1 cidSS{j} s3]);
    end;
    plt('cursor',cidS(nSPl+2),'update',-1);
    axS(nSPl+1) = [];
  end;
  set(axS,'XtickLabel',[]);
end;
setappdata(FIG,'cid',cidS);
if NoCursor plt('hideCur'); end;
if Grid(2)=='n' plt('grid',AX); end;
axes(AX);
set(gcf,'vis','of');

setappdata(FIG,'epopup',plt('pop','choices', ...
{'Properties';           'multiCursor';             'xView slider';  ...
     'Cancel';                                                       ...
      'Range';  'Range\leftrightarrow';  'Range\uparrow\downarrow';  ...
     'Insert'; 'Insert\leftrightarrow'; 'Insert\uparrow\downarrow';  ...
     'Modify'; 'Modify\leftrightarrow'; 'Modify\uparrow\downarrow'}, ...
    'interp','tex','visible','off','callbk','plt click Yedit;'));

setappdata(FIG,'Dedit',{cidS(1) 9 -1 43 8 .5 0 -1});
setappdata(FIG,'NewData',0);
setappdata(FIG,'EditCur',(196-get(0,'screenpix'))/10);
setappdata(FIG,'logTR',1e6);
setappdata(FIG,'snap',[100 100]);
uic = findobj(FIG,'type','uicontrol')';
setappdata(FIG,'uic',uic);
setappdata(FIG,'txt',txt);
setappdata(FIG,'axi',axData);
v = [uic axData amb txt];
for k = 1:length(v) setappdata(v(k),'ty',' xy'); end;

if length(AXISp)
  for k = 1:length(AXISp(:,1))
    p = AXISp(k,:);  m = p(1);  p = p(2:5);
    switch m
      case -2,   if nMenu set(amb,'pos',p); end;
      case -1,   if ishandle(aid) set(aid,'pos',p); end;
      otherwise, if     m<=length(axData)          set(axData(m),'pos',p);
                 elseif m>201 & m-200<=length(uic) set(uic(m-200),'pos',p);
                 elseif m==301                     set(txt,'pos',p(1:2));
                 end;
    end;
  end;
end;

if LineSmooth & exist('isprop') & isprop(Ret1(1),'LineSmoothing')
  a = get(Ret1(1),'LineSmoothing');
  if a(2)=='f' set(Ret1,'LineSmoothing','on'); end;
end;
set(findobj(gcf,'str','O'),'user',LineSmooth);
plt misc tidtop;
h = plt('cursor',CurID,'update',-1);
if xViewOpt     plt('click','Yedit',3); end;
if multiCurOpt  plt('click','Yedit',2); end;
if FigShow set(FIG,'vis','on'); drawnow;
else       set(gcf,'vis','of');
end;
end;

function t = logTicks(lim)
  a = lim(1);  b = lim(2);
  if a<=0 | b<=0  t=lim; return; end;
  ex = floor(log10(a));  p = 10^ex;
  d = floor(a/p);  t = d*p;
  while 1
    d = d+1;
    if d>9 d=1; ex=ex+1;  p=p*10; end;
    v = d*p;   t = [t v];
    if v>=b break; end;
  end;

function setlim(ax,prop,lim);
  if lim(1) <= 0 | lim(2) <= 0
    s = get(ax,[prop(1) 'scale']);
    if s(2) == 'o' lim = abs(lim(2))*[0.001 1]; end;
  end;
  set(ax,prop,lim);

function fixMark()
  markR = findobj(gcf,'tag','markR')';
  if isempty(markR) return; end;
  for l = markR
    u = get(l,'user');
    t = u{1};  ax = u{2};  axr = u{3};  ylim = u{4};  rlim = u{5};
    if ishandle(t) p = get(t,'pos');  else p = [0 0]; end;
    y = [get(l,'y') p(2)];
    y = rlim(1) + diff(rlim) * (y - ylim(1)) / diff(ylim);
    rlim = get(axr,'ylim');   ylim = get(ax,'ylim');
    y = ylim(1) + diff(ylim) * (y - rlim(1)) / diff(rlim);
    set(l,'y',y(1),'user',[u(1:3) {ylim rlim}]);
    if ishandle(t) p(2) = y(2);  set(t,'pos',p); end;
  end;

function edg(txt,clr)
  v6 = version;  v6 = (v6(1)=='6');
  if nargin>1
    if v6 if strcmp(clr,'none') return; end;
          a = get(txt,'par');
          setappdata(txt,'li',line(0,0,'par',a,'color',clr,'clip','off'));
    else  set(txt,'edge',clr);
    end;
  elseif v6  li = getappdata(txt,'li');
             if ishandle(li)
               e = get(txt,'extent');  x=e(1); y=e(2); x2=x+e(3); y2=y+e(4)*1.1;
               set(li,'x',[x x2 x2 x x],'y',[y y y2 y2 y]);
             end;
  end;

function v = auxLoc(w,h,au)
     f = findobj('type','fig','tag',get(gcf,'tag'));  nfig = length(f);
     if nargin==3
       area = 0;
       for k=1:nfig
         p = get(f(k),'pos'); p = prod(p(3:4)); if p>area area=p; gf=f(k); end;
       end;
     else gf=gcf; au=0;
     end;
     set(0,'unit','pix');  sz = get(0,'screens');
     szw = sz(3) - w - 4;
     ppos  = get(gf,'pos');
     if ppos(4)<1 sp = sz(3:4); ppos = ppos .* [sp sp]; end;
     x = min(ppos(1)+ppos(3)+9,szw) + au;
     y = ppos(2) + ppos(4) - h - 30*(nfig-1);
     v = [x y w h];

function r2 = plt2nd(v)
  [r1 r2] = plt(v{:});

function r = get8192(h,prop)
  r = get(h,prop);

function evalQ(a)
  if ischar(a)     a = strrep(a,'"','''');
                   eval(a);
  elseif iscell(a) feval(a{:});
  else             feval(a);
  end;

function evalRep(a,rep)
  if ischar(a)     for k=1:2:length(rep) a = strrep(a,rep{k},rep{k+1}); end;
                   a = strrep(a,'"','''');
                   eval(a);
  elseif iscell(a) feval(a{:});
  else             feval(a);
  end;

function r = evalRep2(a,rep)
  if ischar(a)     for k=1:2:length(rep) a = strrep(a,rep{k},rep{k+1}); end;
                   a = strrep(a,'"','''');
                   r = eval(a);
  elseif iscell(a) r = feval(a{:});
  else             r = feval(a);
  end;

function s = ctriple(val)
  s = strrep(prin(' {%3w!  }',val),' 0.','.');
  if s(1)==' ' s = s(2:end); end;
  return;

function b = isequalNaN(x,y)
  x(isnan(x)) = 123456789;
  y(isnan(y)) = 123456789;
  b = isequal(x,y);

function v = s2d(s)
   v = sscanf(s,'%f');

function v = s2i(s)
   v = sscanf(s,'%d');
