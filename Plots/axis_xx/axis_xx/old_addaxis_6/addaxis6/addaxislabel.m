function varargout = addaxislabel(varargin)
%ADDAXISLABEL adds axis labels to axes made with ADDAXIS.m
%
%  handle_to_text = addaxislabel(axis_number, label);
%
%  See also
%  ADDAXISPLOT, ADDAXIS, SPLOT

if isstr(varargin{1}),
    axnum = varargin{2};
    label = varargin{1};
else
    label = varargin{2};
    axnum = varargin{1};
end

%  get current axis
cah = gca;
%  axh = get(cah,'userdata');
axh = getaddaxisdata(cah,'axisdata');

%  get axis handles
axhand = cah;
postot(1,:) = get(cah,'position');
for I = 1:length(axh)
    axhand(I+1) = axh{I}(1);
    postot(I+1,:) = get(axhand(I+1),'position');
end

%  set current axis to the axis to be labeled
if ishandle(axnum) && strcmpi(get(axnum, 'Type'), 'axes');
    axhand_num = axnum;
else
    axhand_num = axhand(axnum);
end
axes(axhand_num);
htxt = ylabel(label);
set(htxt,'color',get(axhand_num,'ycolor'));

%  set current axis back to the main axis
axes(cah);

if nargout == 1
    varargout{1} = htxt;
end

% Resize function in case axes have changed width..
% Get the figure handle...
hfig = cah;
while ~strcmpi(get(hfig, 'Type'), 'figure')
    hfig = get(hfig, 'Parent'); % In case docked inside a tab/panel..
end

% Resize function in case axes have changed width..
fcn = get(hfig, 'ResizeFcn');
fcn(hfig, 'addaxislabel');

