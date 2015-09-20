% Edit by Luke Plausin, 14/9/14. Works with MATLAB R2013b, untested on
% other platforms. Original AddAxis Package by Harry Lee
function hax = aadaxisresizefcn(hax, OuterPos)
%ADDAXIS  Window resize callback function provided for automatic GUI resize
%
%  usage:
%  This function is internal.
%
%  See also
%  ADDAXIS, ADDAXISPLOT, ADDAXISLABEL, AA_SPLOT

% Main input to function is hidden: OuterPos, this is the outer extent of
% the axis...
HAxUnits = get(hax, 'Units');
set(hax, 'Units', 'Pixels');

hfig = hax;
while ~strcmpi(get(hfig, 'Type'), 'figure')
    hfig = get(hfig, 'Parent'); % In case docked inside a tab/panel..
end

% Default outerpos to figure dimensions..
figunits = get(hfig, 'Units');
set(hfig, 'Units', 'pixels');
if nargin <= 1
    OuterPos = get(hfig, 'Position');
    OuterPos(1:2) = [0, 0]; % Set origin to LH corner
end

% Resize plot area incl. colorbar if present...
spacer = 2; % 2 px border
TForm = [-1 0 1 0; 0 -1 0 1; 0 0 1 0; 0 0 0 1]; % TI to Pos tform matrix
TI = get(hax, 'TightInset');

try
    % This part calculates the border needed to store extra axis labels
    % from CAN plot
    TI_TOT_AA = zeros(1,4);
    
    AAD = getaddaxisdata(hax);
    for i = 1:length(AAD)
        ADH = AAD{i};
        hg_ax = ADH(1); % The axis are a bunch of axis objects and lineseries stuck together
        % Can treat these the same way as ColorBar...
        set(hg_ax, 'Units', 'Pixels');
        AAPosInner = get(hg_ax, 'Position');
        AATI = get(hg_ax, 'TightInset');
        AAPosOuter = AAPosInner + AATI * TForm; % Total current dims of obj.
        % Which side of AXIS TIGHTINSET to add width to?
        if strcmpi(get(hg_ax, 'YAxisLocation'), 'right')
            TI_TOT_AA(3) = TI_TOT_AA(3) + AAPosOuter(3);
        else
            TI_TOT_AA(1) = TI_TOT_AA(1) + AAPosOuter(3);
        end
        % Store for later
        AAPO{i} = AAPosOuter;
        AAPI{i} = AAPosInner;
    end
catch
    TI_TOT_AA = 0;
end

% Colorbar code. This works fine
cbar = findobj(hfig, 'Tag', 'Colorbar');
TI_CB_AX= zeros(1,4);
if ~isempty(cbar)
    if length(cbar) > 1 % Can't think why there'd be two colorbars
        cbar = cbar(1);
    end
    cbunit = get(cbar, 'Units');
    set(cbar, 'Units', 'pixels');
    TICB = get(cbar, 'TightInset');
    TICB(1) = spacer;
    CBInner = get(cbar, 'Position');
    CBInner(3) = 10; % Make it thinner!
    CBOuter = CBInner + TICB * TForm;
    TI_CB_AX = [0 0 CBOuter(3) 0]; % Add width of colorbar to TightInset
end

% Transofrm borders and subtract from outerpos to set innerpos
InnerPos = OuterPos - ((TI + spacer + TI_TOT_AA + TI_CB_AX) * TForm);
InnerPos = round(InnerPos);
InnerPos(InnerPos <= 0) = 10; % Guard against negative width
set(hax, 'Position', InnerPos);

% This position is between Outer and (Inner +TI), to arrange the CAN axis
% around
InsetPos = InnerPos + ((TI) * TForm); % May need the spacer here...

% Reposition all axis labels...
try
    for i = 1:length(AAD)
        ADH = AAD{i};
        hg_ax = ADH(1); % The axis are a bunch of axis objects and lineseries stuck together
        % Can treat these the same way as ColorBar...
        AAPosOuter = AAPO{i};   % Retreive positions..
        AAPosInner = AAPI{i};
        AATI = get(hg_ax, 'TightInset');
        % Adjust inner (active) height to main AX
        AAPosInner(4) = InnerPos(4); % Adjust height equal to ax
        AAPosInner(2) = InnerPos(2); % Set YPos equal to ax
        
        % Which side of MAIN AXIS is it on? Set XPos accordingly
        if strcmpi(get(hg_ax, 'YAxisLocation'), 'right')
            AAPosInner(1) = InsetPos(1) + InsetPos(3) + AATI(1);
            InsetPos(3) = InsetPos(3) + AAPosOuter(3);
        else
            AAPosInner(1) = InsetPos(1) - AAPosOuter(3) + AATI(1);
            InsetPos(1) = InsetPos(1) - AAPosOuter(3);
            InsetPos(3) = InsetPos(3) + AAPosOuter(3);
        end
        
        % Finally, set the new position of CAN Axis
        set(hg_ax, 'Position', AAPosInner);
    end
    
    % Reposition the colorbar to a convenient place
    if ~isempty(cbar)
        CBInner(4) = InnerPos(4); % Adjust height equal to ax
        CBInner(2) = InnerPos(2); % Set YPos equal to ax
        CBInner(1) = OuterPos(1) + OuterPos(3) - CBOuter(3); %Put colorbar to far right
        %InnerPos(1) + InnerPos(3) + TI(3) - CBOuter(3) + spacer + TICB(1);
        set(cbar, 'Position', CBInner);
    end
    
    % Reset all units
    if exist('cbunit', 'var')
        set(cbar, 'Units', cbunit);
    end
    set(hax, 'Units', HAxUnits);
    set(hfig, 'Units', figunits);
    
catch ex
    warning('addaxisdata:WindowResize', 'Could not resize extra axis\nMSG: %s', ex.message);
end

end
