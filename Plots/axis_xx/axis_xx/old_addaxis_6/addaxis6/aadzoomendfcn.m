function aadzoomendfcn(Fig, event)
PreYLim = getaddaxisdata(event.Axes, 'PreZoomYLim');
if ~isscalar(PreYLim)
    PostYLim = get(event.Axes, 'YLim');
    % Propagate to all axes...
    zoompts = (PostYLim - PreYLim(1))/diff(PreYLim);
    
    AAD = getaddaxisdata(event.Axes);
    for i = 1:length(AAD)
        H = AAD{i};
        AXL = get(H(1), 'YLim');
        NYL = zoompts * diff(AXL) + AXL(1);
        set(H(1), 'YLim', NYL);
    end
    setaddaxisdata(event.Axes, 0, 'PreZoomYLim');
    % Resize function in case axes have changed width..
    fcn = get(Fig, 'ResizeFcn');
    fcn(Fig, event);
end
