
function aadzoomstartfcn(Fig, event)
AD = getaddaxisdata(event.Axes, 'PreZoomYLim');
if isscalar(AD) && AD == 0
    setaddaxisdata(event.Axes, get(event.Axes, 'YLim'), 'PreZoomYLim');
end