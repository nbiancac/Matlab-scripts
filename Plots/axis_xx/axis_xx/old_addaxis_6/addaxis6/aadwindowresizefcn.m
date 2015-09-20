% Edit by Luke Plausin, 14/9/14. Works with MATLAB R2013b, untested on
% other platforms. Original AddAxis Package by Harry Lee
function aadwindowresizefcn(hObject, eventdata)
%ADDAXIS  Window resize callback function provided for automatic GUI resize
%
%  usage:
%  This function is internal.
%
%  See also
%  ADDAXIS, ADDAXISPLOT, ADDAXISLABEL, AA_SPLOT
ax_manage = getaddaxisdata(hObject, 'ActiveWindowResizeH');
if ax_manage
    for i = 1:length(ax_manage)
        aadaxisresizefcn(ax_manage(i));
    end
end
end
