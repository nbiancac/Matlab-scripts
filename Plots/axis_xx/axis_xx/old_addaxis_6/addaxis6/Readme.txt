%  Removed debugging line 11/15/2005

%  Added zoom support for version R14SP2 11/24/2005

%  Updated installation notes 11/24/2005

%  In addaxisplot, only automatically change color to the axis color, if no
%  other color specified  11/28/2005

%  changed implementation to use appdata instead of userdata. (setaddaxisdata.m, getaddaxisdata.m)  
%  Also, added addaxisreset to delete added axes and reset the parent axis
%  and addaxisset which is used to set the ylimits of any axis after plotting  12/23/05

%  changed splot.m to aa_splot.m to avoid the same name as splot.m that already 
%  exists in another toolbox.  2/15/06

%  changed addaxislabel so it doesn't matter if you use it addaxislabel(axis_number,label)
%  or addaxislabel(label,axis_number) 2/15/06

%  Included new files aadzoomendfcn, aadzoomstartfcn, aadDataTipText, 
%  aadaxisresizefcn and aadwindowresizefcn to handle MATLAB's built in zoom, 
%  pan and datatip tools in a more natural way. Works as of R2013b (if you
%  are having troubles with older versions of MATLAB, try Harry Lee's original
%  addaxis5 from the file exchange). The new version packs axis in more tightly
%  and arranges around the colorbar, so that you have lots of room left on the
%  plot. Also removed buttonupfcn callback and two zoom functions which are
%  now redundant.     Luke Plausin 16/9/2014.

Installation Notes:
-------------------------------------------

Include the AddAxis6 folder in your MATLAB search path, you can do this 
through the path tool by typing 
>> pathtool   

at the command line. Now you are ready to go. Use this command to launch an
example

>> edit addaxisdatatest.m
