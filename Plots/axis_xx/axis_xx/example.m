% Luke Plausin 2015
%
% This is an example to show you how to use the axisxx and axisyy classes
%
% Example: Multiple axis
% Create the figure, and some test functions...
figure(3); clf;
title('Test 3');
x = 0:.1:4*pi;
plot(x, sin(x));
xlabel('X1 Angle');
ylabel('Y1 Sin(x)');

radii = mod(1:11,2) + 1;
angles = linspace(0,2*pi,11) + pi/2;
star_xpts = real(radii.*exp(i*angles));
star_ypts = imag(radii.*exp(i*angles));

angles = linspace(0,2*pi,150);
circ_xpts = real(exp(i*angles));
circ_ypts = imag(exp(i*angles));

h_ax = gca;
% To plot something on a second Y axis, use the function axisyy. It takes 
% any arguments that plot(...) does, and will return a handle to the 
% additional axis object which is created.

Y2 = axisyy(x, 3*cos(x), 'r-');

% To change properties of the additional axis, you can set properties of
% the axis object. The lineseries and axis object can be manipulated
% through the h_axis and h_axis_line properties. Whenever you change YData,
% limits or anything else the axis will be reset for you.

% Example - change ydata of new line object
set(Y2.h_axis, 'YLim', [-4 4]);
set(Y2.h_axis_line, 'YData', 3*cos(x) + 1);

% A shortcut to YLabel are provided directly through the addaxis object
Y2.YLabel = 'Y2 (Multiple Lines)';

% You can draw multiple objects on the same additional axis by calling the
% plot method on the returned object.
Y2.plot(x, -(3*cos(x)+1), 'r-');

% You can also add plots which share the Y axis with the main plot. XLim
% and XLabel can be input as arguments to the addaxis function
X2 = axisxx(circ_xpts, circ_ypts, 'm-', 'XLim', [-1 1], 'XLabel', 'X2 (Circle)');

% You can have as many axes as you want in your plot. You can also use any
% function you want to plot. Here we will create a third Y axis containing
% a patch object
X3 = axisxx('PlotFcn', @patch, star_xpts/2, star_ypts/2, 'k-', 'FaceColor', 'y', 'LineWidth', 3, ...
    'FaceAlpha', 0.3, 'XLim', [-1 1], 'XLabel', 'X3 (Patch Star)');
set(X3.h_axis, 'xcolor', [0 0 1]);

% The resize configuration also supports colorbars. You may need to force a
% resize for it to display nicely
h_cb = colorbar;
ylabel(h_cb, 'Colorbar compatible!');
addaxis.resizeParentAxis(h_ax);

% If you want to remove one of your axis, just delete the object.
%delete(Y2);

% This package is compatiable with the Zoom, Pan and datatip utilities.
% Have a go!

%% To come in the next update.
% Subplots
% Logarithmic axis
% XY axis (lineseries with both axis different
% 3D plots

