classdef axisxx < addaxis
    properties(Constant)
        dimension = 'X';
    end
    properties(Dependent=true)
        XLabel;
    end
    methods
        function obj = axisxx(varargin)
            obj = obj@addaxis(varargin{:});
        end
        function obj = plot(obj, varargin)
            obj = obj.plot@addaxis(varargin{:});
            % Which side will it be added to?
            numax = length(addaxis.getAxes(obj.h_parent));
            if mod(numax, 2)
                side = 'top';
            else
                side = 'bottom';
            end
            set(obj.h_axis,'xaxislocation',side);
            set(obj.h_axis,'ytick',[]);
            set(obj.h_axis,'box','off');
            % Resize axis in case ticks have changed size
            obj.resizeParentAxis(obj.h_parent);
        end
        function label = get.XLabel(obj)
            % Just get the current text..
            label = obj.getLabel();
        end
        function obj = set.XLabel(obj, val)
            obj = obj.setLabel(val);
        end
        function obj = applyTransform(obj, handle, event)
            % Apply linear transform to data in y dimension. Used for
            % listener callbacks
            try
                p_lim = get(obj.h_parent, 'XLim');
                c_lim = get(obj.h_axis, 'XLim');
                m = diff(p_lim)/diff(c_lim);
                c = p_lim(1) - m*c_lim(1);
                xdata = get(obj.h_axis_line, 'XData');
                xdata = xdata * m + c;
                set(obj.h_parent_line, 'XData', xdata);
            catch ex
                warning('Error in addaxis transformation fcn.\nMessage %s\n Report %s\n', ex.message, ex.getReport);
            end
        end
    end
    methods(Access=protected)
        function CentralPosition = resizePeerAxis(obj, CentralPosition)
            % CentralPosition = TightPosition of parent axis + labels,
            % ticks and other addaxis objects.
            
            % Place the axis in the desired spot, returning space available
            set(obj.h_axis, 'Units', 'Pixels');
            TForm = [-1 0 1 0; 0 -1 0 1; 0 0 1 0; 0 0 0 1]; % TI to Pos tform matrix
            TI = get(obj.h_axis, 'TightInset');
            PositionInner = get(obj.h_axis, 'Position');
            PositionTight = PositionInner + TI * TForm; % Total current dims of obj.
            
            ParentInner = get(obj.h_parent, 'Position');
            
            % Align inner x position and x height to match parent AX
            PositionInner(3) = ParentInner(3);
            PositionInner(1) = ParentInner(1);
            
            % Which side of AXIS TIGHTINSET to add width to?
            if strcmpi(get(obj.h_axis, 'XAxisLocation'), 'top')
                % Place the axis on the top of central axis
                PositionInner(2) = CentralPosition(2) + CentralPosition(4) + TI(2);
                % Increase central box to include the size of this axis
                CentralPosition(4) = CentralPosition(4) + PositionTight(4);
            else
                % Place the axis on the bottom of the central axis
                PositionInner(2) = CentralPosition(2) - TI(4) - PositionInner(4);
                % Move bottom bound and increase height of the centre, box to
                % include this axis
                CentralPosition(2) = CentralPosition(2) - PositionTight(4);
                CentralPosition(4) = CentralPosition(4) + PositionTight(4);
            end
            % Set the new position of Axis
            set(obj.h_axis, 'Position', PositionInner);
        end
        function RequiredSpace = getRequiredSpace(obj)
            % Calculate inner and outer positions required for axis
            TForm = [-1 0 1 0; 0 -1 0 1; 0 0 1 0; 0 0 0 1]; % TI to Pos tform matrix
            set(obj.h_axis, 'Units', 'Pixels');
            PositionInner = get(obj.h_axis, 'Position');
            % Here is a convenient place to set the width...
            PositionInner(4) = 10;
            set(obj.h_axis, 'Position', PositionInner);
            % Calculate space required with tight inset
            TI = get(obj.h_axis, 'TightInset');
            PositionTight = PositionInner + TI * TForm; % Total current dims of obj.
            if strcmpi(get(obj.h_axis, 'XAxisLocation'), 'top')
                RequiredSpace = [0 0 0 PositionTight(4)];
            else
                RequiredSpace = [0 PositionTight(4) 0 0];
            end
        end
    end
end
