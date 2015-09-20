classdef axisyy < addaxis
    properties(Constant)
        dimension = 'Y';
    end
    properties(Dependent=true)
        YLabel;
    end
    methods
        function obj = axisyy(varargin)
            obj = obj@addaxis(varargin{:});
        end
        function obj = plot(obj, varargin)
            obj = obj.plot@addaxis(varargin{:});
            
            % Which side will it be added to?
            numax = length(addaxis.getAxes(obj.h_parent));
            if mod(numax, 2)
                side = 'right';
            else
                side = 'left';
            end
            set(obj.h_axis,'yaxislocation',side);
            set(obj.h_axis,'box','off');
            set(obj.h_axis,'xtick',[]);
            % Resize axis in case ticks have changed size
            obj.resizeParentAxis(obj.h_parent);
        end
        function label = get.YLabel(obj)
            % Just get the current text..
            label = obj.getLabel;
        end
        function obj = set.YLabel(obj, val)
            obj = obj.setLabel(val);
        end
        function obj = applyTransform(obj, handle, event)
            % Apply linear transform to data in y dimension. Used for
            % listener callbacks
            for i = 1:length(obj.h_axis_line)
                try
                    p_lim = get(obj.h_parent, 'YLim');
                    c_lim = get(obj.h_axis, 'YLim');
                    m = diff(p_lim)/diff(c_lim);
                    c = p_lim(1) - m*c_lim(1);
                    ydata = get(obj.h_axis_line(i), 'YData');
                    ydata = ydata * m + c;
                    set(obj.h_parent_line(i), 'YData', ydata);
                catch ex
                    warning('Error in addaxis transformation fcn.\nMessage %s\n Report %s\n', ex.message, ex.getReport);
                end
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
            
            % Align inner y position and y height to match parent AX
            PositionInner(4) = ParentInner(4);
            PositionInner(2) = ParentInner(2);
            
            % Which side of AXIS TIGHTINSET to add width to?
            if strcmpi(get(obj.h_axis, 'YAxisLocation'), 'right')
                % Place the axis on the right of central axis
                PositionInner(1) = CentralPosition(1) + CentralPosition(3) + TI(1);
                % Increase central box to include the size of this axis
                CentralPosition(3) = CentralPosition(3) + PositionTight(3);
            else
                % Place the axis on the left of the central axis
                PositionInner(1) = CentralPosition(1) - TI(3) - PositionInner(3);
                % Move left bound and increase width of the centra, box to
                % include this axis
                CentralPosition(1) = CentralPosition(1) - PositionTight(3);
                CentralPosition(3) = CentralPosition(3) + PositionTight(3);
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
            PositionInner(3) = 10;
            set(obj.h_axis, 'Position', PositionInner);
            % Calculate space required with tight inset
            TI = get(obj.h_axis, 'TightInset');
            PositionTight = PositionInner + TI * TForm; % Total current dims of obj.
            if strcmpi(get(obj.h_axis, 'YAxisLocation'), 'right')
                RequiredSpace = [0 0 PositionTight(3) 0];
            else
                RequiredSpace = [PositionTight(3) 0 0 0];
            end
        end
    end
end
