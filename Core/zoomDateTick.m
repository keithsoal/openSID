function zoomDateTick(obj,event_obj)

nticks = 10;                             % How many tick marks to use

limits = get(event_obj.Axes,'XLim');    % Get x limits after zooming

newticks = linspace(limits(1),limits(2),nticks); % Create n ticks

set(event_obj.Axes,'XTick',newticks);   % Set x tick marks in axes

% Change format using "datetick" but preserve custom ticks: 

datetick(event_obj.Axes,'x','dd.mm HH:MM:SS','keepticks')