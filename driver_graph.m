function driver_graph()
    % Creates main figure with open file button and save graphs button.
    
    % main window
    hFig = uifigure();
    hFig.WindowState = 'maximized';
    set(hFig, 'Color', [0.4 0.4 0.4]);

    % eforce logo
    logo = uiimage(hFig);
    logo.ImageSource = 'logo.png';
    logo.Position = [20 990 250 40];

    openFileButt = uibutton (hFig , 'Position', [1680 , 20 , 100 , 25] , ...
    'Text', 'Open file');
    openFileButt.ButtonPushedFcn = @(src, event) openFileAndProcessData(hFig);

    saveGraphsButt = uibutton (hFig , 'Position', [1800 , 20 , 100 , 25] , ...
    'Text', 'Save graphs');
    saveGraphsButt.ButtonPushedFcn = @ ( src , event ) saveGraphs(hFig);
end
 
% taken from https://www.mathworks.com/matlabcentral/answers/368925-can-i-save-uiaxes-plot-in-anyway
function myCopyObj(axFrom, axTo)
    % Custom version of copyobj() which works to copy a UIAxes object to
    % an old-school axes object.
    
    % Copy children (lines).
    copyobj(axFrom.Children, axTo);
    % Copy titles and labels.
    copyobj([axFrom.Title axFrom.XLabel axFrom.YLabel], axTo)
    
    % Transfer other properties.
    uiAxParams = get(axFrom);
    uiAxParamNames = fieldnames(uiAxParams);
    % Get list of editable params in new axis
    editableParams = fieldnames(set(axTo));
    % Remove the UIAxes params that aren't editable in the new axes (add others you don't want)
    badFields = uiAxParamNames(~ismember(uiAxParamNames, editableParams));
    badFields = [badFields; 'Parent'; 'Children'; 'XAxis'; 'YAxis'; 'ZAxis'; ...
    'Position'; 'OuterPosition'; 'XLabel'; 'YLabel'; 'ZLabel'; 'Title'];
    uiAxGoodParams = rmfield(uiAxParams,badFields);
    % set editable params on new axes
    set(axTo, uiAxGoodParams)
end

function saveGraphs(hFig)
% Saves graphs in default folder with name graph_N,
% where N is number of graph in @hFig (left to right)
    allGraphs = hFig.UserData;
    
    for n = 1 : length(allGraphs)
        f = figure('Visible', false);
        axNew = axes(f);
        % Copy
        myCopyObj(allGraphs(n), axNew);
        % Print
        print(f, sprintf('graph_%d', n), '-dpng');
        delete(f);
    end    
    
end

function updateGauge(sld, ax, x, y, z, label, hFig)
     % Updates graph in @ax, text in @label and stores graph for later
     % saving in userData in @hFig

     % Creates minimum hull and is used as source for later slicing
     d = delaunayTriangulation(x(:), y(:), z(:));
    
     % Slicing the delaunay triangulation at given speed
     px = sliceDelaunay(d, 'z', sld.Value);
     
     % The area when slicing at given speed
     x_axis = px(1,:);
     y_axis = px(2,:);
     label.Text = sprintf("Area at given slice: %.2f",polyarea(x_axis, y_axis));
     
     % Border of slice for nice plot
     k = boundary(x_axis.',y_axis.');
     
     % Set limits so data are more visually clear
     xlim(ax, [min(x), max(x)]);
     ylim(ax, [min(y), max(y)]);
     
     set(ax, 'Color', [0.2 0.2 0.2]);
     plot(ax,x_axis(k), y_axis(k));
     title(ax, sprintf('Speed: %.2f m/s', sld.Value));
     xlabel(ax,'Lateral G (m/s^2)');
     ylabel(ax,'Longitudinal G (m/s^2)');
     
     allGraphs = hFig.UserData;
     allGraphs(3) = ax;
     set(hFig, 'UserData', allGraphs);
end

function openFileAndProcessData(hFig)
% Opens dialog box to open file with .mat filename extension,
% then creates graphs from data loaded from selected file

    uiopen('*.mat');

    ax = uiaxes(hFig, 'Position', [20 450, 900 500]);
    ax2 = uiaxes(hFig, 'Position', [940 450, 900 500]);
    ax3 = uiaxes(hFig, 'Position', [20 100, 600 300]);
    %ax4 = uiaxes(hFig, 'Position', [640 60, 600 300]);
    %ax5 = uiaxes(hFig, 'Position', [1260 60, 600 300]);
    
    % TODO add filtering
    x = IMU_Acceleration_Lat; 
    y = IMU_Acceleration_Long; 
    z = a1_WspeedRR;

    % Create label, where area of slice will be displayed
    label = uilabel(hFig);
    label.Position = [700 150 800 200];
    label.FontSize = 42;
    label.Text = "";
    label.FontColor = [1 1 1];
    
    % Create slider (for changing speed at which is sliced)
    sld = uislider(hFig, ...
        'ValueChangedFcn',@(sld,event) updateGauge(sld, ax3, x, y, z, label, hFig));
    sld.Position = [20 40 600 300];
    sld.Value = 15;
    sld.Limits = [min(z) max(z)];
    
    % Trigger slider action with default value
    updateGauge(sld, ax3, x, y, z, label, hFig);
    
    % GGV diagram with points
    scatter3(ax, x, y, z, 4, z,'filled');
    set(ax, 'Color', [0.2 0.2 0.2]);
    xlabel(ax,'Lateral force (m/s^2)');
    ylabel(ax,'Longitudinal force (m/s^2)');
    zlabel(ax,'Speed [m/s]');
    title(ax,"GGv diagram");
    
    % GGV diagram as solid object created by delaunay triangulation
    DT = delaunay(x, y);
    h = trisurf(DT, x, y, z, 'Parent', ax2);
    set(ax2, 'Color', [0.2 0.2 0.2]);
    set(h, 'EdgeColor', 'none');
    xlabel(ax2, 'Lateral force (m/s^2)');
    ylabel(ax2, 'Longitudinal force (m/s^2)');
    zlabel(ax2, 'Speed [m/s]');
    title(ax2, "GGV diagram");
    
    % Save all graphs to user data in case user wants to save them as
    % images
    allGraphs = [ax, ax2];
    set(hFig, 'UserData', allGraphs);
end