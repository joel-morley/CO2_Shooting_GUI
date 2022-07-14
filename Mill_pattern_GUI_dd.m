function Mill_pattern_GUI(varargin)
    
            % Create UIFigure and hide until all components are created
            UIFigure = uifigure('Visible', 'off');
            UIFigure.Position = [100 100 742 350];
            UIFigure.Name = 'MATLAB UIFigure';
            
    milldata = guidata(UIFigure);
    
    setting_file = 'mill_settings.mat';
    load(setting_file, 'milldata');

%     milldata.freq = [4 5 6 4 6 6];
%     milldata.rad = [10 20 30 40 60 80];
%     milldata.rot = [0.05 0.1 0.15 0.2 0.25 0.3];
%     milldata.n_rings = 4;
%       milldata.ratio = 1;
%       milldata.ecc_angle = 0;
%     milldata.x2 = ones(1,6);
%     milldata.x3 = ones(1,6);
%     milldata.x4 = ones(1,6);
%     milldata.x5 = ones(1,6);
%     milldata.x6 = ones(1,6);
%     milldata.x1 = ones(1,6);
%     milldata.y2 = ones(1,6);
%     milldata.y3 = ones(1,6);
%     milldata.y4 = ones(1,6);
%     milldata.y5 = ones(1,6);
%     milldata.y6 = ones(1,6);
%     milldata.y1 = ones(1,6);

    guidata(UIFigure,milldata)


        % Button pushed function: SaveButton
        function SaveButtonPushed(UIFigure, event)
            save(setting_file, 'milldata');
        end
    
        % Value changed function: RingDropDown
        function RingDropDownValueChanged(UIFigure, event)
            value = RingDropDown.Value;
            milldata = guidata(UIFigure);
            switch value
                case {1, 2, 3, 4, 5, 6}
                   FrequencySpinner.Value = milldata.freq(value);
                   RadiusSpinner.Value = milldata.rad(value);
                   RotationSpinner.Value = milldata.rot(value);
                   milldata.ring = value;
                   guidata(UIFigure,milldata);   
            end
        end

        % Value changed function: NoringsSpinner
        function NoringsSpinnerValueChanged(UIFigure, event)
            value = NoringsSpinner.Value;
            milldata = guidata(UIFigure);
            milldata.n_rings;
            guidata(UIFigure,milldata)
            [millx, milly] = dotmillingschematic(RingDropDown.Value);
            plot_milling();            
        end

        % Value changed function: FrequencySpinner
        function FrequencySpinnerValueChanged(UIFigure, event)
            value = FrequencySpinner.Value;
            ring_num = RingDropDown.Value;
            milldata = guidata(UIFigure);
            switch ring_num
                case {1, 2, 3, 4, 5, 6}
                    milldata.freq(ring_num) = value;
                    guidata(UIFigure,milldata);
            end
            [millx, milly] = dotmillingschematic(RingDropDown.Value);
            plot_milling(); 
        end

        % Value changed function: RadiusSpinner
        function RadiusSpinnerValueChanged(UIFigure, event)
            value = RadiusSpinner.Value;
            ring_num = RingDropDown.Value;
            milldata = guidata(UIFigure);
            switch ring_num
                case {1, 2, 3, 4, 5, 6}
                    milldata.rad(ring_num) = value;
                    guidata(UIFigure,milldata);
            end
            [millx, milly] = dotmillingschematic(RingDropDown.Value);
            plot_milling();            
        end

        % Value changed function: RotationSpinner
        function RotationSpinnerValueChanged(UIFigure, event)
            value = RotationSpinner.Value;
            milldata = guidata(UIFigure);
            milldata.rot = value;
            guidata(UIFigure,milldata);
            [millx, milly] = dotmillingschematic(RingDropDown.Value);
            plot_milling();            
        end

        % Value changed function: Consecutive
        function ConsecutiveCheckBoxValueChanged(UIFigure, event)
            value = ConsecutiveCheckBox.Value;
            milldata = guidata(UIFigure);
            milldata.consecutive = value;
            guidata(UIFigure,milldata);
            [millx, milly] = dotmillingschematic(RingDropDown.Value);
            plot_milling();            
        end

        % Value changed function: Radial
        function RadialCheckBoxValueChanged(UIFigure, event)
            value = RadialCheckBox.Value;
            milldata = guidata(UIFigure);
            milldata.radial = value;
            guidata(UIFigure,milldata);
            [millx, milly] = dotmillingschematic(RingDropDown.Value);
            plot_milling();            
        end

        % Value changed function: clockwise
        function ClockwiseCheckBoxValueChanged(UIFigure, event)
            value = ClockwiseCheckBox.Value;
            milldata = guidata(UIFigure);
            milldata.clockwise = value;
            guidata(UIFigure,milldata);
            [millx, milly] = dotmillingschematic(RingDropDown.Value);
            plot_milling();             
        end

        % Value changed function: RatioSpinner
        function RatioSpinnerValueChanged(UIFigure, event)
            value = RatioSpinner.Value;
            milldata = guidata(UIFigure);
            milldata.ratio = value;
            guidata(UIFigure,milldata);
            for idx = 1:RingDropDown.Value
                [millx, milly] = dotmillingschematic(idx);
                plot_milling();
            end             
        end

        % Value changed function: AngleSpinner
        function AngleSpinnerValueChanged(UIFigure, event)
            value = AngleSpinner.Value;
            milldata = guidata(UIFigure);
            milldata.ecc_angle = value;
            guidata(UIFigure,milldata);
            for idx = 1:RingDropDown.Value
                [millx, milly] = dotmillingschematic(idx);
                plot_milling();
            end
        end

    %% Helper functions
    
    function [xp, yp] = dotmillingschematic(ring_num)
    milldata = guidata(UIFigure);
    
    for j = 1:milldata.n_rings
        
    r = milldata.rad(j); % Fiber radius in mm
    deg_of_ecc = milldata.ratio;
    phi = milldata.ecc_angle; % angle of ellipse
    rx = r.*deg_of_ecc;
    ry = r;
    
    r = @(theta) (rx*ry*(r_list(r_idx)^2))./(sqrt(((rx*r_list(j))^2)*sin((phi*pi/180)+theta).^2+((ry*r_list(j))^2)*cos((phi*pi/180)+theta).^2)); % Distance from origin of point (radius) based on the angle of rotation about the origin (constant if no eccentricity)
    theta(j,:)=linspace(0+(milldata.rot(j)*j),...% Angle of rotation about origin required between each point on ring
                        (2*pi+(milldata.rot(j)*j))-((2*pi+(milldata.rot(j)*j))/milldata.freq(j)),...
                        milldata.freq(j));
    xp = (rx.*cos(theta(j,:)).*cos(phi*pi/180))-(ry.*sin(theta(j,:)).*sin(phi*pi/180));
    yp = (rx.*cos(theta(j,:)).*sin(phi*pi/180))+(ry.*sin(theta(j,:)).*cos(phi*pi/180));
    milldata.x_pos{j} = xp;
    milldata.y_pos{j} = yp;
    
    finalarray_x=[];
    finalarray_y=[];

	if milldata.consecutive == 0 % Make milling pattern random or maximally separated
        else
        
            ring_freq = length(milldata.x_pos{j});
            section =ring_freq/3;
            idx = ones(1,length(milldata.x_pos{j}));
            for i = 2:length(milldata.x_pos{j})
                if idx(i-1)+section > ring_freq
                    if i > ring_freq 
                        if i > 2*ring_freq
    %                         idx(i) = idx(i-ring)+ring;
    %                         else
                            idx(i) = idx(i-ring_freq)+ring_freq;
                         end
                        else
                        idx(i) = i-(2*((i-1)/3));
                    end
                    else
                    idx(i) = idx(i-1)+section;
                end
            	
            end
            milldata.x_pos{j} = milldata.x_pos{j}(idx);
            milldata.y_pos{j} = milldata.y_pos{j}(idx);
%         section = length(hdata.mill.x_dot)/3;
%         for i = 2:length(hdata.mill.x_dot)
%             if idx(i-1)+section > length(hdata.mill.x_dot)
%                 idx(i) = i-(2*((i-1)/3));
%             else
%                 idx(i) = idx(i-1)+section;
%             end
%         end

        end
    end
        for idx = 1:NoringsSpinner.Value
        finalarray_x = [finalarray_x, milldata.x_pos{idx}];
        finalarray_y = [finalarray_y, milldata.y_pos{idx}];
        milldata.x_dot = finalarray_x;
        milldata.y_dot = finalarray_y; 
    end    
    guidata(UIFigure,milldata);   
    end

    function plot_milling()
        milldata = guidata(UIFigure);
        totalshots = length(milldata.x_dot);
        f = linspace(1,10,totalshots);
        xc=0;
        yc=0;
        milldata.mill_plot = scatter(UIAxes, milldata.x_dot, milldata.y_dot',[],'+','LineWidth',2,'CData',f);
        guidata(UIFigure,milldata)
    end
    
    %% GUI Layout

        % Create FrequencySpinnerLabel
        FrequencySpinnerLabel = uilabel(UIFigure);
        FrequencySpinnerLabel.HorizontalAlignment = 'right';
        FrequencySpinnerLabel.Position = [18 165 62 22];
        FrequencySpinnerLabel.Text = 'Frequency';
        % Create FrequencySpinner
        FrequencySpinner = uispinner(UIFigure, 'ValueChangedFcn',@FrequencySpinnerValueChanged);
        FrequencySpinner.Position = [90 154 82 42];
        FrequencySpinner.Limits = [1 100];
        FrequencySpinner.Step = 3;
        FrequencySpinner.Value = milldata.freq(1);

        % Create RadiusSpinnerLabel
        RadiusSpinnerLabel = uilabel(UIFigure);
        RadiusSpinnerLabel.HorizontalAlignment = 'right';
        RadiusSpinnerLabel.Position =  [31 105 43 22];
        RadiusSpinnerLabel.Text = 'Radius (um)';
        % Create RadiusSpinner
        RadiusSpinner = uispinner(UIFigure, 'ValueChangedFcn', @RadiusSpinnerValueChanged);
        RadiusSpinner.Position = [90 95 82 42];
        RadiusSpinner.Limits = [1 100];
        RadiusSpinner.Value = milldata.rad(1);

        % Create RotationSpinnerLabel
        RotationSpinnerLabel = uilabel(UIFigure);
        RotationSpinnerLabel.HorizontalAlignment = 'right';
        RotationSpinnerLabel.Position = [24 43 50 22];
        RotationSpinnerLabel.Text = 'Rotation (deg)';
        % Create RotationSpinner
        RotationSpinner = uispinner(UIFigure, 'ValueChangedFcn', @RotationSpinnerValueChanged);
        RotationSpinner.Position = [91 33 82 42];
        RadiusSpinner.Limits = [0 100];
        RotationSpinner.Value = milldata.rot(1);

        % Create NoringsSpinnerLabel
        NoringsSpinnerLabel = uilabel(UIFigure);
        NoringsSpinnerLabel.HorizontalAlignment = 'right';
        NoringsSpinnerLabel.Position = [17 295 54 22];
        NoringsSpinnerLabel.Text = 'No. rings';
        % Create NoringsSpinner
        NoringsSpinner = uispinner(UIFigure, 'ValueChangedFcn', @NoringsSpinnerValueChanged);
        NoringsSpinner.Position = [89 285 82 42];
        RadiusSpinner.Limits = [1 100];
        NoringsSpinner.Value = milldata.n_rings;

        % Create ShootingOrderPanel
        ShootingOrderPanel = uipanel(UIFigure);
        ShootingOrderPanel.TitlePosition = 'centertop';
        ShootingOrderPanel.Title = 'Shooting Order';
        ShootingOrderPanel.Position = [591 216 129 107];

        % Create ConsecutiveCheckBox
        ConsecutiveCheckBox = uicheckbox(ShootingOrderPanel, 'ValueChangedFcn', @ConsecutiveCheckBoxValueChanged);
        ConsecutiveCheckBox.Text = 'Consecutive';
        ConsecutiveCheckBox.Position = [16 59 80 22];

        % Create RadialCheckBox
        RadialCheckBox = uicheckbox(ShootingOrderPanel, 'ValueChangedFcn', @RadialCheckBoxValueChanged);
        RadialCheckBox.Text = 'Radial';
        RadialCheckBox.Position = [16 35 87 22];

        % Create ClockwiseCheckBox
        ClockwiseCheckBox = uicheckbox(ShootingOrderPanel, 'ValueChangedFcn', @ClockwiseCheckBoxValueChanged);
        ClockwiseCheckBox.Text = 'Clockwise';
        ClockwiseCheckBox.Position = [16 9 87 22];

        % Create EccentricityPanel
        EccentricityPanel = uipanel(UIFigure);
        EccentricityPanel.TitlePosition = 'centertop';
        EccentricityPanel.Title = 'Eccentricity';
        EccentricityPanel.Position = [591 96 129 102];

        % Create RatioSpinnerLabel
        RatioSpinnerLabel = uilabel(EccentricityPanel);
        RatioSpinnerLabel.HorizontalAlignment = 'right';
        RatioSpinnerLabel.Position = [19 43 34 22];
        RatioSpinnerLabel.Text = 'Ratio';
        % Create RatioSpinner
        RatioSpinner = uispinner(EccentricityPanel, 'ValueChangedFcn', @RatioSpinnerValueChanged);
        RatioSpinner.Position = [57 43 57 22];
        RatioSpinner.Value = milldata.ratio;

        % Create AngleSpinnerLabel
        AngleSpinnerLabel = uilabel(EccentricityPanel);
        AngleSpinnerLabel.HorizontalAlignment = 'right';
        AngleSpinnerLabel.Position = [15 9 36 22];
        AngleSpinnerLabel.Text = 'Angle';
        % Create AngleSpinner
        AngleSpinner = uispinner(EccentricityPanel, 'ValueChangedFcn', @AngleSpinnerValueChanged);
        AngleSpinner.Position = [56 9 58 22];
        AngleSpinner.Value = milldata.ecc_angle;

        % Create RingDropDownLabel
        RingDropDownLabel = uilabel(UIFigure);
        RingDropDownLabel.HorizontalAlignment = 'right';
        RingDropDownLabel.Position = [50 230 30 22];
        RingDropDownLabel.Text = 'Ring';
        % Create RingDropDown
        RingDropDown = uidropdown(UIFigure, 'ValueChangedFcn', @RingDropDownValueChanged);
        RingDropDown.Items = {'1', '2', '3', '4', '5', '6'};
        RingDropDown.ItemsData = [1 2 3 4 5 6];
        RingDropDown.Position = [95 230 100 22];
        RingDropDown.Value = 1;

        % Create SaveButton
        SaveButton = uibutton(UIFigure, 'push', 'ButtonPushedFcn', @SaveButtonPushed);
        SaveButton.Position = [592 29 129 44];
        SaveButton.Text = 'Save';

        % Create UIAxes
        UIAxes = uiaxes(UIFigure);
        UIAxes.XLim = [-80 80];
        UIAxes.YLim = [-80 80];
        UIAxes.XAxisLocation = 'origin';
        UIAxes.XTick = [-60 -40 -20 0 20 40 60];
        UIAxes.XTickLabel = '';
        UIAxes.YAxisLocation = 'origin';
        UIAxes.YTick = [-60 -40 -20 0 20 40 60];
        UIAxes.YTickLabel = '';
        UIAxes.Position = [240 34 300 300];
        UIAxes.YTickLabel = '';
        plot_milling

        % Show the UIFigureure after all components are created
        UIFigure.Visible = 'on';

        % Construct UIFigure
        function UIFigure = Mill_pattern

            % Create UIUIFigureure and components
            createComponents(UIFigure)

            % Register the UIFigure with UIFigure Designer
            registerUIFigure(UIFigure, UIFigure)

            if nargout == 0
                clear UIFigure
            end
        end

        % Code that executes before UIFigure deletion
        function delete(UIFigure)

            % Delete UIUIFigureure when UIFigure is deleted
            delete(UIFigure)
        end
end