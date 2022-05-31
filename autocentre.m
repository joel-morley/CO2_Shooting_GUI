function autocentre(varargin)

    % Properties that correspond to app components
        UIFigure=matlab.ui.Figure;
        UpButton=matlab.ui.control.Button;
        LiveButton= matlab.ui.control.StateButton;
        Image =  matlab.ui.control.Image;
        DownButton = matlab.ui.control.Button;
        LeftButton =  matlab.ui.control.Button;
        RightButton =  matlab.ui.control.Button;
        OutButton= matlab.ui.control.Button;
        InButton = matlab.ui.control.Button;
        AutoCentreButton = matlab.ui.control.Button;
        CloseButton = matlab.ui.control.Button;
        StepsizemmButtonGroup = matlab.ui.container.ButtonGroup;
        ssBig = matlab.ui.control.ToggleButton;
        ssSmall = matlab.ui.control.ToggleButton;
        ssTiny= matlab.ui.control.ToggleButton;
        ssLarge = matlab.ui.control.ToggleButton;

        step_size = 1;
        C885 = 1;
        AtCube = 1;
        

        function [I] = grabImage(app)
            [err, M, w, h] = uEye_camera(1); % fetch image
            M = reshape(M+128, 3, w*h);
            I = reshape(M, 3, w, h);
            I = permute(I, [3 2 1])/255;
        end

        % Code that executes after component creation
        function startupFcn(app)
            uEye_camera(0); % initializing camera
            uEye_camera(8, 10); %exposure time

            % For PI
            clear C885
            % For attocube
            clear AtCube
            PI_Controller = PI_GCS_Controller ();
            PI_Controller.Destroy();
            app.C885 = PI_control();
            app.AtCube = AtCube_control();
            opt = processPropertyValuePairs({'loadsettings'}, varargin, 'struct', true);
            loadsettings = true;
            if isfield(opt, 'loadsettings') && ~opt.loadsettings
                loadsettings = false;  
            end
            
            app.AtCube.turnOn()
        end

        % Value changed function: LiveButton
        function LiveButtonValueChanged(app, event)
            value = app.LiveButton.Value;    
            while value == true
                    I = grabImage(app);
                    app.Image.ImageSource = I;
                    pause(0.1);
                end
        end

        % Callback function
        function CloseButtonValueChanged(app, event)
            uEye_camera(4);
            close all
        end

        % Button pushed function: UpButton
        function UpButtonPushed(app, event)
             y_pos = app.C885.getPosition_y();
             app.C885.move_y(y_pos-app.step_size);
        end

        % Selection changed function: StepsizemmButtonGroup
        function ButtonValueChanged(app, event)
            if app.ssLarge.value == true
                app.step_size = 5;
            end
            if app.ssBig.value == true
                app.step_size = 1;
            end
            if app.ssSmall.value == true
                app.step_size = 0.5;
            end
            if app.ssTiny.value == true
                app.step_size = 0.1;
            end
        end

        % Callback function
        function Button_2ValueChanged(app, event)
            value = app.Button_2.Value;
            app.step_size = 0.5;
        end

        % Callback function
        function Button_3ValueChanged(app, event)
            value = app.Button_3.Value;
            app.step_size = 1;
        end

        % Callback function
        function Button_4ValueChanged(app, event)
            value = app.Button_4.Value;
            app.step_size = 5;
        end

        % Button pushed function: CloseButton
        function CloseButtonPushed(app, event)
            uEye_camera(4);
            app.C885.delete()
            app.AtCube.delete()
            delete(app)
        end

        % Callback function
        function UIFigureCloseRequest(app, event)
            
            
        end

        % Button pushed function: DownButton
        function DownButtonPushed(app, event)
             y_pos = app.C885.getPosition_y();
             app.C885.move_y(y_pos+app.step_size);
        end

        % Button pushed function: RightButton
        function RightButtonPushed(app, event)
            x_pos = app.C885.getPosition_x();
            app.C885.move_x(x_pos+app.step_size);
        end

        % Button pushed function: InButton
        function InButtonPushed(app, event)
            z_pos = app.AtCube.getPosition_z();
            app.AtCube.move_z(z_pos+app.step_size);
        end

        % Button pushed function: OutButton
        function OutButtonPushed(app, event)
            z_pos = app.AtCube.getPosition_z();
            app.AtCube.move_z(z_pos-app.step_size);
        end

        % Button pushed function: LeftButton
        function LeftButtonPushed(app, event)
            x_pos = app.C885.getPosition_x();
            app.C885.move_x(x_pos-app.step_size);
        end
    end

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Color = [0.94 0.94 0.94];
            app.UIFigure.Position = [100 100 612 364];
            app.UIFigure.Name = 'MATLAB App';

            % Create UpButton
            app.UpButton = uibutton(app.UIFigure, 'push');
            app.UpButton.ButtonPushedFcn = @UpButtonPushed;
            app.UpButton.Position = [469 213 70 38];
            app.UpButton.Text = 'Up';

            % Create LiveButton
            app.LiveButton = uibutton(app.UIFigure, 'state');
            app.LiveButton.ValueChangedFcn = createCallbackFcn(app, @LiveButtonValueChanged, true);
            app.LiveButton.Text = 'Live';
            app.LiveButton.Position = [23 315 99 34];

            % Create Image
            app.Image = uiimage(app.UIFigure);
            app.Image.Position = [13 15 379 286];

            % Create DownButton
            app.DownButton = uibutton(app.UIFigure, 'push');
            app.DownButton.ButtonPushedFcn = createCallbackFcn(app, @DownButtonPushed, true);
            app.DownButton.Position = [469 103 70 38];
            app.DownButton.Text = 'Down';

            % Create LeftButton
            app.LeftButton = uibutton(app.UIFigure, 'push');
            app.LeftButton.ButtonPushedFcn = createCallbackFcn(app, @LeftButtonPushed, true);
            app.LeftButton.Position = [413 158 70 38];
            app.LeftButton.Text = 'Left';

            % Create RightButton
            app.RightButton = uibutton(app.UIFigure, 'push');
            app.RightButton.ButtonPushedFcn = createCallbackFcn(app, @RightButtonPushed, true);
            app.RightButton.Position = [529 158 70 38];
            app.RightButton.Text = 'Right';

            % Create OutButton
            app.OutButton = uibutton(app.UIFigure, 'push');
            app.OutButton.ButtonPushedFcn = createCallbackFcn(app, @OutButtonPushed, true);
            app.OutButton.Position = [424 35 70 38];
            app.OutButton.Text = 'Out';

            % Create InButton
            app.InButton = uibutton(app.UIFigure, 'push');
            app.InButton.ButtonPushedFcn = createCallbackFcn(app, @InButtonPushed, true);
            app.InButton.Position = [509 35 70 38];
            app.InButton.Text = 'In';

            % Create AutoCentreButton
            app.AutoCentreButton = uibutton(app.UIFigure, 'push');
            app.AutoCentreButton.FontWeight = 'bold';
            app.AutoCentreButton.FontColor = [0.0745 0.6235 1];
            app.AutoCentreButton.Position = [150 315 99 34];
            app.AutoCentreButton.Text = {'Auto-Centre'; ''};

            % Create CloseButton
            app.CloseButton = uibutton(app.UIFigure, 'push');
            app.CloseButton.ButtonPushedFcn = createCallbackFcn(app, @CloseButtonPushed, true);
            app.CloseButton.Position = [276 315 99 34];
            app.CloseButton.Text = 'Close';

            % Create StepsizemmButtonGroup
            app.StepsizemmButtonGroup = uibuttongroup(app.UIFigure);
            app.StepsizemmButtonGroup.SelectionChangedFcn = createCallbackFcn(app, @ButtonValueChanged, true);
            app.StepsizemmButtonGroup.TitlePosition = 'centertop';
            app.StepsizemmButtonGroup.Title = 'Step size (mm)';
            app.StepsizemmButtonGroup.Position = [391 265 213 64];

            % Create ssBig
            app.ssBig = uitogglebutton(app.StepsizemmButtonGroup);
            app.ssBig.Text = '1';
            app.ssBig.Position = [113 11 33 22];
            app.ssBig.Value = true;

            % Create ssSmall
            app.ssSmall = uitogglebutton(app.StepsizemmButtonGroup);
            app.ssSmall.Text = '0.5';
            app.ssSmall.Position = [63 11 33 22];

            % Create ssTiny
            app.ssTiny = uitogglebutton(app.StepsizemmButtonGroup);
            app.ssTiny.Text = '0.1';
            app.ssTiny.Position = [13 11 33 22];

            % Create ssLarge
            app.ssLarge = uitogglebutton(app.StepsizemmButtonGroup);
            app.ssLarge.Text = '5';
            app.ssLarge.Position = [163 11 33 22];

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end



        % Construct app
        function app = app1

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.UIFigure)

            % Execute the startup function
            runStartupFcn(app, @startupFcn)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
        end

