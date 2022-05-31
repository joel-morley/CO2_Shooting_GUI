function simpleApp2
%% Startup Processes

uEye_camera(0); % initializing camera
uEye_camera(8, 10); %exposure time
pause(0.5)

% For PI
% addpath ( 'C:\Users\Public\PI\PI_MATLAB_Driver_GCS2' );
clear C885
% For attocube
clear AtCube

PI_Controller = PI_GCS_Controller ();
PI_Controller.Destroy();

% handles to X, Y, Z motor axes
if ( ~exist ( 'C885', 'var' ) || ~isa ( C885, 'PI_control' ) )
    C885 = PI_control();
else
    error('There exis the var. "C885"')
end

% C885.turnOnXY()
if ( ~exist ( 'AtCube', 'var' ) || ~isa ( AtCube, 'AtCube_control' ) )
    AtCube = AtCube_control();
else
    error('There exis the var. "AtCube"')
end
AtCube.turnOn()

%beamsplitter position
C885.move_z(23.9);

fig = uifigure;
fig.Color = [0.94 0.94 0.94];
fig.Position = [100 400 773 400];
fig.Name = 'Auto-Centre';
Im = grabImage();
hdata = guidata(fig);
hdata.step_size = 1;
hdata.Image = uint8(Im);
guidata(fig,hdata);

%% Callbacks

        % Button pushed function: SaveButton
        function SaveImagecallback(src,event)
        	persistent lastdir;
        	while SaveButton.Value == true
                if lastdir
                	[filename, pathname] = uiputfile('*.tif', 'Save Image', lastdir);
                else
                    [filename, pathname] = uiputfile('*.tif', 'Save Image');
                end
                if filename
                    save(fullfile(pathname, filename), Image)
                    lastdir = pathname;
                end
        	end
        end

        % Button pushed function: CloseButton
        function CloseButtonPushed(fig, event)
            stop(tm);
            uEye_camera(4);
            C885.delete()
            AtCube.delete()
            closereq(); 
        end
    
        % Button pushed function: AutoCentreButton
        function AutoCentreButtonPushed(fig, event)
            Im_0 = hdata.Image;
            IR_filt=imgaussfilt(Im_0 ,5,'FilterSize',21); %% Filter
            [centers, radii] = imfindcircles(IR_filt, [180 200]);

            if isempty(centers)
                disp('No fiber found')
            else
                if size(centers(1,:))>1
                    disp('Centre may not be correct, try autofocus first')
                end
            x_start = C885.getPosition_x();
            y_start = C885.getPosition_y();

            pixel=0.000329; % Pixel FOV in mm
            diff_x = pixel*(968-centers(1,1));
            diff_y = pixel*(548-centers(1,2));
            C885.move_x(x_start-diff_x)
            C885.move_y(y_start-diff_y)
            end
        end

        % Button pushed function: AutoFocusButton
        function AutoFocusButtonPushed(fig, event)
            step_n = 0.014;
            range_n = 0.210 ;
            focus = AtCube.getPosition_z();
            start0 = focus-(range_n/2);
            zdata = linspace(start0,start0+range_n,range_n/step_n);
            AtCube.move_z(start0);
            pause(0.2);

            FM = focus_scan(3,zdata);

            % beam waist, x0, fitting
            [fitpar_x, zfit_x, vfit_x,z_highres, x_highres] = fplane_wide_Fit(zdata, FM);
            % plot data with fitted curve and print w0, z0
            figure        
            plot(zdata, FM, 'o-')
            hold on
            plot(z_highres, x_highres);
            hold off

            step_n = 0.0014;
            range_n = 0.021;
            % [y_min,min_index] = min(FM);
            % y_min_pos = zdata(min_index);
            focus = fitpar_x(2);
            start0 = focus-(range_n/2);
            zdata = linspace(start0,start0+range_n,range_n/step_n);
            AtCube.move_z(start0);
            pause(0.2);
    
            FM = focus_scan(6,zdata);

            % beam waist, x0, fitting
            [fitpar_x, zfit_x, vfit_x,z_highres, x_highres] = fplane_wide_Fit(zdata, FM);
            % plot data with fitted curve and print w0, z0
            figure        
            plot(zdata, FM, 'o-')
            hold on
            plot(z_highres, x_highres);
            hold off
            
        end
    
        % Button pushed function: ReconstructSurfaceButton
        function ReconstructSurfaceButtonPushed(fig, event)
        end
    
        % Value changed function: GridSwitch
        function GridSwitchValueChanged(app, event)
            value = GridSwitch.Value;
            disp(value)
            if value == false
                set(gridax,'visible', 'on');
            else
                set(gridax,'visible', 'off');
            end
        end
    
        % Selection changed function: StepsizemmButtonGroup
        function ButtonValueChanged(fig, event)
                    if ssLarge.Value == true
                        hdata = guidata(fig);
                        hdata.step_size = 5;
                        disp(hdata.step_size)
                        guidata(fig,hdata)
                    end
                    if ssBig.Value == true
                        hdata = guidata(fig);
                        hdata.step_size = 1;
                        guidata(fig,hdata)
                    end
                    if ssSmall.Value == true
                        hdata = guidata(fig);
                        hdata.step_size = 0.5;
                        guidata(fig,hdata)
                    end
                    if ssTiny.Value == true
                        hdata = guidata(fig);
                        hdata.step_size = 0.1;
                        guidata(fig,hdata)
                    end
                end
    
        % Button pushed function: UpButton
        function UpButtonPushed(fig, event)
             y_pos = C885.getPosition_y();
             disp(y_pos);
             C885.move_y(y_pos+hdata.step_size);
        end

        % Button pushed function: DownButton
        function DownButtonPushed(fig, event)
             y_pos = C885.getPosition_y();
             C885.move_y(y_pos-hdata.step_size);
        end

        % Button pushed function: RightButton
        function RightButtonPushed(fig, event)
            x_pos = C885.getPosition_x();
            C885.move_x(x_pos-hdata.step_size);
        end
    
        % Button pushed function: LeftButton
        function LeftButtonPushed(fig, event)
            x_pos = C885.getPosition_x();
            C885.move_x(x_pos+hdata.step_size);
        end

        % Button pushed function: InButton
        function InButtonPushed(fig, event)
            z_pos = AtCube.getPosition_z();
            AtCube.move_z(z_pos-hdata.step_size);
        end

        % Button pushed function: OutButton
        function OutButtonPushed(fig, event)
            z_pos = AtCube.getPosition_z();
            AtCube.move_z(z_pos+hdata.step_size);
        end

%% Helper functions
function ImBG_autoCont = grabImage()
    [err, M, w, h] = uEye_camera(1); % fetch image
        if err>0
            error('Capturing an image failed!')
        end
    M(M<0)=M(M<0)+256;
    M = reshape(M, 3, w*h);
    I = reshape(M, 3, w, h);
    I = permute(I, [3 2 1]);
    Im_BG = 0.7*I(:,:,2)+0.3*I(:,:,3);
    ImBG_autoCont = rescale(Im_BG,0,256);
end

function updateImage(src,evt)
    hdata = guidata(fig);
    Im_opt_cont = grabImage();
    hdata.Image = uint8(Im_opt_cont);
%     I = grabImage();
%     fig.UserData.Image = grabImage();
%     fig.UserData.test = fig.UserData.test+1;
    imshow(hdata.Image, 'Parent', Imageax);
    guidata(fig,hdata);
end

function [FM] = focus_scan(frames,zdata)
    for i = 1:size(zdata,2)
      for j = 1:frames
        ImBG_autoCont(:,:,j) = grabImage();
      end
        Imcont3 = mean(ImBG_autoCont,3);
        IR_filt=imgaussfilt(Imcont3,4,'FilterSize',11); %% Filter
        AtCube.move_z(zdata(i));
        FM(i) = fmeasure(Imcont3,'CONT',[545 968 500 500]);
    end
end

%% GUI layout
            % Create Image, must use plot axes for the image in order to also draw the
            % alignment lines
            Imageax = uiaxes(fig);
            gridax = uiaxes(fig);
            Imageax.Position = [23 63 479 277];
            Imageax.XLim = [0 1936]; % Set limits of axes
            Imageax.YLim = [0 1096];
            gridax.Position = [23 63 479 277];
            gridax.XLim = [0 1936]; % Set limits of axes
            gridax.YLim = [0 1096];
            gridax.Color = 'None';
            Imageax.TickLength = [0 0];
            Imageax.XTick = [968];
            Imageax.XTickLabel = '';
            Imageax.YTick = [548];
            Imageax.YTickLabel = '';
            gridax.TickLength = [0 0];
            gridax.XTick = [968];
            gridax.XTickLabel = '';
            gridax.YTick = [548];
            gridax.YTickLabel = '';
            gridax.XGrid = 'on';
            gridax.YGrid = 'on';
            gridax.GridAlpha = 1;
            gridax.GridLineStyle = '--';
            drawcircle(gridax,'Center',[968,548],'Radius',250,'Color','k','FaceAlpha',0.01,'LineWidth',0.1,'MarkerSize',0.1);
            tm = timer('TimerFcn', @updateImage,'Period', 0.1,'ExecutionMode','fixedSpacing','BusyMode','drop');
            start(tm);

               
            % Create SaveButton
            SaveButton = uibutton(fig,'state','ValueChangedFcn',@SaveImagecallback);
            SaveButton.Text = 'Save';
            SaveButton.Position = [34 349 99 34];
            
            % Create AutoCentreButton
            AutoCentreButton = uibutton(fig, 'push','ButtonPushedFcn', @AutoCentreButtonPushed);
            AutoCentreButton.FontWeight = 'bold';
            AutoCentreButton.FontColor = [0.0745 0.6235 1];
            AutoCentreButton.Position = [160 349 99 34];
            AutoCentreButton.Text = {'Auto-Centre'; ''};
            
            % Create AutoFocusButton
            AutoFocusButton = uibutton(fig, 'push','ButtonPushedFcn', @AutoFocusButtonPushed);
            AutoFocusButton.FontWeight = 'bold';
            AutoFocusButton.FontColor = [0.0745 0.6235 1];
            AutoFocusButton.Position = [286 349 99 34];
            AutoFocusButton.Text = {'Auto-Focus'; ''};
            
            % Create ReconstructSurfaceButton
            ReconstructSurfaceButton = uibutton(fig,'ButtonPushedFcn', @ReconstructSurfaceButtonPushed);
            ReconstructSurfaceButton.Position = [47 26 447 33];
            ReconstructSurfaceButton.Text = 'Reconstruct Surface';
            
            % Create GridSwitchLabel
            GridSwitchLabel = uilabel(fig);
            GridSwitchLabel.HorizontalAlignment = 'center';
            GridSwitchLabel.Position = [539 239 28 22];
            GridSwitchLabel.Text = 'Grid';

            % Create GridSwitch
            GridSwitch = uiswitch(fig, 'slider','ValueChangedFcn',@GridSwitchValueChanged);
            GridSwitch.Position = [531 267 45 20];

            % Create CloseButton
            CloseButton = uibutton(fig, 'push','ButtonPushedFcn', @CloseButtonPushed);
            CloseButton.Position = [411 349 99 34];
            CloseButton.Text = 'Close';
            
            % Create UpButton
            UpButton = uibutton(fig, 'push','ButtonPushedFcn', @UpButtonPushed);
            UpButton.Position = [617 249 70 38];
            UpButton.Text = 'Up';

            % Create DownButton
            DownButton = uibutton(fig, 'push','ButtonPushedFcn', @DownButtonPushed);
            DownButton.Position = [617 139 70 38];
            DownButton.Text = 'Down';
            
            % Create LeftButton
            LeftButton = uibutton(fig, 'push','ButtonPushedFcn', @LeftButtonPushed);
            LeftButton.Position = [561 194 70 38];
            LeftButton.Text = 'Left';
            
            % Create RightButton
            RightButton = uibutton(fig, 'push','ButtonPushedFcn', @RightButtonPushed);
            RightButton.Position = [677 194 70 38];
            RightButton.Text = 'Right';

            % Create OutButton
            OutButton = uibutton(fig, 'push','ButtonPushedFcn', @OutButtonPushed);
            OutButton.Position = [572 71 70 38];
            OutButton.Text = 'Out';

            % Create InButton
            InButton = uibutton(fig, 'push','ButtonPushedFcn', @InButtonPushed);
            InButton.Position = [657 71 70 38];
            InButton.Text = 'In';
            
            % Create StepsizemmButtonGroup
            StepsizemmButtonGroup = uibuttongroup(fig,'SelectionChangedFcn',@ButtonValueChanged);
            StepsizemmButtonGroup.TitlePosition = 'centertop';
            StepsizemmButtonGroup.Title = 'Step size (mm)';
            StepsizemmButtonGroup.Position = [539 301 213 64];

            % Create ssTiny
            ssTiny = uitogglebutton(StepsizemmButtonGroup);
            ssTiny.Text = '0.1';
            ssTiny.Position = [13 11 33 22];
            
            % Create ssSmall
            ssSmall = uitogglebutton(StepsizemmButtonGroup);
            ssSmall.Text = '0.5';
            ssSmall.Position = [63 11 33 22];
             
            % Create ssBig
            ssBig = uitogglebutton(StepsizemmButtonGroup);
            ssBig.Text = '1';
            ssBig.Position = [113 11 33 22];
            ssBig.Value = true;

            % Create ssLarge
            ssLarge = uitogglebutton(StepsizemmButtonGroup);
            ssLarge.Text = '5';
            ssLarge.Position = [163 11 33 22];
            
end

