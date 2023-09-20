classdef LiquidIngressSimulator_exported < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                        matlab.ui.Figure
        TabGroup                        matlab.ui.container.TabGroup
        MicrostructuresettingTab        matlab.ui.container.Tab
        PoresettingPanel                matlab.ui.container.Panel
        PoreisotropicalityEditField     matlab.ui.control.NumericEditField
        Poreisotropicity0110Label       matlab.ui.control.Label
        UpdateButton                    matlab.ui.control.Button
        PoresizedistributionDropDown    matlab.ui.control.DropDown
        PoresizedistributionDropDownLabel  matlab.ui.control.Label
        StandarddeviationumEditField    matlab.ui.control.NumericEditField
        StandarddeviationumEditFieldLabel  matlab.ui.control.Label
        MeanradiusumEditField           matlab.ui.control.NumericEditField
        MeanradiusumEditFieldLabel      matlab.ui.control.Label
        Panel                           matlab.ui.container.Panel
        ColormapDropDown                matlab.ui.control.DropDown
        ColormapDropDownLabel           matlab.ui.control.Label
        ControlboxPanel                 matlab.ui.container.Panel
        LiquidfrontthreadholdEditField  matlab.ui.control.NumericEditField
        ThreadholdLabel                 matlab.ui.control.Label
        TimemsEditField                 matlab.ui.control.NumericEditField
        TimemsEditFieldLabel            matlab.ui.control.Label
        ResetButton                     matlab.ui.control.Button
        SaveButton                      matlab.ui.control.Button
        StopButton_2                    matlab.ui.control.Button
        StartButton_2                   matlab.ui.control.Button
        LoadstructureButton             matlab.ui.control.Button
        UIAxes23                        matlab.ui.control.UIAxes
        UIAxes22                        matlab.ui.control.UIAxes
        UIAxes21                        matlab.ui.control.UIAxes
        AssingdatainworkspaceButton     matlab.ui.control.Button
        MicrostructuredetailsPanel      matlab.ui.container.Panel
        StructureheightmmEditField      matlab.ui.control.NumericEditField
        StructureheightmmEditFieldLabel  matlab.ui.control.Label
        CapillaryheightumEditField      matlab.ui.control.NumericEditField
        CapillaryheightumEditFieldLabel  matlab.ui.control.Label
        CapillaryDistributionPanel      matlab.ui.container.Panel
        UIAxes11                        matlab.ui.control.UIAxes
        UIAxes14                        matlab.ui.control.UIAxes
        UIAxes13                        matlab.ui.control.UIAxes
        UIAxes12                        matlab.ui.control.UIAxes
        MicrostructuresettingPanel      matlab.ui.container.Panel
        InterlayerconnectivityEditField  matlab.ui.control.NumericEditField
        Interlayerconnectivity50500Label  matlab.ui.control.Label
        LayernumberEditField            matlab.ui.control.NumericEditField
        Layernumber105000Label          matlab.ui.control.Label
        CapillarysizenumberEditField    matlab.ui.control.NumericEditField
        CapillarysizenumberEditFieldLabel  matlab.ui.control.Label
        CapillaryarraynumberEditField   matlab.ui.control.NumericEditField
        CapillaryarraynumberEditFieldLabel  matlab.ui.control.Label
        DisintegrationmeidumsettingPanel  matlab.ui.container.Panel
        ContactangledegEditField        matlab.ui.control.NumericEditField
        ContactangledegEditFieldLabel   matlab.ui.control.Label
        SurfaceTensionmNmEditField      matlab.ui.control.NumericEditField
        SurfaceTensionmNmEditFieldLabel  matlab.ui.control.Label
        WatertemperatureCEditField      matlab.ui.control.NumericEditField
        WatertemperatureCLabel          matlab.ui.control.Label
        DynamicViscositymPasEditField   matlab.ui.control.NumericEditField
        DynamicViscosityPasLabel        matlab.ui.control.Label
        demoTab                         matlab.ui.container.Tab
        StopButton                      matlab.ui.control.Button
        StartButton                     matlab.ui.control.Button
        UIAxes                          matlab.ui.control.UIAxes
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Liquid Ingress Simulation
% Created by Jongmin Lee, University of Cambridge, UK. 2023
% Contact Information: jl2112@cam.ac.uk
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    properties (Access = private)
        RandTimer          % Timer object
        PlotLine           % Line object
        Structure          % Porous medium structure
        Hydration          % Hydration status
        WetMatMap1         % Wet Matrix Map 1, ecch cell
        WetMatMap2         % Wet Matrix Map 2, layer
        WetMatMap3         % Liquid ingress plot
        SimTimer           % Simulation timer
        SimTime            % Simulation time
    end

    methods (Access = private)
    
        function RandTimerFcn(app,~,~)
            % Generate a random number
            randnum = rand;
    
            % Update YData in plot
            ydata = app.PlotLine.YData;
            ydata = circshift(ydata,1);
            ydata(1) = randnum;
            app.PlotLine.YData = ydata;
       end        
       
       function SimTimerFcn(app,~,~)
            wetMat = app.Hydration.wetMat;
            wetVec = app.Hydration.wetVec;
            liquidFrontTime = app.Hydration.liquidFrontTime;
            simTime = app.SimTime;

            app.updateWetMat;
    
            % Update CData in plot
            app.WetMatMap1.CData = wetMat;
            app.WetMatMap2.YData = wetVec;
            app.WetMatMap3.XData = liquidFrontTime;

            simTime = simTime + 0.1;
            app.SimTime = simTime;

            app.TimemsEditField.Value = simTime;
        end
        
        function liquidUpdate(app)
            % Initial dissolution medium setting
            dynViscosityArray=[1.6735 1.619 1.5673 1.5182 1.4715 1.4271 1.3847 1.3444 1.3059 1.2692 1.234 1.2005 1.1683 1.1375 1.1081 1.0798 1.0526 1.0266 ... 
                1.0016 0.9775 0.9544 0.9321 0.9107 0.89 0.8701 0.8509 0.8324 0.8145 0.7972 0.7805 0.7644 0.7488 0.7337 0.7191 0.705 0.6913 0.678 0.6652 0.6527];
            surfaceTensionArray=[75.3671 75.2259 75.0841 74.9417 74.7688 74.6552 74.5111 74.3663 74.221 74.0752 73.9287 73.7816 73.634 73.4858 73.337 73.1877 ...
                73.0377 72.8872 72.7361 72.5845 72.4323 72.2795 72.1261 71.9722 71.8177 71.6627 71.5071 71.3509 71.1942 71.0369 70.879 70.7206 70.5616 70.4021 70.2421 ...
                70.0815 69.9203 69.7586 69.5963];
            % Data ranged from temp. 2 to 40
            % Reference: IAPWS

            liquidTemp = app.WatertemperatureCEditField.Value;
            dynViscosity = dynViscosityArray(liquidTemp-1); % in mPas
            app.DynamicViscositymPasEditField.Value = dynViscosity;
            surTension = surfaceTensionArray(liquidTemp-1); % in mN/m
            app.SurfaceTensionmNmEditField.Value = surTension;
            
        end
        
        function createConnectionMatrix(app)
            capArrayNum = app.CapillaryarraynumberEditField.Value;
            layerNum = app.LayernumberEditField.Value;
            interCon = app.InterlayerconnectivityEditField.Value/100;

            conMat = cell(layerNum,capArrayNum); % upstream connection
            conMatDown = cell(layerNum,capArrayNum); % downstream connection

            for i = 1:layerNum-1
                for j = 1:capArrayNum
                    conVec = single(ones(1,capArrayNum));
                    conVec_P = rand(1,capArrayNum);
                    conVec(conVec_P > interCon/capArrayNum) = 0;
                    conVec = conVec.*(1:capArrayNum);
                    conMat(i,j) = {nonzeros(conVec)};
                end
            end

            app.Structure.conMat = conMat;

            for i = 1: layerNum-1
                for j = 1:capArrayNum
                    if ~isempty(conMat{i,j})
                        for addIdx = conMat{i,j}'
                            addVec = conMatDown{i+1,addIdx};
                            addVec = [addVec; j];
                            conMatDown{i+1,addIdx} = addVec;
                        end
                    end
                end
            end

            app.Structure.conMatDown = conMatDown;

            statMat = zeros(layerNum,capArrayNum);
            wetMat = zeros(layerNum,capArrayNum);
            statMat(1,:) = 1;
            
            app.Hydration.statMat = single(statMat);
            app.Hydration.wetMat = wetMat;
        end
        
        function updateStatMat(app)
            statMat = app.Hydration.statMat;
        end
        
        function updateWetMat(app)
            wetMat = app.Hydration.wetMat;
            statMat = app.Hydration.statMat;
            lfThreashold = app.LiquidfrontthreadholdEditField.Value/100;
            liquidFrontTime = app.Hydration.liquidFrontTime;

            conMat = app.Structure.conMat;
            conMatDown = app.Structure.conMatDown;
            incVec = app.Structure.increment; % 0.01ms increment
            

            %update function
            for i = 1:10
                updateMat = zeros(size(wetMat));
                incMat = zeros(size(wetMat))+incVec;

                incMat(statMat==0) = 0;
                incMat(statMat==4) = 0;
                addMat = incMat;
                incMat(statMat==1) = 0;
                incMat(statMat==2) = 0;
                addMat = addMat + incMat;
                wetMat = wetMat + addMat;

                updateMat(wetMat>1) = 1;
                [y,x] = find(updateMat);
                
                if ~isempty(y)
                    for idx_row = y
                        for idx_column = x
                            if isequal(statMat(idx_row,idx_column),1)
                                    if ~isempty(conMat{idx_row,idx_column})
                                        for openIdx = conMat{idx_row,idx_column}'
                                            if isequal(statMat(idx_row+1,openIdx),0)
                                                   statMat(idx_row+1,openIdx) = 1;
                                            elseif isequal(statMat(idx_row+1,openIdx),2)
                                                   statMat(idx_row+1,openIdx) = 3;
                                            end

                                            if ~isempty(conMatDown{idx_row+1,openIdx})
                                                for downIdx = conMatDown{idx_row+1,openIdx}'
                                                    if isequal(statMat(idx_row,downIdx),0)
                                                        statMat(idx_row,downIdx) = 2;
                                                    elseif isequal(statMat(idx_row,downIdx),1)
                                                            statMat(idx_row,downIdx) = 3;
                                                    end
                                                end
                                            end
                                        end
                                    end

                            elseif isequal(statMat(idx_row,idx_column),2)
                                    if ~isempty(conMatDown{idx_row,idx_column})
                                        for openIdx = conMatDown{idx_row,idx_column}'
                                            if isequal(statMat(idx_row-1,openIdx),0)
                                                statMat(idx_row-1,openIdx) = 2;
                                            elseif isequal(statMat(idx_row-1,openIdx),1)
                                                    statMat(idx_row-1,openIdx) = 3;
                                            end

                                            if ~isempty(conMat{idx_row-1,openIdx})
                                                for upIdx = conMat{idx_row-1,openIdx}'
                                                    if isequal(statMat(idx_row,upIdx),0)
                                                        statMat(idx_row,upIdx) = 1;
                                                    elseif isequal(statMat(idx_row,upIdx),2)
                                                            statMat(idx_row,upIdx) = 3;
                                                    end
                                                end
                                            end
                                        end
                                    end
                            end
                        end
                    end
                end

                wetMat(wetMat>1) = 1;
                statMat(wetMat==1) = 4;
            end
            
            volVec = app.Structure.volume;
            totalVol = sum(volVec);
            wetVec = sum(wetMat.*volVec,2)/totalVol;
            
            updateLFTime = ones(size(liquidFrontTime));
            updateLFTime(wetVec < lfThreashold) = 0;
            updateLFTime(liquidFrontTime ~= 0 ) = 0;
            liquidFrontTime(updateLFTime == 1) = app.SimTime;
            
            app.Hydration.wetVec = wetVec;
            app.Hydration.wetMat = wetMat;
            app.Hydration.statMat = statMat;
            app.Hydration.liquidFrontTime = liquidFrontTime;
        end
        
        function resetSim(app)
            app.SimTime = 0;
            app.TimemsEditField.Value = app.SimTime;
            
            wetMat = app.Hydration.wetMat;
            wetMat(:,:) = 0;
            statMat = zeros(size(wetMat));
            statMat(1,:) = 1;
            
            app.Hydration.statMat = single(statMat);
            app.Hydration.wetMat = wetMat;
            LoadstructureButtonPushed(app);
            
        end
    end


    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app)
            % Configure x- and y- axis
            app.UIAxes.XLim = [0 60];
            app.UIAxes.XDir = 'reverse';
            app.UIAxes.YLim = [0 1];
            app.SimTime = 0;
            
            % Initial dissolution medium update
            app.liquidUpdate;

            % Initial plot is all zeros
            app.PlotLine = plot(app.UIAxes,0:60,zeros(1,61));
            
            % Create timer object
            app.RandTimer = timer(...
                'ExecutionMode', 'fixedRate', ...    % Run timer repeatedly
                'Period', 0.1, ...                     % Period is 0.1 second
                'BusyMode', 'queue',...              % Queue timer callbacks when busy
                'TimerFcn', @app.RandTimerFcn);      % Specify callback function

            app.SimTimer = timer(...
                'ExecutionMode', 'fixedRate', ...    % Run timer repeatedly
                'Period', 0.1, ...                     % Period is 0.1 second
                'BusyMode', 'queue',...              % Queue timer callbacks when busy
                'TimerFcn', @app.SimTimerFcn);      % Specify callback function
        end

        % Button pushed function: StartButton
        function StartButtonPushed(app, event)
            % If timer is not running, start it
            if strcmp(app.RandTimer.Running, 'off')
               start(app.RandTimer);
            end
        end

        % Button pushed function: StopButton
        function StopButtonPushed(app, event)
            % Stop the timer
            stop(app.RandTimer);
        end

        % Close request function: UIFigure
        function UIFigureCloseRequest(app, event)
            % Stop timer, then delete timer and app
            stop(app.RandTimer);
            delete(app.RandTimer);
            stop(app.SimTimer);
            delete(app.SimTimer);
            delete(app);
        end

        % Value changed function: WatertemperatureCEditField
        function WatertemperatureCEditFieldValueChanged(app, event)
            value = app.WatertemperatureCEditField.Value;
            app.liquidUpdate;
        end

        % Button pushed function: UpdateButton
        function UpdateButtonPushed(app, event)
            Rmean = app.MeanradiusumEditField.Value;
            Rsigma = app.StandarddeviationumEditField.Value;
            capArrayNum = app.CapillaryarraynumberEditField.Value;
            layerNum = app.LayernumberEditField.Value;
            interCon = app.InterlayerconnectivityEditField.Value;

            x=linspace(Rmean-4*Rsigma,Rmean+4*Rsigma,capArrayNum);
            Rpdf = normpdf(x,Rmean,Rsigma); % capillary radius distribution function
            dynVis = app.DynamicViscositymPasEditField.Value*10^-3;
            surTen = app.SurfaceTensionmNmEditField.Value*10^-3;
            conAng = deg2rad(app.ContactangledegEditField.Value);
            
            % Interpolate the radius function for the given array
            xq = linspace(0,max(cumtrapz(Rpdf)),capArrayNum+1);
            xq = (xq(1:end-1) + xq(2:end))/2;
            capRadius = interp1(cumtrapz(Rpdf),x,xq);
            capHeight = mean(capRadius)*2*app.PoreisotropicalityEditField.Value;
            matrixHeight = capHeight*app.LayernumberEditField.Value/1000;

            app.CapillaryheightumEditField.Value = capHeight;
            app.StructureheightmmEditField.Value = matrixHeight;

            capVolume = 2*pi*capRadius*capHeight;
            capTimeCoef = 2*dynVis*(capHeight*10^-6).^2/(surTen*10^-3*cos(conAng));
            capTime = capTimeCoef./(capRadius*10^-6);
            %capTime = 2*dynVis*(capHeight*10^-6)^2./(surTen*10^-3*capRadius*10^-6);


            app.Structure.radius = capRadius;
            app.Structure.volume = capVolume;
            app.Structure.time = capTime;
            app.Structure.increment = 0.01e-3./capTime;

            app.createConnectionMatrix;

            % Plot porous medium properties
            ux1 = app.UIAxes11;
            ux2 = app.UIAxes12;
            ux3 = app.UIAxes13;
            ux4 = app.UIAxes14;

            plot(ux1,x,Rpdf)
            plot(ux2,(1:size(xq,2)),capRadius);
            plot(ux3,(1:size(xq,2)),capVolume);
            plot(ux4,(1:size(xq,2)),capTime*1000);


        end

        % Button pushed function: AssingdatainworkspaceButton
        function AssingdatainworkspaceButtonPushed(app, event)
            assignin('base',"poreVectors",app.Structure);
            assignin('base',"hydrationMatrix",app.Hydration);
        end

        % Button pushed function: LoadstructureButton
        function LoadstructureButtonPushed(app, event)

            try
                statMat = app.Hydration.statMat;
                structureHeight = app.StructureheightmmEditField.Value;
            catch ME
                fig = app.UIFigure;
                uialert(fig,'Update the microstructure settings first!','Warning');
                return;
            end

            wetMat = app.Hydration.wetMat;
            volVec = app.Structure.volume;
            maxTime = max(app.Structure.time);
            layerNumber = app.LayernumberEditField.Value;
            totalTime = maxTime*layerNumber*1000;
            totalHeight = app.StructureheightmmEditField.Value;
            lfThreshold = app.LiquidfrontthreadholdEditField.Value/100;
            liquidFrontDepth = linspace(0,totalHeight,layerNumber);
            liquidFrontTime = zeros(layerNumber,1);
            totalVol = sum(volVec);

            ax1 = app.UIAxes21;
            ax2 = app.UIAxes22;
            ax3 = app.UIAxes23;
            
            clo(ax1)
            axis(ax1,"tight")
            cmap = app.ColormapDropDown.Value;
            clim(ax1,[0,1])
            colormap(ax1,cmap)
            colorbar(ax1)
            
            app.WetMatMap1 = imagesc(ax1,wetMat);
            ax1.YDir = 'normal';

            wetVec = sum(wetMat.*volVec,2)/totalVol;
            app.WetMatMap2 = barh(ax2,wetVec);
            xline(ax2,lfThreshold,'--');
            axis(ax2,"tight")
            ax2.YDir = 'normal';
            ax2.XLim = [0 1];

            app.WetMatMap3 = plot(ax3,liquidFrontTime,liquidFrontDepth,'.');
            axis(ax3,"tight")
            grid(ax3,"on");
            ax3.XLim = [0 totalTime];
            ax3.YLim = [0 totalHeight];

            app.Hydration.liquidFrontTime = liquidFrontTime;
            app.Hydration.wetVec = wetVec;
        end

        % Button pushed function: StartButton_2
        function StartButton_2Pushed(app, event)
            % If timer is not running, start it
            if strcmp(app.SimTimer.Running, 'off')
               start(app.SimTimer);
            end
        end

        % Button pushed function: StopButton_2
        function StopButton_2Pushed(app, event)
            % Stop the timer
            stop(app.SimTimer);
        end

        % Button pushed function: ResetButton
        function ResetButtonPushed(app, event)
            % If timer is not running, reset it
            if strcmp(app.SimTimer.Running, 'off')
               app.resetSim;
            end
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Position = [100 100 1651 847];
            app.UIFigure.Name = 'Random Number Generator';
            app.UIFigure.Resize = 'off';
            app.UIFigure.CloseRequestFcn = createCallbackFcn(app, @UIFigureCloseRequest, true);

            % Create TabGroup
            app.TabGroup = uitabgroup(app.UIFigure);
            app.TabGroup.Position = [14 18 1626 810];

            % Create MicrostructuresettingTab
            app.MicrostructuresettingTab = uitab(app.TabGroup);
            app.MicrostructuresettingTab.Title = 'Microstructure setting';

            % Create DisintegrationmeidumsettingPanel
            app.DisintegrationmeidumsettingPanel = uipanel(app.MicrostructuresettingTab);
            app.DisintegrationmeidumsettingPanel.Title = 'Disintegration meidum setting';
            app.DisintegrationmeidumsettingPanel.Position = [16 427 233 162];

            % Create DynamicViscosityPasLabel
            app.DynamicViscosityPasLabel = uilabel(app.DisintegrationmeidumsettingPanel);
            app.DynamicViscosityPasLabel.HorizontalAlignment = 'right';
            app.DynamicViscosityPasLabel.Position = [8 79 148 22];
            app.DynamicViscosityPasLabel.Text = 'Dynamic Viscosity [mPa*s]';

            % Create DynamicViscositymPasEditField
            app.DynamicViscositymPasEditField = uieditfield(app.DisintegrationmeidumsettingPanel, 'numeric');
            app.DynamicViscositymPasEditField.ValueDisplayFormat = '%5.3f';
            app.DynamicViscositymPasEditField.Editable = 'off';
            app.DynamicViscositymPasEditField.Position = [167 79 57 22];

            % Create WatertemperatureCLabel
            app.WatertemperatureCLabel = uilabel(app.DisintegrationmeidumsettingPanel);
            app.WatertemperatureCLabel.HorizontalAlignment = 'right';
            app.WatertemperatureCLabel.Position = [7 111 124 22];
            app.WatertemperatureCLabel.Text = 'Water temperature [C]';

            % Create WatertemperatureCEditField
            app.WatertemperatureCEditField = uieditfield(app.DisintegrationmeidumsettingPanel, 'numeric');
            app.WatertemperatureCEditField.Limits = [2 40];
            app.WatertemperatureCEditField.ValueDisplayFormat = '%.0f';
            app.WatertemperatureCEditField.ValueChangedFcn = createCallbackFcn(app, @WatertemperatureCEditFieldValueChanged, true);
            app.WatertemperatureCEditField.Position = [167 111 57 22];
            app.WatertemperatureCEditField.Value = 37;

            % Create SurfaceTensionmNmEditFieldLabel
            app.SurfaceTensionmNmEditFieldLabel = uilabel(app.DisintegrationmeidumsettingPanel);
            app.SurfaceTensionmNmEditFieldLabel.HorizontalAlignment = 'right';
            app.SurfaceTensionmNmEditFieldLabel.Position = [7 48 133 22];
            app.SurfaceTensionmNmEditFieldLabel.Text = 'Surface Tension [mN/m]';

            % Create SurfaceTensionmNmEditField
            app.SurfaceTensionmNmEditField = uieditfield(app.DisintegrationmeidumsettingPanel, 'numeric');
            app.SurfaceTensionmNmEditField.ValueDisplayFormat = '%5.3f';
            app.SurfaceTensionmNmEditField.Editable = 'off';
            app.SurfaceTensionmNmEditField.Position = [167 48 57 22];

            % Create ContactangledegEditFieldLabel
            app.ContactangledegEditFieldLabel = uilabel(app.DisintegrationmeidumsettingPanel);
            app.ContactangledegEditFieldLabel.HorizontalAlignment = 'right';
            app.ContactangledegEditFieldLabel.Position = [6 17 110 22];
            app.ContactangledegEditFieldLabel.Text = 'Contact angle [deg]';

            % Create ContactangledegEditField
            app.ContactangledegEditField = uieditfield(app.DisintegrationmeidumsettingPanel, 'numeric');
            app.ContactangledegEditField.Limits = [0 180];
            app.ContactangledegEditField.ValueDisplayFormat = '%4.2f';
            app.ContactangledegEditField.Position = [167 17 57 22];
            app.ContactangledegEditField.Value = 30;

            % Create MicrostructuresettingPanel
            app.MicrostructuresettingPanel = uipanel(app.MicrostructuresettingTab);
            app.MicrostructuresettingPanel.Title = 'Microstructure setting';
            app.MicrostructuresettingPanel.Position = [16 597 234 167];

            % Create CapillaryarraynumberEditFieldLabel
            app.CapillaryarraynumberEditFieldLabel = uilabel(app.MicrostructuresettingPanel);
            app.CapillaryarraynumberEditFieldLabel.HorizontalAlignment = 'right';
            app.CapillaryarraynumberEditFieldLabel.Position = [6 113 127 22];
            app.CapillaryarraynumberEditFieldLabel.Text = 'Capillary array number';

            % Create CapillaryarraynumberEditField
            app.CapillaryarraynumberEditField = uieditfield(app.MicrostructuresettingPanel, 'numeric');
            app.CapillaryarraynumberEditField.Limits = [0 1000];
            app.CapillaryarraynumberEditField.ValueDisplayFormat = '%.0f';
            app.CapillaryarraynumberEditField.Position = [177 113 48 22];
            app.CapillaryarraynumberEditField.Value = 50;

            % Create CapillarysizenumberEditFieldLabel
            app.CapillarysizenumberEditFieldLabel = uilabel(app.MicrostructuresettingPanel);
            app.CapillarysizenumberEditFieldLabel.HorizontalAlignment = 'right';
            app.CapillarysizenumberEditFieldLabel.Position = [5 83 121 22];
            app.CapillarysizenumberEditFieldLabel.Text = 'Capillary size number';

            % Create CapillarysizenumberEditField
            app.CapillarysizenumberEditField = uieditfield(app.MicrostructuresettingPanel, 'numeric');
            app.CapillarysizenumberEditField.Limits = [0 50];
            app.CapillarysizenumberEditField.ValueDisplayFormat = '%.0f';
            app.CapillarysizenumberEditField.Position = [177 83 48 22];
            app.CapillarysizenumberEditField.Value = 10;

            % Create Layernumber105000Label
            app.Layernumber105000Label = uilabel(app.MicrostructuresettingPanel);
            app.Layernumber105000Label.HorizontalAlignment = 'right';
            app.Layernumber105000Label.Position = [5 52 142 22];
            app.Layernumber105000Label.Text = 'Layer number (10 - 5000)';

            % Create LayernumberEditField
            app.LayernumberEditField = uieditfield(app.MicrostructuresettingPanel, 'numeric');
            app.LayernumberEditField.Limits = [10 5000];
            app.LayernumberEditField.ValueDisplayFormat = '%.0f';
            app.LayernumberEditField.Position = [177 52 48 22];
            app.LayernumberEditField.Value = 100;

            % Create Interlayerconnectivity50500Label
            app.Interlayerconnectivity50500Label = uilabel(app.MicrostructuresettingPanel);
            app.Interlayerconnectivity50500Label.HorizontalAlignment = 'right';
            app.Interlayerconnectivity50500Label.Position = [5 21 177 22];
            app.Interlayerconnectivity50500Label.Text = 'Interlayer connectivity (50 - 500)';

            % Create InterlayerconnectivityEditField
            app.InterlayerconnectivityEditField = uieditfield(app.MicrostructuresettingPanel, 'numeric');
            app.InterlayerconnectivityEditField.Limits = [50 500];
            app.InterlayerconnectivityEditField.ValueDisplayFormat = '%.0f';
            app.InterlayerconnectivityEditField.Position = [191 21 34 22];
            app.InterlayerconnectivityEditField.Value = 200;

            % Create CapillaryDistributionPanel
            app.CapillaryDistributionPanel = uipanel(app.MicrostructuresettingTab);
            app.CapillaryDistributionPanel.Title = 'Capillary Distribution';
            app.CapillaryDistributionPanel.Position = [262 15 510 748];

            % Create UIAxes12
            app.UIAxes12 = uiaxes(app.CapillaryDistributionPanel);
            title(app.UIAxes12, 'Radius')
            xlabel(app.UIAxes12, 'Bin')
            ylabel(app.UIAxes12, 'Radius [um]')
            zlabel(app.UIAxes12, 'Z')
            app.UIAxes12.LineWidth = 0.2;
            app.UIAxes12.Box = 'on';
            app.UIAxes12.FontSize = 10;
            app.UIAxes12.Position = [15 371 480 150];

            % Create UIAxes13
            app.UIAxes13 = uiaxes(app.CapillaryDistributionPanel);
            title(app.UIAxes13, 'Volume')
            xlabel(app.UIAxes13, 'Bin')
            ylabel(app.UIAxes13, 'Volume [um^3]')
            zlabel(app.UIAxes13, 'Z')
            app.UIAxes13.LineWidth = 0.2;
            app.UIAxes13.Box = 'on';
            app.UIAxes13.FontSize = 10;
            app.UIAxes13.Position = [15 199 480 150];

            % Create UIAxes14
            app.UIAxes14 = uiaxes(app.CapillaryDistributionPanel);
            title(app.UIAxes14, 'Saturation time')
            xlabel(app.UIAxes14, 'Bin')
            ylabel(app.UIAxes14, 'Time [ms]')
            zlabel(app.UIAxes14, 'Z')
            app.UIAxes14.LineWidth = 0.2;
            app.UIAxes14.Box = 'on';
            app.UIAxes14.FontSize = 10;
            app.UIAxes14.Position = [15 27 480 150];

            % Create UIAxes11
            app.UIAxes11 = uiaxes(app.CapillaryDistributionPanel);
            title(app.UIAxes11, 'Pore radius probability density function (4 sigma)')
            xlabel(app.UIAxes11, 'Radius [um]')
            ylabel(app.UIAxes11, 'Probablity')
            zlabel(app.UIAxes11, 'Z')
            app.UIAxes11.LineWidth = 0.2;
            app.UIAxes11.Box = 'on';
            app.UIAxes11.FontSize = 10;
            app.UIAxes11.Position = [15 543 480 171];

            % Create MicrostructuredetailsPanel
            app.MicrostructuredetailsPanel = uipanel(app.MicrostructuresettingTab);
            app.MicrostructuredetailsPanel.Title = 'Microstructure details';
            app.MicrostructuredetailsPanel.Position = [16 93 233 99];

            % Create CapillaryheightumEditFieldLabel
            app.CapillaryheightumEditFieldLabel = uilabel(app.MicrostructuredetailsPanel);
            app.CapillaryheightumEditFieldLabel.HorizontalAlignment = 'right';
            app.CapillaryheightumEditFieldLabel.Position = [7 44 116 22];
            app.CapillaryheightumEditFieldLabel.Text = 'Capillary height (um)';

            % Create CapillaryheightumEditField
            app.CapillaryheightumEditField = uieditfield(app.MicrostructuredetailsPanel, 'numeric');
            app.CapillaryheightumEditField.ValueDisplayFormat = '%5.2f';
            app.CapillaryheightumEditField.Editable = 'off';
            app.CapillaryheightumEditField.Position = [167 44 57 22];

            % Create StructureheightmmEditFieldLabel
            app.StructureheightmmEditFieldLabel = uilabel(app.MicrostructuredetailsPanel);
            app.StructureheightmmEditFieldLabel.HorizontalAlignment = 'right';
            app.StructureheightmmEditFieldLabel.Position = [7 15 122 22];
            app.StructureheightmmEditFieldLabel.Text = 'Structure height (mm)';

            % Create StructureheightmmEditField
            app.StructureheightmmEditField = uieditfield(app.MicrostructuredetailsPanel, 'numeric');
            app.StructureheightmmEditField.ValueDisplayFormat = '%5.2f';
            app.StructureheightmmEditField.Editable = 'off';
            app.StructureheightmmEditField.Position = [167 15 57 22];

            % Create AssingdatainworkspaceButton
            app.AssingdatainworkspaceButton = uibutton(app.MicrostructuresettingTab, 'push');
            app.AssingdatainworkspaceButton.ButtonPushedFcn = createCallbackFcn(app, @AssingdatainworkspaceButtonPushed, true);
            app.AssingdatainworkspaceButton.Position = [40 36 190 32];
            app.AssingdatainworkspaceButton.Text = 'Assing data in workspace';

            % Create Panel
            app.Panel = uipanel(app.MicrostructuresettingTab);
            app.Panel.Title = 'Panel';
            app.Panel.Position = [789 15 815 748];

            % Create UIAxes21
            app.UIAxes21 = uiaxes(app.Panel);
            title(app.UIAxes21, 'Capillary saturation')
            xlabel(app.UIAxes21, 'Bin')
            ylabel(app.UIAxes21, 'Layer')
            zlabel(app.UIAxes21, 'Z')
            app.UIAxes21.LineWidth = 0.2;
            app.UIAxes21.Box = 'on';
            app.UIAxes21.FontSize = 10;
            app.UIAxes21.Position = [8 200 525 517];

            % Create UIAxes22
            app.UIAxes22 = uiaxes(app.Panel);
            title(app.UIAxes22, 'Saturation (accumulated)')
            xlabel(app.UIAxes22, 'Saturation')
            ylabel(app.UIAxes22, 'Layer')
            zlabel(app.UIAxes22, 'Z')
            app.UIAxes22.LineWidth = 0.2;
            app.UIAxes22.Box = 'on';
            app.UIAxes22.FontSize = 10;
            app.UIAxes22.Position = [554 200 236 517];

            % Create UIAxes23
            app.UIAxes23 = uiaxes(app.Panel);
            title(app.UIAxes23, 'Liquid Ingress')
            xlabel(app.UIAxes23, 'Time (ms)')
            ylabel(app.UIAxes23, 'Depth (mm)')
            zlabel(app.UIAxes23, 'Z')
            app.UIAxes23.LineWidth = 0.2;
            app.UIAxes23.Box = 'on';
            app.UIAxes23.FontSize = 10;
            app.UIAxes23.Position = [436 11 354 182];

            % Create LoadstructureButton
            app.LoadstructureButton = uibutton(app.Panel, 'push');
            app.LoadstructureButton.ButtonPushedFcn = createCallbackFcn(app, @LoadstructureButtonPushed, true);
            app.LoadstructureButton.Position = [225 151 173 30];
            app.LoadstructureButton.Text = 'Load structure';

            % Create ControlboxPanel
            app.ControlboxPanel = uipanel(app.Panel);
            app.ControlboxPanel.Title = 'Control box';
            app.ControlboxPanel.Position = [40 35 381 104];

            % Create StartButton_2
            app.StartButton_2 = uibutton(app.ControlboxPanel, 'push');
            app.StartButton_2.ButtonPushedFcn = createCallbackFcn(app, @StartButton_2Pushed, true);
            app.StartButton_2.Position = [13 12 83 32];
            app.StartButton_2.Text = 'Start';

            % Create StopButton_2
            app.StopButton_2 = uibutton(app.ControlboxPanel, 'push');
            app.StopButton_2.ButtonPushedFcn = createCallbackFcn(app, @StopButton_2Pushed, true);
            app.StopButton_2.Position = [104 12 83 32];
            app.StopButton_2.Text = 'Stop';

            % Create SaveButton
            app.SaveButton = uibutton(app.ControlboxPanel, 'push');
            app.SaveButton.Position = [286 12 83 32];
            app.SaveButton.Text = 'Save';

            % Create ResetButton
            app.ResetButton = uibutton(app.ControlboxPanel, 'push');
            app.ResetButton.ButtonPushedFcn = createCallbackFcn(app, @ResetButtonPushed, true);
            app.ResetButton.Position = [195 12 83 32];
            app.ResetButton.Text = 'Reset';

            % Create TimemsEditFieldLabel
            app.TimemsEditFieldLabel = uilabel(app.ControlboxPanel);
            app.TimemsEditFieldLabel.HorizontalAlignment = 'right';
            app.TimemsEditFieldLabel.Position = [222 55 59 22];
            app.TimemsEditFieldLabel.Text = 'Time (ms)';

            % Create TimemsEditField
            app.TimemsEditField = uieditfield(app.ControlboxPanel, 'numeric');
            app.TimemsEditField.Limits = [0 Inf];
            app.TimemsEditField.ValueDisplayFormat = '%.1f';
            app.TimemsEditField.Editable = 'off';
            app.TimemsEditField.FontSize = 16;
            app.TimemsEditField.Position = [286 53 81 26];

            % Create ThreadholdLabel
            app.ThreadholdLabel = uilabel(app.ControlboxPanel);
            app.ThreadholdLabel.HorizontalAlignment = 'right';
            app.ThreadholdLabel.Position = [18 54 147 22];
            app.ThreadholdLabel.Text = 'Liquid front threadhold (%)';

            % Create LiquidfrontthreadholdEditField
            app.LiquidfrontthreadholdEditField = uieditfield(app.ControlboxPanel, 'numeric');
            app.LiquidfrontthreadholdEditField.Limits = [0 100];
            app.LiquidfrontthreadholdEditField.ValueDisplayFormat = '%.0f';
            app.LiquidfrontthreadholdEditField.Position = [169 54 37 22];
            app.LiquidfrontthreadholdEditField.Value = 70;

            % Create ColormapDropDownLabel
            app.ColormapDropDownLabel = uilabel(app.Panel);
            app.ColormapDropDownLabel.HorizontalAlignment = 'right';
            app.ColormapDropDownLabel.Position = [39 155 58 22];
            app.ColormapDropDownLabel.Text = 'Colormap';

            % Create ColormapDropDown
            app.ColormapDropDown = uidropdown(app.Panel);
            app.ColormapDropDown.Items = {'parula', 'jet', 'copper', 'bone', 'hot'};
            app.ColormapDropDown.Position = [112 155 100 22];
            app.ColormapDropDown.Value = 'parula';

            % Create PoresettingPanel
            app.PoresettingPanel = uipanel(app.MicrostructuresettingTab);
            app.PoresettingPanel.Title = 'Pore setting';
            app.PoresettingPanel.Position = [16 215 233 207];

            % Create MeanradiusumEditFieldLabel
            app.MeanradiusumEditFieldLabel = uilabel(app.PoresettingPanel);
            app.MeanradiusumEditFieldLabel.HorizontalAlignment = 'right';
            app.MeanradiusumEditFieldLabel.Position = [12 108 98 22];
            app.MeanradiusumEditFieldLabel.Text = 'Mean radius [um]';

            % Create MeanradiusumEditField
            app.MeanradiusumEditField = uieditfield(app.PoresettingPanel, 'numeric');
            app.MeanradiusumEditField.Limits = [0 Inf];
            app.MeanradiusumEditField.ValueDisplayFormat = '%5.2f';
            app.MeanradiusumEditField.Position = [163 108 61 22];
            app.MeanradiusumEditField.Value = 5;

            % Create StandarddeviationumEditFieldLabel
            app.StandarddeviationumEditFieldLabel = uilabel(app.PoresettingPanel);
            app.StandarddeviationumEditFieldLabel.HorizontalAlignment = 'right';
            app.StandarddeviationumEditFieldLabel.Position = [11 80 132 22];
            app.StandarddeviationumEditFieldLabel.Text = 'Standard deviation [um]';

            % Create StandarddeviationumEditField
            app.StandarddeviationumEditField = uieditfield(app.PoresettingPanel, 'numeric');
            app.StandarddeviationumEditField.Limits = [0 Inf];
            app.StandarddeviationumEditField.ValueDisplayFormat = '%5.2f';
            app.StandarddeviationumEditField.Position = [163 80 61 22];
            app.StandarddeviationumEditField.Value = 1;

            % Create PoresizedistributionDropDownLabel
            app.PoresizedistributionDropDownLabel = uilabel(app.PoresettingPanel);
            app.PoresizedistributionDropDownLabel.HorizontalAlignment = 'right';
            app.PoresizedistributionDropDownLabel.Position = [14 156 117 22];
            app.PoresizedistributionDropDownLabel.Text = 'Pore size distribution';

            % Create PoresizedistributionDropDown
            app.PoresizedistributionDropDown = uidropdown(app.PoresettingPanel);
            app.PoresizedistributionDropDown.Items = {'Normal distribution', 'Function load'};
            app.PoresizedistributionDropDown.Position = [49 135 154 22];
            app.PoresizedistributionDropDown.Value = 'Normal distribution';

            % Create UpdateButton
            app.UpdateButton = uibutton(app.PoresettingPanel, 'push');
            app.UpdateButton.ButtonPushedFcn = createCallbackFcn(app, @UpdateButtonPushed, true);
            app.UpdateButton.Position = [20 9 201 30];
            app.UpdateButton.Text = 'Update';

            % Create Poreisotropicity0110Label
            app.Poreisotropicity0110Label = uilabel(app.PoresettingPanel);
            app.Poreisotropicity0110Label.HorizontalAlignment = 'right';
            app.Poreisotropicity0110Label.Position = [12 52 143 22];
            app.Poreisotropicity0110Label.Text = 'Pore isotropicity (0.1 - 10)';

            % Create PoreisotropicalityEditField
            app.PoreisotropicalityEditField = uieditfield(app.PoresettingPanel, 'numeric');
            app.PoreisotropicalityEditField.Limits = [0.1 10];
            app.PoreisotropicalityEditField.ValueDisplayFormat = '%3.2f';
            app.PoreisotropicalityEditField.Position = [176 52 48 22];
            app.PoreisotropicalityEditField.Value = 1;

            % Create demoTab
            app.demoTab = uitab(app.TabGroup);
            app.demoTab.Title = 'demo';

            % Create UIAxes
            app.UIAxes = uiaxes(app.demoTab);
            xlabel(app.UIAxes, 'Seconds')
            ylabel(app.UIAxes, 'Random number calculated')
            app.UIAxes.XTickLabelRotation = 0;
            app.UIAxes.YTickLabelRotation = 0;
            app.UIAxes.ZTickLabelRotation = 0;
            app.UIAxes.Box = 'on';
            app.UIAxes.XGrid = 'on';
            app.UIAxes.YGrid = 'on';
            app.UIAxes.Position = [19 234 1019 451];

            % Create StartButton
            app.StartButton = uibutton(app.demoTab, 'push');
            app.StartButton.ButtonPushedFcn = createCallbackFcn(app, @StartButtonPushed, true);
            app.StartButton.Position = [66 36 100 22];
            app.StartButton.Text = 'Start';

            % Create StopButton
            app.StopButton = uibutton(app.demoTab, 'push');
            app.StopButton.ButtonPushedFcn = createCallbackFcn(app, @StopButtonPushed, true);
            app.StopButton.Position = [217 36 100 22];
            app.StopButton.Text = 'Stop';

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = LiquidIngressSimulator_exported

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
    end
end