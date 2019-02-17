function plotOrbSim(y, t, plotSettings, bodyRadii, plotColors, tVec, filePath)
%% function plot_system(y, t, plotSettings, bodyRadii, plotColors, tVec, filePath)
% Version 1.1.1 
% 
% Produces and saves visual plot of body orbits and trajectories and animates plot
%
% Called by script: run_system.m       
% 
% INPUTS:
% y                 Position and velocity matrix of bodies over time
% t                 Time vector corresponding to y
% plotSettings      Binary vector containing:
%    nBodies           Number of bodies analyzed by sim
%    dispPlot          Display plot or no
%    aniPlot           Animate plot or no
%    plotStep          Ratio of time steps per plot point               
%    scaleBodies       Binary to scale bodies or not
%    savePlot          Binary to save plot or not
% bodyRadii         Vector of body radii, or 0
% plotColors        Plot colors, or 0
% tVec              Vector containing timing info
% filePath          Path and name type of plot to be saved if selected 
% 
% OUTPUTS:
% Plot              Animated plot of bodies in orbit.  Plot of orbits
%                   can be saved if desired.
% 
% 
% Cy Haukdal
% Haukdal@vt.edu
% 
% Copyright (c) 2019, Cy Haukdal
% All rights reserved.


%% Parsing Inputs

iniTime = tVec(1);
tStep = tVec(2);
iniStep = iniTime/tStep+1;
nBodies = plotSettings(1);
dispPlot = plotSettings(2);
aniPlot = plotSettings(3);
plotStep = plotSettings(4);
scaleBodies = plotSettings(5);
savePlot = plotSettings(6);

%% Plotting and Saving

% Designed for 1080p display, alter position vector for other screen resolutions 
f = figure('Name', 'Simulation Results', 'NumberTitle', 'off', 'renderer',...
    'openGL', 'position', [1920/4 1080/4 1920/2 1080/2]);
if ~iscell(plotColors)
    plotColors = cell(1,nBodies);
    for iter = 1:nBodies
        plotColors{iter} = [0,0,0];
    end
end

% Display plot check
if dispPlot == 0
   set(0,'DefaultFigureVisible','off'); 
end

% Plotting body trajectories
for iter = 1:nBodies-1 % doesnt plot sun
    % Body position minus moving sun position for stability
    f(iter) = plot3(y(iniStep:end,iter*6+1)-y(iniStep:end,1), y(iniStep:end,iter*6+2)...
        -y(iniStep:end,2), y(iniStep:end,iter*6+3)-y(iniStep:end,3),...
        'linewidth', 0.3, 'color', plotColors{iter});
    hold on
end

% Axis labels and scaling and background to black (we are in space after all)
xlabel('x [km]')
ylabel('y [km]')
zlabel('z [km]')
axis equal

% If you are in this deep and want a black background:
% set(gca, 'Color', 'k')

% Save orbit plot if requested
if savePlot == 1
    saveStatus = msgbox(strcat('Plot file saved: ', ' ', filePath, '.fig'));
    saveas(f, filePath, 'fig')    
end

if aniPlot == 1
    % Animated lot options from user
    placeHolder = ones(1,nBodies-1)*5000; % Matlab gets mad without this
    if bodyRadii == 0
        sz = log(placeHolder./.00001);  
    elseif scaleBodies == 0
        sz = bodyRadii*2;
    else
        sz=log(bodyRadii*2./.00001); 
    end

    % Animating body positions
    for iter1 = iniStep:plotStep:length(t)
        for iter2 = 1:nBodies-1
            % Body position minus moving sun position for stability
            f(iter2) = plot3(y(iter1,iter2*6+1)-y(iter1,1), ...
                y(iter1,iter2*6+2)-y(iter1,2), y(iter1,iter2*6+3)-y(iter1,3),...
                '.', 'color', plotColors{iter2}, 'markersize', sz(iter2));
        end
        pause(.01) % Required for animation
        delete(f) % Required for animation
    end
end

%     Closing message box
if savePlot == 1
    if ishandle(saveStatus)
        delete(saveStatus)
    end
end

% Re-setting plot default setting in matlab to visible
set(0,'DefaultFigureVisible','on');

% funciton end
end
