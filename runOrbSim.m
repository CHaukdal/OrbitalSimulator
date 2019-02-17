%% runOrbSim.m    
% Version 1.1.1
% 
% Simulates the motion of the solar system. Designed for analysis of
% system of satellites in orbit around Jupiter's Moons Europa and Ganymede.
% Easily modified for inclusion of other bodies and satellites.
% 
% CALLS FUNCTIONS:
% solveOrbSim.m
% plotOrbSim.m
% 
% INSTRUCTIONS:
% 1. Read descriptions and edit information as needed in sections 
%    "Information Cell Array for Selected Bodies", "Run Options", 
%    "Advanced Run Options", and "Plot Options"
% 2. Click Run
% 
% IMPORTANT OUTPUTS:
% t                 Time vector, can be saved if desired
% y                 Positions and velocites matrix, can be saved
% [Plot]            Chosen plot and data body orbits, can be saved
% [Analysis]        Chosen plots and data from analysis, can be saved 
% 
% 
% Cy Haukdal
% Haukdal@vt.edu
% 
% Copyright (c) 2019, Cy Haukdal
% All rights reserved.

clear
clc
close all


%% Information Cell Array for Selected Bodies
% change this structure to change the bodies simulated. Currently includes
% sun, 8 planets, and pluto.
% mass = body mass [kg]
% name = chosen body descriptor
% mdb = mean distance from solar system barycenter [km]
% xEph = x ephemeris position, cartesian coordinate from solar system
% baycenter [km]
% yEph = y "    "
% zEph = z "    "
% vxEph = x ephemeris velocity, cartesian coordiante from solar system
% barycenter [km/s]
% yEph = y "    "
% zEph = z "    "
%
% In the form bodyInfo = {mass1, 'name1', mdb1,...
%                       xEph1, yEph1, zEph1,... 
%                       vxEph1, vyEph1, vzEph1;...
%                       mass2, 'name2', mdb2,... 
%                       xEph2, yEph2, zEph2,...
%                       vxEph2, vyEph2, vzEph2;... etc.};

% Ephemeris data can be sourced from JPL HORIZONS System, availible at: 
% https://ssd.jpl.nasa.gov/?horizons
% Example sourced for 2019-feb-14 00:00:00.0000 Zulu time
% This example can be used for any time after 2019-feb-14 but updating
% ephemeris data will increase the accuracy of the simulation
%
% Initial sim time can be set to a later date than 2019-feb-14 (or date of
% new emhemeris data) by changing iniTime to the number of earth years
% since the data's timestamp
% For example, if a start date of 2019-feb-20 at 12:00:00.0000 Zulu is
% desired, iniTime would be 6 days 12 hrs = 6.5 days / 365.25 days per year
% or iniTime = .017796
% It is STRONGLY reccomended the ephemeris data is changed to the desired
% start date, as using this method still requires the function to
% simulate time from 0 regardeless of iniTime chosen and will take
% additional time

% Change these [earth years]
iniTime = 0;
tStep = .01; 
simLength = 300; 

% Leave these be
simTime = simLength+iniTime;
tVec = [iniTime tStep simTime];

% Change this
% ONLY ONE body with initial position zero and/or initial scalar velocity 
% zero allowed, and you must put it in the fisrt row of the cell array.
% Having zero bodies at pos = [0 0 0], vel = [0 0 0] will work.
% Other bodies like pos = [0 0 1], vel = [0 1 0] will work.
% Having multiple bodies at any combo of |pos| = 0 |vel| = 0, 
% |pos| =/= 0 |vel| = 0, or |pos| = 0 |vel| =/= 0 will NOT work.
bodyInfo = {1.9885e30, 'Sun', 0,...
    0, 0, 0,...
    0, 0, 0;...
    3.3011e23, 'Mercury', 57909050,...
    5.002621341152667E+07, 1.091922423877655E+07, -3.811228705514119E+06,...
    -1.876683641789480E+01, 4.999144734136460E+01, 5.805451887547726E+00;...
    4.8675e24, 'Venus', 108208000,...
    -9.246153485952877E+07, -5.491702725827458E+07, 4.550011361320481E+06,...
    1.791768733299224E+01, -3.009509677131724E+01, -1.447413967883188E+00;...
    5.97237e24, 'Earth', 149095000,...
    -1.207683279357202E+08, 8.640890342186935E+07, -1.044142235558853E+04,...
    -1.768792864922047E+01, -2.443459403117340E+01, 7.043517675402455E-04;...
    6.4171e23, 'Mars', 227939200,...
    9.273287397753367E+07, 2.068930938569493E+08, 2.024910379655212E+06,...
    -2.118094270812971E+01, 1.203448521992707E+01, 7.718505294098783E-01;...
    1.8982e27, 'Jupiter', 778.57e6,...
    -2.738907460398777E+08, -7.488014414760873E+08, 9.232225403508514E+06,...
    1.211290578888045E+01, -3.864117898988086E+00, -2.549106545797410E-01;...
    5.6834e26, 'Saturn', 1433.53e6,...
    3.269438158221755E+08, -1.467534600645217E+09, 1.250235259336644E+07,...
    8.896967013344858E+00, 2.070824816589930E+00, -3.905657555353882E-01;...
    8.6810e25, 'Uranus', 2742e6,...
    2.531811882268407E+09, 1.554370058840762E+09, -2.702695884676933E+07,...
    -3.612824012488903E+00, 5.486081641811738E+00, 6.713412711299993E-02;...
    1.02413e23, 'Neptune', 4.50e9,...
    4.340571784849433E+09, -1.099031244421821E+09, -7.740037958078402E+07,...
    1.297727562652028E+00, 5.301832905882191E+00, -1.383835127416295E-01;...
    1.303e22, 'Pluto', 5.90638e9,...
    1.797258000104297E+09, -4.715485121966079E+09, -1.528978071004462E+07,...
    5.182764472252960E+00, 7.690111858234998E-01, -1.588510432535072E+00};


%% Run Options
% Settings for system solver

% For convenience (do not change):
nBodies = size(bodyInfo,1);

% Choose which bodies will exert gravitational influence in sim:
% For example, sun only would be selBodies = 1; (This saves time but is
% less accurate and not accurate if satellites of a planet are simulated)
% This simulates all bodies:
selVec = 1:nBodies;

% Do you want to save t (time step vector), y (corresponding positions and
% velocities), and bodyInfo to text files? If so, saveData = 1, if not, saveData = 0.
% If saveData = 1, choose dataNameY, dataNameT, and dataLocation as seen 
% below. If saveData = 0, dataNameY, dataNameT, and dataLocation will be ignored.
saveData = 1;
% saveData = 0;
dataNameY = 'exampleNameY';
dataNameT = 'exampleNameT';
dataNameB = 'exampleNameBodyInfo';
dataLocation = 'C:\test'; 

%% Advanced Run Options (can be left alone):
% Changes ode solver settings

% Relative error tolerance settings (Matlab default 1e-3, can be left alone)
% Quick run (may diverge):
% relErrTol = 1; 
% Accurate run suggestion (may be dropped further)
relErrTol = 1e-5;

% Absolute error tolerance settings (Matalb default 1e-6, runs slowly)
% Split due to expected magnitude differences in position and velocity 
% Normalized RMS error * coeffecient used
% Note: Currently only position is normalized, velocity is standard RMS. 
% Velocity will be updated to normalized RMS in next version.
% Quick run (may diverege or produce inaccurate results):
% posCoeff = 1e-2;
% velCoeff = 1e-4;
% Accurate run suggestions (may be dropped further)
posCoeff = 1e-6;
velCoeff = 1e-8; %1e-9 would be nice but runs slowly

% Additional ODE controls such as 'NormControl' or 'Stats' can be applied
% manually in line 124 of solve_system.m, 'OutputFcn' is already defined as 
% a local function and can be modified at line 189 of solve_system.m

%% Plot Options
% Settings for plotting function

% Do you want generate a plot and/or save the orbit plot? If so, 
% plotOn = 1, if not, plotOn = 0
plotOn = 1;
% plotOn = 0;

% Do you want to display the plot generated? If so, dispPlot = 1, if not,
% dispPlot = 0.
dispPlot = 1;
% dispPlot = 0;

% Do you want to animate the plot generated? If so, aniPlot = 1, if
% not, aniPlot = 0
% aniPlot = 1;
aniPlot = 0;

% Ratio of time steps per plot point (increase to speed up animation)
plotStep = 10;
% plotStep = 100;

% Do you want to see the bodies at their proper size? If so, enter their
% radii in bodyRadii as seen below. Otherwise, set bodyRadii = 0;
bodyRadii = [2440 6052 6378 3397 71492 60269 25559 24766 1150];
% bodyRadii = 0;

% Do you want to log scale the bodies? 0 for an accurate scale, 1 so you
% can see them all at the same time while zoomed out.
scaleBodies = 1;
% scalebodies = 0;

% Do you want to save the orbits produced? If so, savePlot = 1, if not
% savePlot = 0. If savePlot = 1, choose fileName and fileLocation as seen 
% below, otherwise fileName and fileLocation will be ignored. 3D plot only
% support saving as .fig file and is done by default. 
% savePlot = 1;
savePlot = 0;
fileName = 'exampleName';
fileLocation = 'C:\test'; 

% If you want to be fancy, you can color bodies for plotting 
% here. If so, enter the rgb percentage for each body as a fraction seen 
% below. If not, set plotColors = 0. Currently set to solar system.
plotColors = {[.31, .31, .31], [.4, .35, .25], [0, .20, .7],...
    [.99, .5, .4], [.38, .36, .32], [.25, .22, .16], [.29, .42, .46]...
    [.05, .43,.46], [.27, .25, .22]};
% plotColors = 0;

% Conglomerated plot settings (do not chage):
plotSettings = [nBodies dispPlot aniPlot plotStep scaleBodies savePlot];

%% Running Solver

errVec = [relErrTol posCoeff velCoeff];
[t, y] = solveOrbSim(bodyInfo, selVec, tVec, errVec); 

%% Plotting and Saving if Desired

if saveData == 1
    SavePathY = strcat(dataLocation, '\', dataNameY, '.txt');
    SavePathT = strcat(dataLocation, '\', dataNameT, '.txt');
%     SavePathB = strcat(dataLocation, '\', dataNameB, '.mat')
    dlmwrite(SavePathY, y)
    dlmwrite(SavePathT, t)
%     save(SavePathB, bodyInfo, 'c')
end

if plotOn == 1
    filePath = strcat(fileLocation, '\', fileName);
    plotOrbSim(y, t, plotSettings, bodyRadii, plotColors, tVec, filePath)
end


