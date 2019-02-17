%% analyzeOrbSim.m    
% Version 1.1.1
% 
% Takes position vs time data produced by run_system.m and produces plots
% and raw data of communicaiton and sun seeing windows
% 
% INSTRUCTIONS:
% 1. Load data or use data from workspace for t and y in "Data Input" section
% 2. Edit analysis settings as needed in "Settings" section
% 3. Click run
% 
% INPUTS:
% t                 Time step vector 
% y                 Position and veloicty matrix
% bodyInfo          information about bodies 
% 
% IMPORTANT OUTPUTS:
% bodyWindows       Array where rows are pairing of bodies noted in column one,
%                   columns 2-n are indicies of t where bodies can see eachother
% bodyPlot(n)       Binary plots of window vs time for each pairing selected
% distPlot(n)       Plots of distance between selected pairings over time
% comboPlot(n)      Binary plots of times when bodies are visible to
%                   eachother and within a given "communication" distance
% commWindows       Array where rows are pairing of bodies noted in column one,
%                   columns 2-n are indicies of t where bodies can see
%                   eachother and are within a set communcation distance
% 
% 
% Cy Haukdal
% Haukdal@vt.edu
% 
% Copyright (c) 2019, Cy Haukdal
% All rights reserved.

clc

% [NOTE:] Indicates a future change will be made

%% Data Input
% Load data or use workspace data

tStep = .01;
% NOTE: need to streamline this input

% Is you data in the workspace and named t (simTime/tStep x 1 double), y
% (simTime/tStep x nBodies*6 double), and bodyInfo (nBodies x 9 cell array)
% If so, importData = 0. If your data is not in the workspace and needs to
% be imported, importData = 1. Using workspace varaibles initially or after
% data is loaded for subsequent runs is preferred for speed.
% importData = 1;
importData = 0;

% If your data is not in the workspace and you need to load the data, enter
% the file path in loadPath and the names of the files corresponding to t,
% y, and bodyInfo in tName, yName, and bName below.
loadPath = 'C:\test';
tName = 'exampleNameT';
yName = 'exampleNameY';
bName = 'exampleNameBodyInfo';

%% Settings

% Do you want to check for line of sight, communications widow, or both?
% sightLine and/or commLine = 1 will produce data for each choice, = 0 will
% not. 
sightLine = 1;
% sightLine = 0;
commLine = 1;
% commLine = 0;

% Do you want to plot line of sight, communications windows, or both? =1
% will plot, =0 will not
losPlot = 1;
% losPlot = 0;
commPlot = 1;
% commPlot = 0;

% Body pairs to be checked. Format bodyPairs, an n x 2 matrix of n pairs of
% bodies to be checked. Each body has an associated index, as inputted by 
% the user in bodyInfo. Example given below corresponding to example 
% bodyinfo in runOrbSim that checks sun to earth, earth to mars, and sun to jupiter.
bodyPairs = [1 2; 4 5; 1 6];
% bodyPairs = [4 5];

% Enter the radii of all bodies in the form seen below:
% NOTE: need to streamline this input
bodyRadii = [695508 2440 6052 6378 3397 71492 60269 25559 24766 1150]*1e3; %[m]

% If you are checking communications, what is your maximum
% communication range? Set commRange to max range in km
commRange = 100e9; %[km]

% Set a tolerance in km in which a communications singal will be disrupted
% past the surface of the body. For example, if a communications signal 
% from a Mars orbiting satellite must pass at least 100 km above the
% surface of mars to not be blocked by the atmosphere, sigTol = 100.
sigTol = 100; %[km]

%% Data Import and Parsing

% Loading data and formatting if needed
if importData == 1
    tTab = readtable(strcat(loadPath, '\', tName, '.txt'));
    t = table2array(tTab);
    yTab = readtable(strcat(loadPath, '\', yName, '.txt'));
    y = table2array(yTab);
%     NOTE: this needs to be changed for saving a cell array
%     bodyInfoTab = readtable(strcat(loadPath, '\', bName, '.txt'));
%     bodyInfo = table2array(bodyInfoTab);
end

% Extracting needed y data
nBodies = size(bodyInfo,1);
bodiesVec = 1:nBodies;
numPoints = size(y,1);

% Trimming y to only include position data
yPos = zeros(numPoints, nBodies*3);
yPosNorm = zeros(numPoints, nBodies);
for iter1 = 0:nBodies-1
    yPos(:, iter1*3+1:iter1*3+3) = y(:, iter1*6+1:iter1*6+3);
    
end

% Determining the max and min positions of each body from the barycenter
% Note: is there a way to get rid of this nested for loop
posMax = zeros(1, nBodies);
posMin = posMax;
for iter1 = 0:nBodies-1
    for iter2 = 1:numPoints
        yPosNorm(iter2,iter1+1) = norm(yPos(iter2,iter1*3+1:iter1*3+3));
    end
    posMax(iter1+1) = max(yPosNorm(:,iter1+1));
    posMin(iter1+1) = min(yPosNorm(:,iter1+1));
end


%% Checking for Line-of-Sight

% Position vector between bodies, line of sight matrix, comm matrix, and 
% scalar position delta are preallocated
% posVecBdys = zeros(numPoints,3);
% unitPosVecBdys = zeros(numPoints,3);
% posVecCheck = zeros(numPoints,3);
% unitPosVecCheck = zeros(numPoints,3);
numPairs = size(bodyPairs,1);
losCheck = NaN(numPoints, nBodies, numPairs);
commCheck = NaN(numPoints, nBodies, numPairs);

% Analysis Waitbar
% anaTime = 0;
% anaStatus = waitbar(anaTime, 'Analysis Progress');
% NOTE: This slows down the run considerably, need to determine the correct implementation

% Each body pair is checked 
% NOTE: Need to put this into a subfunctions for each comm, line of sight,
% and both to eliminate nested ifs and save computational time
for iter1 = 1:numPairs
    
    % Select bodies that are not currently in pair to check against
    bod2Check = find(bodiesVec ~= bodyPairs(iter1,1) & bodiesVec ~= bodyPairs(iter1,2));
    
    % Min and Max distances from barycenter to most extreme point of body
    % pair data
    minBdyDist = min(posMin(bodyPairs(iter1,1)), posMin(bodyPairs(iter1,2)));
    maxBdyDist = max(posMax(bodyPairs(iter1,1)), posMax(bodyPairs(iter1,2)));

    % Looping over each time step 
    for iter2 = 1:numPoints
        
        % Determining the position vector between body pair at each
        % time step
        posVecBdys = yPos(iter2, bodyPairs(iter1,2)*3-2:bodyPairs(iter1,2)*3)-...
            yPos(iter2, bodyPairs(iter1,1)*3-2:bodyPairs(iter1,1)*3);
        
        % Unit position vector between bodies
        unitPosVecBdys = posVecBdys/norm(posVecBdys);  
        
        % Looping over bodies to check against
        for iter3 = bod2Check      
            
            % Position vector between first body in pair and body to be
            % checked against
            posVecCheck = yPos(iter2, bodyPairs(iter1,2)*3-2:bodyPairs(iter1,2)*3)...
                - yPos(iter2, bodyPairs(iter1,1)*3-2:bodyPairs(iter1,1)*3);
            
            % Position vector between bodies scaled to same length as
            % position vector between first body in pair and checking body
            posVecBdysAdj = unitPosVecBdys*norm(posVecCheck);
            
            % Tip distance of scaled bodies position vector and body 1 to
            % checking body position vector
            posDelta = norm(posVecBdysAdj - yPos(iter2,iter3*3-2:iter3*3));
            
            % If min distance of check body from barycenter is less than 
            % the min distnce of both the two bodies to be checked from the
            % barycenter, or the max distnace of check body is greater than
            % the max distance of both the two bodies, it cannot interfere
            % NOTE: assumes center of sun is barycenter. This is NOT
            % ACCURATE for low-orbiting satellites. Need to update.
            
            % Logical to determine if checking body is outside the orbits 
            % of the the pair bodies
            distLogical = (posMin(iter3) < minBdyDist) || (posMax(iter3) > maxBdyDist);
            
            % If line of sight is to be checked
            if sightLine == 1
                
                % Checking if line of sight does not pass through checking body
                % or checking body orbit is not within the orbits of body pair
                if (posDelta >= bodyRadii(iter3)) || distLogical
                    losCheck(iter2, iter3, iter1) = 1;
                    
                % If pair position vector is within radius of body
                elseif posDelta <= bodyRadii(iter3)
                    losCheck(iter2, iter3, iter1) = 0;
                end
            end    
            
            % If comm capability is to be checked
            if commLine == 1
                
                % Checking if line of sight does not pass through checking body
                % plus the signal tolerance or checking bosy orbit is not
                % within the orbits of the body pair, and within
                % communicaiton distance
                if ((posDelta >= bodyRadii(iter3)+sigTol) || distLogical) && (norm(posVecBdys) <= commRange)
                    commCheck(iter2, iter3, iter1) = 1;
                    
                % If pair position vector is within radius of body + signal
                % tolerance or outside of comm range
                elseif (posDelta <= bodyRadii(iter3)+sigTol) || (norm(posVecBdys) >= commRange)
                    commCheck(iter2, iter3, iter1) = 0;
                end
            end           
        end
        % Updating waitbar
%             if ishandle(anaStatus)
%                 anaTime = iter2*iter1+(iter1-1)*iter2;
%                 totAnaTime = numPairs*numPoints;
%                 waitbar(anaTime/totAnaTime,anaStatus)
%             end
    end
end

% if ishandle(anaStatus)
%     delete(anaStatus)
% end

%% Plotting

tHrs = t*tStep*365.25*24;
% NOTE: this includes one extra step so is a tiny bit off, need to
% restructure loops in solve

% Line of Sight
if losPlot == 1
    fLos = figure('Name', 'Line of Sight Window', 'NumberTitle', 'off', 'renderer',...
        'openGL', 'position', [1920/4 1080/4 1920/2 1080/2]);
    
    % Setting LOS data to binary rather than including which checking body
    % is impeding the line of sight
    losBin = NaN(numPoints, numPairs);
    for iter = 1:numPoints
        minPerStep = min(losCheck(iter,:,:));
        losBin(iter,:) = minPerStep;
    end
    % LOS plotting
    legendNames = cell(1,numPairs);
    for iter = 1:numPairs
        fLos(iter) = plot(tHrs,losBin(:,iter)*iter);
        legendNames{iter} = strcat(bodyInfo{bodyPairs(iter,1),2},'-',bodyInfo{bodyPairs(iter,2),2});
        hold on
    end
    axis([0 tHrs(end) 0 numPairs+1])
    set(gca,'ytick',[])
    annotation('textbox',[.04 0 .2 .14],'String','No Window','EdgeColor','none')
    ylabel('<----- Good Windows ----->')
    xlabel('Time [hrs]')
    title('Line-of-Sight Windows')
    legend(legendNames)
end

hold off

% Comm 
if commPlot == 1
    fComm = figure('Name', 'Communications Window', 'NumberTitle', 'off', 'renderer',...
        'openGL', 'position', [1920/4 1080/4 1920/2 1080/2]);
    
    % Setting comm data to binary rather than including which checking body
    % is impeding the comm
    commBin = NaN(numPoints, numPairs);
    for iter = 1:numPoints
        minPerStep = min(commCheck(iter,:,:));
        commBin(iter,:) = minPerStep;
    end
    % comm plotting
    legendNames = cell(1,numPairs);
    for iter = 1:numPairs
        fComm(iter) = plot(tHrs,commBin(:,iter)*iter);
        legendNames{iter} = strcat(bodyInfo{bodyPairs(iter,1),2},'-',bodyInfo{bodyPairs(iter,2),2});
        hold on
    end
    axis([0 tHrs(end) 0 numPairs+1])
    set(gca,'ytick',[])
    annotation('textbox',[.04 0 .2 .14],'String','No Window','EdgeColor','none')
    ylabel('<----- Good Windows ----->')
    xlabel('Time [hrs]')
    title('Communications Windows')
    legend(legendNames)
end





    
    
    
    
    
    
    
