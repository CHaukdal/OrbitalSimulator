function [t, y] = solveOrbSim(bodyInfo, selBodies, tVec, errVec)
%% function [t, y] = solveOrbSim(bodyInfo, selBodies, tVec, errVec)    
% Version 1.1.1
% 
% Takes body information and other settings and produces and saves raw 
% position and velocity data over time step vector
%
% Called by script: analyze_system.m        
% 
% INPUTS:
% bodyInfo          Cell array of user inputted info about bodies to simulate
% selBodies         Which bodies considered in gravity calculations                 
% tVec              Vector containing timing info
% errVec            Vector containing ode solver error info    
% 
% OUTPUTS:
% t                 Time vector, length of (simTime/tStep)+1       
% y                 Position x,y,z (y(:,1*b), y(:,2*b), y(:,3*b)) and velocity ux,uy,uz
%                   (y(:,4*b), y(:,5*b), y(:,6*b)) of body number b
%                   corresponding to time step of vector t
% nBodies           Number of bodies analyzed (for plotting)  
%                      
%  
% Cy Haukdal
% Haukdal@vt.edu
% 
% Copyright (c) 2019, Cy Haukdal
% All rights reserved.
 

%% Global Variable Definitions

global numBodies;
global G;
global mu;
global incBodies;
global totTime;
global odeStatus;
G = 6.674e-11;


%% Sorting Input Data

% Unpacking tVec
iniTime = tVec(1);
tStep = tVec(2);
simTime = tVec(3);
totTime = (simTime-iniTime)*365.35*24*60*60;

% Unpacking errVec
relErrTol = errVec(1);
posCoeff = errVec(2);
velCoeff = errVec(3);

% Indexing chosen bodies to analyze
infoSize = size(bodyInfo);
numBodies = infoSize(1);

incBodies = selBodies; % Globals are annoying but necessary

% Time in seconds rather than fractional steps
realTime = (0:tStep:simTime)*365.25*24*60*60;

% Initializing vectors for speed
mu = zeros(1,numBodies);

bodies = cell(1,numBodies);

semiMajor = zeros(1,numBodies);

x = zeros(1,numBodies);
y = zeros(1,numBodies);
z = zeros(1,numBodies);

ux = zeros(1,numBodies);
uy = zeros(1,numBodies);
uz = zeros(1,numBodies);

% Parsing bodyInfo
for iter = 1:numBodies
    
    % Gravitational parameters for each body [m^3 kg^-1 s^-2]
    mu(iter) = bodyInfo{iter,1}*G;
    
    % Names for labeling
    bodies{iter} = bodyInfo{iter,2};
    
    % Average distange from solar system barycenter [m]
    semiMajor(iter) = bodyInfo{iter,3}*1e3;
    
    % Positions [m]
    x(iter) = bodyInfo{iter,4}*1e3;
    y(iter) = bodyInfo{iter,5}*1e3;
    z(iter) = bodyInfo{iter,6}*1e3;
    
    % Velocities [m/s]
    ux(iter) = bodyInfo{iter,7}*1e3;
    uy(iter) = bodyInfo{iter,8}*1e3;
    uz(iter) = bodyInfo{iter,9}*1e3;   
    
end


%% Solving System

% Ephemeris vector and absolute tolerance vector in form to feed ode113 
% preallocated for speed
ephVec = zeros(1,6*numBodies);
absTol = zeros(1,6*numBodies);
for iter=1:numBodies
    ephVec((iter-1)*6+1:(iter-1)*6+6) = [x(iter) y(iter) z(iter) ux(iter) uy(iter) uz(iter)];
    % Normalized RMS error will divide by zero for bodies with a zero initial
    % position or initial velocity
    if iter == 1
        absTol((iter-1)*6+1:(iter-1)*6+6) = [posCoeff posCoeff posCoeff velCoeff velCoeff velCoeff];
    else
        absErrTolPos = sqrt(x(iter)^2+y(iter)^2+z(iter)^2)/semiMajor(iter)*posCoeff;
        absErrTolVel = sqrt(ux(iter)^2+uy(iter)^2+uz(iter)^2)*velCoeff;
        absTol((iter-1)*6+1:(iter-1)*6+6) = [absErrTolPos absErrTolPos absErrTolPos absErrTolVel absErrTolVel absErrTolVel];
    end
end

% Determining options for ode solver and starting watibar
odeTime = 0;
odeStatus = waitbar(odeTime, 'ODE Solver Progress');
options=odeset('RelTol', relErrTol, 'AbsTol', absTol, 'OutputFcn', @statusBar); 

% Solving using ODE113
% ODE113 chosen for its sutiability for orbital mechanics
[t,y]=ode113(@odeIni,realTime,ephVec,options);

if ishandle(odeStatus)
    delete(odeStatus)
end

% Funciton end
end


%% Local Functions

function Derivatives = odeIni(~,y)
% Creating time derivatives for ODE113 input
% This part actually does the physics 

global numBodies;
global mu;
global incBodies;

% Preallocating [dx,dy,dz,dux,duy,duz] for speed
% Column vector for ode113 input
Derivatives = zeros(6*(numBodies),1);

for iter = 1:numBodies
    
    % Only using bodies that have been chosen to exert gravitational influence 
    ind = incBodies;

    % Preventing body from calculating against itself and dividing by zero
    if ind(iter) == iter
        ind(iter) = [];
    end

    % Determining change in positon by looking at velocity 
    % x
    Derivatives((iter-1)*6+1) = y((iter-1)*6+4);
    % y
    Derivatives((iter-1)*6+2) = y((iter-1)*6+5);
    % z
    Derivatives((iter-1)*6+3) = y((iter-1)*6+6);
    
    % Denominator for inverse square law for all bodies at the same time,
    % left squared to decrease number of operations
    denom = (y((iter-1)*6+1)-y((ind-1)*6+1)).^2+...%x
        (y((iter-1)*6+2)-y((ind-1)*6+2)).^2+... %y
        (y((iter-1)*6+3)-y((ind-1)*6+3)).^2; %z
    
    % Inverse square law for each body to determine acceleration 
    % Initial Velocity - acceleration - calibration 
    % x
    Derivatives((iter-1)*6+4) = Derivatives((iter-1)*6+4) - sum(mu(ind)'./denom.*(y((iter-1)*6+1) - y((ind-1)*6+1))./sqrt(denom));
    % y
    Derivatives((iter-1)*6+5) = Derivatives((iter-1)*6+5) - sum(mu(ind)'./denom.*(y((iter-1)*6+2) - y((ind-1)*6+2))./sqrt(denom));
    % z
    Derivatives((iter-1)*6+6) = Derivatives((iter-1)*6+6) - sum(mu(ind)'./denom.*(y((iter-1)*6+3) - y((ind-1)*6+3))./sqrt(denom));  
end

% Function end
end

function imOnlyPuttingThisVariableHereSomMatalbDoesntGetMad = statusBar(t,~,~)
% Updates waitbar for ode113 progress
global totTime;
global odeStatus;
    
if size(t) == 1 & ishandle(odeStatus) % I'd love to use a short circuit && matalb but you give me an error if i do ¯\_(?)_/¯
    odeTime = t;
    waitbar(odeTime/totTime,odeStatus)
end
imOnlyPuttingThisVariableHereSomMatalbDoesntGetMad = 0;

% Function end
end
