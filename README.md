# OrbitalSimulator
Simulator for multiple bodies in an orbit or trajectory.  Produces communication and body visible (e.g. sun seeing) windows.


Uses gravitational forces between bodies and initial ephemeris data to calcualte the orbits and trajectories of multiple bodies in a system.


Run_system calls solve_system and plot_system and produces:

Raw data of body positions and a corresponding time vector

Plot of orbits/trajectories and animated visualization of bodies over chosen time span


Analyze_system requires data produced by run_system and produces:

Plots and raw data of visual contact betwween bodies (designed to give sun seeing windows)

Plots and raw data of communications window between two bodies (visual contact and within comm range)

Future versions will take into account travel time for comm signals, and will introduce spacecraft dynamics
so that antenna pointing requirements/effeciency can be taken into account.


Use this link to generate ephemeris data: 

https://ssd.jpl.nasa.gov/horizons.cgi

Future versions will include automatic import capability.


Direct questions or suggestions to:

Cy Haukdal

B.S., M.S. Aerospace Engineering

Haukdal@vt.edu


Copyright (c) 2019, Cy Haukdal

All rights reserved.
