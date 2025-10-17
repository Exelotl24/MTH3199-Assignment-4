%this function computes the orbit of a planet about a sun
%the sun is assumed to located at the origin (0,0)
%and motion is restricted to the x-y plane
%INPUTS:
%t_list: a list of times to compute the position & velocity of the planet
%V0: initial condition. V0 is a column vector consisting
% of the initial position and velocity of the planet:
% V0 = [x(0); y(0); dx/dt(0); dy/dt(0)]
%orbit_params: a struct describing the system parameters
% orbit_params.m_sun: mass of the sun
% orbit_params.m_planet: mass of the planet
% orbit_params.G: gravitational constant
% Force = -m_planet*m_sun*G/rˆ2
%OUTPUTS:
%V_list: the state of the planet at each time in t_list
% if t_list is a list, then V_list is a Nx4 MATRIX
% where each ROW has the form [x_i,y_i,dx/dt_i,dy/dt_i]
% corresponding to the values at the ith time
% if t_list is a SCALAR (i.e. t_list = t),
% then V_list is a COLUMN VECTOR of the form:
% [x(t); y(t); dx/dt(t); dy/dt(t)]
%NOTES:
%This function needs all the other functions in this file to run
%I HIGHLY RECOMMEND JUST SAVING THIS FUNCTION IN ITS OWN FILE
%DON'T COPY AND PASTE INTO YOUR CODE! IT'S NOT WORTH IT!
%
%USAGE EXAMPLE:
%At the bottom of this file is a function usage_example()
%which shows how to use compute_planetary_motion(...)
%You can start from there and then tweak it.
function V_list = compute_planetary_motion(t_list,V0,orbit_params)
%Orion's special code here
end