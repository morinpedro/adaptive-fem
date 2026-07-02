% file init_data.m
% in this file all problem data and adaptivity parameter are initialized
% this file is called from afem.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  problem  data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% folder or directory where the domain mesh is described
domain = 'L_shape_dirichlet';
% initial global refinements
global_refinements = 1;

% we declare a "struct" for storing the problem data
global prob_data

% diffusion coefficient (a) of the equation
prob_data.a = 1;
% convection coefficient (b) of the equation (row vector)
prob_data.b = [0.0  0.0];
% reaction coefficient (c) of the equation
prob_data.c = 0.0;

% right-hand side function f
%prob_data.f = @(x) sin(pi*x(1))*sin(pi*x(2))/4/pi;
%prob_data.f = @(x) 2*(x(1)>0.5);
%prob_data.f = @(x) exp(-100*sum((x-1/2).*(x-1/2)));
%prob_data.f = @(x) 20*exp(-10*norm(x)^2)*(2-20*norm(x)^2);
prob_data.f = @(x) 0;

% Dirichlet data, function g_D
prob_data.gD = @(x) u_ex3(x);

% Neumann data, function g_N
prob_data.gN = @(x) 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  data for a posteriori estimators and adaptive strategy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% We declare a data structure for storing all 
% the adaptive strategy parameters
global adapt

% weight in front of interior residual
adapt.C(1) = 1.0;
% weight in front of jump residual
adapt.C(2) = 1.0;

% tolerance for the adaptive strategy
adapt.tolerance = 1e-7;

% maximum number of iterations of adaptive strategy
adapt.max_iterations = 10;

% marking_strategy, possible options are
% GR: global (uniform) refinement,  
% MS: maximum strategy,  
% GERS: guaranteed error reduction strategy (D\"orfler's)
% ES: equidistribution strategy,           (not implemented yet) 
% MES: modified equidistribution strategy, (not implemented yet)
adapt.strategy = 'MS';

% n_refine, number of refinements of each marked element
adapt.n_refine = 2;

% parameters of the different marking strategies
% MS: Maximum strategy
adapt.MS_gamma = 0.5;

% GERS: guaranteed error reduction strategy (D\"orfler's)
adapt.GERS_theta_star = 0.8;
adapt.GERS_nu = 0.1;
