function [prob_data, adapt] = init_data()
% function [prob_data, adapt] = init_data()
% In this function all problem data and adaptivity parameters are initialized.
% This function is called from afem.m

% example = 'big-square-Dirichlet'
% example = 'square-Dirichlet'
% example = 'square-mixed'
example = 'L-shape'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  problem  data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% we declare a "struct" for storing the problem data
% diffusion coefficient (a) of the equation
prob_data.a = 1;
% convection coefficient (b) of the equation (row vector)
prob_data.b = [0.0  0.0];
% reaction coefficient (c) of the equation
prob_data.c = 0.0;
% folder or directory where the domain mesh is described
prob_data.domain = 'square_all_dirichlet';
prob_data.initial_global_refinements = 0;
prob_data.f = @(x) 0;
prob_data.gD = @(x) 0;
prob_data.gN = @(x) 0;
% (these default field values can be changed inside each example, below)

switch example
case 'big-square-Dirichlet'
  prob_data.initial_global_refinements = 1;
  prob_data.domain = 'big_square_all_dirichlet';
  prob_data.u_exact = @(x) exp(-10*norm(x)^2);
  prob_data.grad_u_exact = @(x) -20*prob_data.u_exact(x).*x;
  prob_data.f = @(x) 20*exp(-10*norm(x)^2)*(2-20*norm(x)^2);
  prob_data.gD = @(x) prob_data.u_exact(x);
case 'square-Dirichlet'
  prob_data.initial_global_refinements = 1;
  prob_data.domain = 'square_all_dirichlet';
  prob_data.u_exact = @(x) exp(-10*norm(x)^2);
  prob_data.grad_u_exact = @(x) -20*prob_data.u_exact(x).*x;
  prob_data.f = @(x) 20*exp(-10*norm(x)^2)*(2-20*norm(x)^2);
  prob_data.gD = @(x) prob_data.u_exact(x);
case 'square-mixed'
  domain = 'square_mixed';
  prob_data.initial_global_refinements = 1;
  prob_data.domain = 'square_mixed';
  prob_data.u_exact = @(x) 0;
  prob_data.grad_u_exact = @(x) 0*x;
  prob_data.f = @(x) 1;
  prob_data.gD = @(x) 0;
case 'L-shape'
  prob_data.domain = 'L_shape_dirichlet';
  prob_data.initial_global_refinements = 0;
  theta = @(x) atan2(x(2),x(1))+(2*pi*(x(2)<0));
  prob_data.u_exact = @(x) norm(x)^(2/3)*sin(2/3*theta(x));
  prob_data.grad_u_exact = @(x) 2/3*norm(x)^(-1/3)*[-sin(theta(x)/3); cos(theta(x)/3)];
  prob_data.f = @(x) 0;
  prob_data.gD = @(x) prob_data.u_exact(x);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  data for a posteriori estimators and adaptive strategy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
% Doerfler: Doerfler 'bulk' strategy
adapt.strategy = 'Doerfler';

% n_refine, number of refinements of each marked element
adapt.n_refine = 2;

% parameters of the different marking strategies
% MS: Maximum strategy
adapt.MS_gamma = 0.5;

% Doerfler: Doerfler 'bulk' strategy
adapt.Doerfler_theta = 0.5;
adapt.Doerfler_nu = 0.1;

