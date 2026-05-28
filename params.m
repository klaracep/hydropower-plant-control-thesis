function p = params()

% Penstock
p.L = 160; % Penstock length [m]
p.D = 4.5; % Penstock diameter
p.D_end = 2.5; % Penstock diameter after narrowing
p.L_vert = 130; % Length of the vertical penstock section [m]
p.L_cone = 15; % Length of the conical penstock section [m]
p.L_hor = p.L- p.L_vert-p.L_cone; % Length of the horizontal penstock section [m]
p.eta = 0.2; % Elbow loss coefficient
p.a = 1400; % Wave speed in water
p.f = 0.02; % Loss coefficient
p.g = 9.81;
p.rho = 1000;
p.HL = 170; % Head at the reservoir side
p.H_net = 160; % Net head reported for HPP Lipno
p.Q_max = 46; % Maximum flow rate for HPP Lipno
p.Ht = 10; % Turbine tailwater head
p.A = pi*p.D^2/4;

% Turbine
p.k = p.Q_max/sqrt(p.H_net); % Turbine constant, computed from Lipno parameters as k = Q_max/(u_max*sqrt(H_net))
p.Pn_el = 60e6; % Nominal electrical power (Lipno)
p.n_nom = 375; % Nominal speed (Lipno)
p.omega_nom = 2*pi*p.n_nom/60;
p.I = 2e6; % Moment of inertia
p.u0 = 0.1; % Initial guess for u0
p.n0 = 50; % Initial speed, selected so that u0 is positive
p.omega0 = 2*pi*p.n0/60;
p.n_ref_final = p.n_nom; % Target speed
p.omega_ref_final = 2*pi*p.n_ref_final/60;

p.Tw = 0.3; % Filter

% Losses (b_loss, c_loss)
% Mechanical and hydraulic losses equal 7% of nominal electrical power
p.b = 0.02 * p.Pn_el / p.omega_nom^2; % Mechanical losses
p.c = 0.05 * p.Pn_el / p.omega_nom^3; % Windage and ventilation losses

% Hydraulic efficiency
% Pn_el = eta_t * P_T_nom - Ploss
% P_T_nom computed from Lipno parameters: P_T_nom = Q_max*rho*g*H_net = 72.2 MW
% Ploss = 0.07*Pn_el
% eta_t = (Pn_el + 0.07*Pn_el)/P_T_nom

Pn_T = p.Q_max*p.rho*p.g*p.H_net;
p.eta_t = (0.07*p.Pn_el + p.Pn_el)/Pn_T;

% Discretization
p.tsim = 200; % Simulation length
p.Nx   = 401; % Penstock nodes
p.dx = p.L/(p.Nx-1);
p.dt = p.dx/p.a;
p.Nt = round(p.tsim/p.dt);

end