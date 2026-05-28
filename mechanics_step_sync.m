function [P_el, P_T] = mechanics_step_sync(Q_T, H_b, p)
% Turbine mechanics for synchronous operation

rho = p.rho;
g   = p.g;
Ht  = p.Ht;
eta_t = p.eta_t;

omega_s = p.omega_nom;
b = p.b; % b_loss
c = p.c; % c_loss

% Mechanical power of the turbine
P_T = rho * g * Q_T * (H_b - Ht);

% Losses
M_loss = b * omega_s + c * omega_s^2;
P_loss = M_loss*omega_s;

% Electrical power supplied to the grid
P_el = eta_t*P_T - P_loss;

end


