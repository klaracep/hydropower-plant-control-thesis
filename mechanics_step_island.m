function [omega_new, P_T] = mechanics_step_island(omega, Q_T, H_b, Mz, p)

% Turbine mechanics for island operation
rho = p.rho;
g   = p.g;
Ht  = p.Ht;
eta_t = p.eta_t;

I  = p.I;
dt = p.dt;

b = p.b;   % Mechanical loss - b_loss
c = p.c;   % Hydraulic loss - c_loss

% Mechanical power of the turbine
P_T = eta_t *rho * g * Q_T * (H_b - Ht);

% domega = 1/I * (f(omega) - Mz)
% omega(k+1) = omega(k)+ dt/(I*(1-dt*J)/(2*I))*(f(omega) - Mz)
% f(omega)
f = (1/max(omega, 1e-3))*P_T - b*omega - c * omega^2;

% Jacobian
J = -b - 2*c*omega- (1/max(omega, 1e-3)^2)*P_T;

den = I*(1-((dt*J)/(2*I)));

% Tustin
omega_new = omega + dt/den*(f - Mz);

end