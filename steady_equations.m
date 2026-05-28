function F = steady_equations(x, p, omega0, Mz)
% Computation of steady-state equilibrium
% Balance equations
% HL-Ht (total head) = hf (penstock losses) + ht (turbine head)
% Torque balance at steady state: 0 = M_T - M_loss - Mz

V = x(1);
u = x(2);

Q = p.A*V;

% Hydraulics
hf = p.f*p.L/(2*p.D*p.g)*V^2;
ht = (Q/(p.k*u))^2;

F1 = hf + ht - (p.HL - p.Ht);

% Head at the turbine inlet
Hb = p.HL - hf;

% Turbine power
P_T = p.rho*p.g*Q*(Hb - p.Ht);
M_T = p.eta_t * P_T / omega0;

% Mechanical losses b_loss, c_loss
M_loss = p.b*omega0 + p.c*omega0^2;

F2 = M_T - M_loss - Mz;

F = [F1; F2];

end