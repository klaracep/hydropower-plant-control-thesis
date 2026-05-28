clc; clear; close all;
% Comparison of the linear and nonlinear model response to a disturbance in island operation
% Step changes in load torque (percentages of the nominal load) - 5%, 10%, 20%,
% 30%
% Operating point computed for 375 rpm and load corresponding to 90% of the
% nominal electrical power
% Default simulation setting is for closed-loop operation with PID


% Parameters
p = params();
cellfun(@(f) assignin('caller',f,p.(f)), fieldnames(p));

tsim = 30; % Simulation time

Nt = round(tsim/dt);
t  = (0:Nt-1)*dt;

P_el_p = 0.9*Pn_el;
Mz_nom = P_el_p/ omega_nom; % Nominal Mz
omega_ref = omega_nom;

% Operating point

fun = @(x) steady_equations(x, p, omega_ref, Mz_nom);
x0  = [1; 0.3];
sol = fsolve(fun, x0);

V_p = sol(1);
u_p = sol(2);


% Hydraulic initialization
x_space = linspace(0,L,Nx)';


% Linearization

omega_p = omega_ref;

Q_p  = A*V_p;
hf_p = f*L/(2*D*g)*V_p^2;
Hb_p = HL - hf_p;

fprintf('Working point:\n');
fprintf('V_p = %.3f m/s\n', V_p);
fprintf('u_p = %.3f\n', u_p);
fprintf('Q_p = %.3f m3/s\n', Q_p);
fprintf('Hb_p = %.3f m\n', Hb_p);

% Constants of the linear model
B  = (g*A)/a;
Kh = (k*u_p)/(2*sqrt(Hb_p-Ht));
Ku = k*sqrt(Hb_p-Ht);

CQ = (rho*g*eta_t)*(Hb_p-Ht)/omega_p;
CH = (rho*g*eta_t)*Q_p/omega_p;
Cw = (rho*g*eta_t)*Q_p*(Hb_p-Ht)/omega_p^2;

Dw = b + c*omega_p; % Linearization of losses
T  = I/(Cw + Dw);
tau= (Nx-1)*dt; % number_of_segments*dt

s = tf('s');
delay = exp(-2*tau*s);

num_h = (B*CQ - CH) + (B*CQ + CH)*delay;
den_h = (Kh + B) + (B - Kh)*delay;

G_h = Ku * num_h / den_h;
G_f = 1/(Tw*s + 1);
G_mech = 1 / ((Cw + Dw)*(T*s + 1));

Gu = G_h * G_mech * G_f; % u to omega
Gm = G_mech; % Mz to omega (disturbance)

Gud = c2d(Gu, dt, 'zoh');
Gmd = c2d(Gm, dt, 'zoh');

Gd_ss = ss([Gud -Gmd]); % omega = Gu*u - Gm*Mz

A_lin = Gd_ss.A;
B_lin = Gd_ss.B;
C_lin = Gd_ss.C;
D_lin = Gd_ss.D;

nx = size(A_lin,1);

% PID initialization
[Pi, ~] = params_pid();
P = Pi;


t_step = 10; % Disturbance time

% Loop
steps = [0.95 0.9 0.8 0.7];    % percentages of Mz_nom

figure
set(gcf,'Units','centimeters')
set(gcf,'Position',[0 0 12.5 18])

for i = 1:length(steps)

    % Hydraulic initialization
    % Steady-state profile
    % From the continuity and momentum equations, where d/dt is zero
    % Continuity: dV/dx = 0 -> Vp = constant
    % Momentum: g*(dH/dx) + f/(2*D)*Vp^2 = 0 -> integration -> profile

    H = HL - (f/(2*D*g))*V_p^2 .* x_space;
    V = ones(Nx,1)*V_p;

    % Mechanical initialization
    omega_nl  = zeros(1,Nt);
    omega_lin = zeros(1,Nt);
    omega_nl(1)  = omega_p;
    omega_lin(1) = omega_p;

    % Linear model initialization
    x_lin = zeros(nx,1);

    % PID initialization
    w0 = omega_ref;
    y0 = omega_p;

    x0 = Pi.c*w0 - y0;

    P0 = Pi.r0*(Pi.b*w0 - y0);
    D0 = 0;

    I0 = (u_p - P0 - Pi.rd*D0)/Pi.ri;

    Xi_lin = [I0 x0 u_p 0];
    Xi_nl  = [I0 x0 u_p 0];
    
    u_nl = u_p;

    % Load torque
    Mz = Mz_nom * ones(1,Nt);
    Mz(t > t_step) = steps(i) * Mz_nom;

    % Loop
    for k = 2:Nt

        % Nonlinear model
        [u_nl_cmd, Xi_nl] = PID_parallel(omega_ref, omega_nl(k-1), dt, Xi_nl, P);
        
        % Servomechanism
        du = (u_nl_cmd - u_nl)/Tw;
        u_nl = u_nl_cmd + dt*du;

        %u_nl = u_p; % without control % UNCOMMENT for simulation without PID

        [H, V, Q_T, H_b] = MOC_step(H, V, u_nl, p);

        [omega_nl(k), ~] = mechanics_step_island(omega_nl(k-1), Q_T, H_b, Mz(k), p);

        % Linear model
        [u_lin_abs, Xi_lin] = PID_parallel(omega_ref, omega_lin(k-1), dt, Xi_lin, P);

        du  = u_lin_abs - u_p; % deviation
        dMz = Mz(k) - Mz_nom;  % deviation

        %du  = 0; % without control % UNCOMMENT for simulation without PID
        %dMz = Mz(k) - Mz_nom; % UNCOMMENT for simulation without PID

        u_vec = [du; dMz]; % inputs

        x_lin = A_lin*x_lin + B_lin*u_vec;
        domega = C_lin*x_lin + D_lin*u_vec;

        omega_lin(k) = omega_p + domega; % absolute value
    end

    % Error
    err = omega_nl - omega_lin;
    rms_err = sqrt(mean(err.^2));

    fprintf('Drop in Mz %.0f%%  RMS error = %.6e\n',(1-steps(i))*100, rms_err);
    
    ax = subplot(4,1,i);

    ax.FontSize = 10;
    ax.TickLabelInterpreter = 'latex';
    plot(t, omega_nl*60/(2*pi),'LineWidth',1.5)
    hold on
    plot(t, omega_lin*60/(2*pi),'LineWidth',1)
    title(sprintf('Drop in $M_z$ %.0f\\%%', (1-steps(i))*100), 'Interpreter','latex', ...
    'FontSize',11)
    xlabel('t [s]','Interpreter','latex', ...
    'FontSize',10)
    ylabel('n [rpm]','Interpreter','latex', ...
    'FontSize',10)
    %ylim([375 inf]) % UNCOMMENT for simulation without PID

    text(0.02, 0.85, sprintf('RMS = %.2e', rms_err), ...
    'Units','normalized', ...
    'Interpreter','latex', ...
    'FontSize',10, ...          
    'BackgroundColor','white', ...
    'EdgeColor','black', ...
    'Margin',1)    

    lgd = legend({'Nonlinear','Linear'}, ...
             'Interpreter','latex', ...
             'FontSize',10,'Location','southeast');

    lgd.ItemTokenSize = [10 6]; 

    grid on
end

set(findall(gcf,'-property','Interpreter'),'Interpreter','latex')

%exportgraphics(gcf, 'lin_vs_nlin_island.pdf','ContentType','vector');