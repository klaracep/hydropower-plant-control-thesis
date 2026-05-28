clear; clc; close all;

% Comparison of the full operating cycle for the constant cross-section penstock model
% and the segmented penstock model
% Index classic denotes the constant cross-section penstock model
% Index geom denotes the segmented penstock model


% Parameters
p = params();
cellfun(@(f) assignin('caller',f,p.(f)), fieldnames(p));

x = linspace(0, L, Nx);

[Pi, Ps] = params_pid();

tsim = 200; % Simulation time
Nt = round(tsim/dt);

% Steady state
fun = @(x) steady_equations(x, p, omega0, 0);
sol = fsolve(fun,[1;u0]);

V0 = sol(1);
u0 = sol(2);



% Steady-state profile for the constant cross-section penstock
H_classic = HL - (f/(2*D*g))*V0^2 .* x;

V_classic = ones(Nx,1) * V0;

% Profile for the segmented penstock

i_vert_end = round(L_vert/dx) + 1; % end of S1
i_cone_end = i_vert_end + round(L_cone/dx); % end of S3

D_start = D;

A_start = pi*D_start^2/4;
A_end   = pi*D_end^2/4;

Q0 = A_start * V0;



% Segment S1: constant cross-section

x1 = x(1:i_vert_end);

V1 = Q0/A_start;

H1 = HL - (f/(2*D_start*g))*V1^2 .* (x1 - x1(1));

% Segment S3: conical section

Nc = i_cone_end - i_vert_end + 1;

x3 = x(i_vert_end:i_cone_end);

D3 = linspace(D_start,D_end,Nc);

b = (D_end-D_start)/(x3(end)-x3(1));

H3 = H1(end) ...
    - (2*f*Q0^2/(pi^2*g*b)) ...
    * (1/D_start^4 - 1./D3.^4);

% Segment S4: constant cross-section

x4 = x(i_cone_end:Nx);

V4 = Q0/A_end;

H4 = H3(end) ...
    - (f/(2*D_end*g))*V4^2 .* (x4 - x4(1));

% Head profile

H_geom = [H1(1:end-1) H3(1:end-1) H4];
H_geom = H_geom(:);

% Velocity profile

V_geom = zeros(Nx,1);

% S1

V_geom(1:i_vert_end) = V1;

% S3

for j = 1:Nc

    i = i_vert_end + j - 1;

    A_loc = pi*D3(j)^2/4;

    V_geom(i) = Q0/A_loc;

end

% S4

V_geom(i_cone_end:Nx) = V4;

% Initialization
omega_classic = omega0;
omega_geom    = omega0;

x0_i = Pi.c*omega0 - omega0;
I0_i = u0 / Pi.ri;

Xi_classic = [I0_i x0_i u0 0];
Xi_geom    = [I0_i x0_i u0 0];

u_classic = u0;
u_geom    = u0;

t_sync = 60; % Time of transition to synchronous operation
t_disc = 120; % Time of transition to the disconnection test

P_el_classic = 0;
P_el_geom    = 0;

u_end_isl_classic = u0;
u_end_isl_geom    = u0;

u_end_sync_classic = u0;
u_end_sync_geom    = u0;

% Speed reference ramp
omega_ref_hist = zeros(1,Nt);
omega_ref_hist(1) = omega0;

ramp_time = 10;
ramp_rate = (n_nom - n0) / ramp_time;

% Power reference with three ramps
P_ref = zeros(1, Nt);

tconst_p = 20;
k_sync = floor(t_sync/dt);

w1p = 0.3;
w2p = 0.5;
w3p = 0.9;

ramp_p = 0.3;

k_samp_rampy_1 = floor(w1p / ramp_p / dt);
k_samp_rampy_2 = floor((w2p - w1p) / ramp_p / dt);
k_samp_rampy_3 = floor((w3p - w2p) / ramp_p / dt);

k_cp = floor(tconst_p / dt);

P_ref(k_sync+(1:k_cp)) = w1p * Pn_el;
P_ref(k_sync+(1:k_samp_rampy_1)) = ...
    cumsum(ones(1,k_samp_rampy_1))*dt*ramp_p*Pn_el;

P_ref(k_sync+k_cp+(1:k_cp)) = w2p * Pn_el;
P_ref(k_sync+k_cp+(1:k_samp_rampy_2)) = ...
    w1p*Pn_el + cumsum(ones(1,k_samp_rampy_2))*dt*ramp_p*Pn_el;

P_ref(k_sync+2*k_cp+1:Nt) = w3p * Pn_el;
P_ref(k_sync+2*k_cp+(1:k_samp_rampy_3)) = ...
    w2p*Pn_el + cumsum(ones(1,k_samp_rampy_3))*dt*ramp_p*Pn_el;

% History
% t_hist = zeros(1,Nt);
t_hist = (0:Nt-1)*dt;

omega_classic_hist = zeros(1,Nt);
omega_geom_hist = zeros(1,Nt);
omega_classic_hist(1) = omega0;
omega_geom_hist(1) = omega0;

u_classic_hist = zeros(1,Nt);
u_geom_hist = zeros(1,Nt);
u_classic_hist(1) = u_classic;
u_geom_hist(1) = u_geom;

P_el_classic_hist = zeros(1,Nt);
P_el_geom_hist = zeros(1,Nt);

% Storage of hydraulic head
save_every = 5;
Kmax = floor(Nt / save_every);

H_classic_hist = zeros(Nx, Kmax);
H_geom_hist    = zeros(Nx, Kmax);
t_hist_anim    = zeros(1, Kmax);

k_anim = 0;

% Simulation loop
for n = 2:Nt

    t = t_hist(n);

    n_ref = min(n0 + ramp_rate*t, n_ref_final);
    omega_ref = 2*pi*n_ref/60;
    omega_ref_hist(n) = omega_ref;

    if t < t_sync
        % Island
        [u_classic_cmd, Xi_classic] = PID_parallel(omega_ref, omega_classic, dt, Xi_classic, Pi);
        u_classic = u_classic + dt*(u_classic_cmd - u_classic)/Tw;

        [u_geom_cmd, Xi_geom] = PID_parallel(omega_ref, omega_geom, dt, Xi_geom, Pi);
        u_geom = u_geom + dt*(u_geom_cmd - u_geom)/Tw;

        [H_classic, V_classic, Qc, Hb_c] = MOC_step(H_classic, V_classic, u_classic, p);
        [H_geom, V_geom, Qg, Hb_g] = MOC_step_geom(H_geom, V_geom, u_geom, p);

        [omega_classic, ~] = mechanics_step_island(omega_classic, Qc, Hb_c, 0, p);
        [omega_geom, ~]    = mechanics_step_island(omega_geom, Qg, Hb_g, 0, p);

        u_end_isl_classic = u_classic;
        u_end_isl_geom    = u_geom;

    elseif t < t_disc
        % Sync
        Pref_pu = P_ref(n)/Pn_el;

        Pel_pu_c = P_el_classic/Pn_el;
        Pel_pu_g = P_el_geom/Pn_el;

        if abs(t - t_sync) < dt

           % classic
            w0 = Pref_pu;
            y0 = Pel_pu_c;

            x0 = Ps.c*w0 - y0;
            P0 = Ps.r0*(Ps.b*w0 - y0);
            D0 = 0;

            I0 = (u_end_isl_classic - P0 - Ps.rd*D0)/Ps.ri;
            Xs_classic = [I0 x0 u_end_isl_classic 0];

            % geom
            w0 = Pref_pu;
            y0 = Pel_pu_g;

            x0 = Ps.c*w0 - y0;
            P0 = Ps.r0*(Ps.b*w0 - y0);
            D0 = 0;

            I0 = (u_end_isl_geom - P0 - Ps.rd*D0)/Ps.ri;
            Xs_geom = [I0 x0 u_end_isl_geom 0];
        end


        [u_classic_cmd, Xs_classic] = PID_parallel(Pref_pu, Pel_pu_c, dt, Xs_classic, Ps);
        [u_geom_cmd,    Xs_geom]    = PID_parallel(Pref_pu, Pel_pu_g, dt, Xs_geom, Ps);

        u_classic = u_classic + dt*(u_classic_cmd - u_classic)/Tw;
        u_geom    = u_geom    + dt*(u_geom_cmd    - u_geom)/Tw;

        [H_classic, V_classic, Qc, Hb_c] = MOC_step(H_classic, V_classic, u_classic, p);
        [H_geom, V_geom, Qg, Hb_g] = MOC_step_geom(H_geom, V_geom, u_geom, p);

        [P_el_classic, ~] = mechanics_step_sync(Qc, Hb_c, p);
        [P_el_geom, ~]    = mechanics_step_sync(Qg, Hb_g, p);

        omega_classic = omega_nom;
        omega_geom    = omega_nom;

        u_end_sync_classic = u_classic;
        u_end_sync_geom    = u_geom;

    else
        % Disconnection

        if abs(t - t_disc) < dt
            w0 = omega_nom;

            % classic
            y0 = omega_classic;
            x0 = Pi.c*w0 - y0;
            P0 = Pi.r0*(Pi.b*w0 - y0);
            I0 = (u_end_sync_classic - P0)/Pi.ri;
            Xi_classic = [I0 x0 u_end_sync_classic 0];

            % geom
            y0 = omega_geom;
            x0 = Pi.c*w0 - y0;
            P0 = Pi.r0*(Pi.b*w0 - y0);
            I0 = (u_end_sync_geom - P0)/Pi.ri;
            Xi_geom = [I0 x0 u_end_sync_geom 0];
        end

        [u_classic_cmd, Xi_classic] = PID_parallel(omega_ref, omega_classic, dt, Xi_classic, Pi);
        [u_geom_cmd,    Xi_geom]    = PID_parallel(omega_ref, omega_geom, dt, Xi_geom, Pi);

        u_classic = u_classic + dt*(u_classic_cmd - u_classic)/Tw;
        u_geom    = u_geom    + dt*(u_geom_cmd    - u_geom)/Tw;

        [H_classic, V_classic, Qc, Hb_c] = MOC_step(H_classic, V_classic, u_classic, p);
        [H_geom, V_geom, Qg, Hb_g] = MOC_step_geom(H_geom, V_geom, u_geom, p);

        [omega_classic, ~] = mechanics_step_island(omega_classic, Qc, Hb_c, 0, p);
        [omega_geom, ~]    = mechanics_step_island(omega_geom, Qg, Hb_g, 0, p);

        P_el_classic = 0;
        P_el_geom    = 0;
    end

    % Save data to history arrays
    omega_classic_hist(n) = omega_classic;
    omega_geom_hist(n) = omega_geom;
    u_classic_hist(n) = u_classic;
    u_geom_hist(n) = u_geom;
    P_el_classic_hist(n) = P_el_classic;
    P_el_geom_hist(n) = P_el_geom;

    if mod(n, save_every) == 0
        k_anim = k_anim + 1;

        H_classic_hist(:,k_anim) = H_classic;
        H_geom_hist(:,k_anim)    = H_geom;

        t_hist_anim(k_anim) = t;
    end
end

% Plot

% Full cycle plot
figure
set(gcf,'Units','centimeters')
set(gcf,'Position',[0 0 12.5 18])

ax1 = subplot(3,1,1);
ax1.FontSize = 10;
ax1.TickLabelInterpreter = 'latex';

plot(t_hist, omega_classic_hist*60/(2*pi),'LineWidth',1.5)
hold on
plot(t_hist, omega_geom_hist*60/(2*pi),'LineWidth',1)
plot(t_hist, omega_ref_hist*60/(2*pi),'--','Color',[0.5 0.5 0.5],'LineWidth',1)

ylabel('$n$ [rpm]','Interpreter','latex','FontSize',10)

lgd = legend({'$n$','$n_{g}$','$n_{ref}$'}, ...
    'Location','southeast', ...
    'FontSize',10, ...
    'Interpreter','latex');

lgd.ItemTokenSize = [10 6];
title('HPP full cycle - penstock comparison','Interpreter','latex','FontSize',10)

grid on

ax2 = subplot(3,1,2);
ax2.FontSize = 10;
ax2.TickLabelInterpreter = 'latex';

plot(t_hist, u_classic_hist,'LineWidth',1.5)
hold on
plot(t_hist, u_geom_hist,'LineWidth',1)

ylabel('$u$ [-]','Interpreter','latex','FontSize',10)

lgd = legend({'$u$','$u_{g}$'}, ...
    'FontSize',10, ...
    'Interpreter','latex');

lgd.ItemTokenSize = [10 6];

grid on

ax3 = subplot(3,1,3);
ax3.FontSize = 10;
ax3.TickLabelInterpreter = 'latex';

plot(t_hist, P_el_classic_hist/1e6,'LineWidth',1.5)
hold on
plot(t_hist, P_el_geom_hist/1e6,'LineWidth',1)

mask = (t_hist >= t_sync) & (t_hist < t_disc);

plot(t_hist(mask), P_ref(mask)/1e6, '--', ...
    'Color',[0.5 0.5 0.5], ...
    'LineWidth',1)

ylabel('$P_{el}$ [MW]','Interpreter','latex','FontSize',10)
xlabel('$t$ [s]','Interpreter','latex','FontSize',10)

lgd = legend({'$P_{el}$','$P_{el,g}$','$P_{ref}$'}, ...
    'Location','northeast', ...
    'FontSize',10, ...
    'Interpreter','latex');

lgd.ItemTokenSize = [10 6];

grid on

%exportgraphics(gcf,'geom_vs_classic_full_cycle.pdf','ContentType','vector');



% Plots of hydraulic head time response in synchronous operation

idx_mid = round(Nx/2); % Middle of the penstock
idx_end = Nx; % End of the penstock

% Mask for synchronous operation
id_sync = (t_hist_anim >= t_sync) & (t_hist_anim < t_disc);

t_plot = t_hist_anim(id_sync) - t_sync;

% Plot of hydraulic head time response in the middle of the penstock
figure
set(gcf,'Units','centimeters')
set(gcf,'Position',[0 0 12.5 8.5])

ax1 = gca;
ax1.FontSize = 10;
ax1.TickLabelInterpreter = 'latex';

hold on

plot(t_plot, H_classic_hist(idx_mid, id_sync), ...
    'LineWidth',1.5)

plot(t_plot, H_geom_hist(idx_mid, id_sync), ...
    'LineWidth',1.5)

xlabel('$t$ [s]','Interpreter','latex','FontSize',10)
ylabel('$H$ [m]','Interpreter','latex','FontSize',10)

lgd = legend({'$H$','$H_g$'}, ...
    'Interpreter','latex', ...
    'FontSize',10, ...
    'Location','southeast');

lgd.ItemTokenSize = [10 6];
title('Hydraulic head evolution in synchronous mode - x = 0.5L','Interpreter','latex','FontSize',10)

grid on

%exportgraphics(gcf,'H_mid_sync.pdf','ContentType','vector')

% Plot of hydraulic head time response at the end of the penstock
figure
set(gcf,'Units','centimeters')
set(gcf,'Position',[0 0 12.5 8.5])

ax2 = gca;
ax2.FontSize = 10;
ax2.TickLabelInterpreter = 'latex';

hold on

plot(t_plot, H_classic_hist(idx_end, id_sync), ...
    'LineWidth',1.5)

plot(t_plot, H_geom_hist(idx_end, id_sync), ...
    'LineWidth',1.5)

xlabel('$t$ [s]','Interpreter','latex','FontSize',10)
ylabel('$H$ [m]','Interpreter','latex','FontSize',10)

lgd = legend({'$H_b$','$H_{b,g}$'}, ...
    'Interpreter','latex', ...
    'FontSize',10, ...
    'Location','southeast');

lgd.ItemTokenSize = [10 6];
title('Hydraulic head evolution in synchronous mode - x = L','Interpreter','latex','FontSize',10)

grid on

%exportgraphics(gcf,'H_end_sync.pdf','ContentType','vector')