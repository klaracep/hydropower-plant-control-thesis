clear; clc; close all;
% Full power plant cycle - nonlinear model

% Parameters
p = params();
cellfun(@(f) assignin('caller',f,p.(f)), fieldnames(p));

x = linspace(0, L, Nx);

% PID parameters
% Pi - island operation
% Ps - synchronous operation
[Pi, Ps] = params_pid();
w0 = omega0;
y0 = omega0;

x0_i = Pi.c*w0 - y0;
I0_i = u0 / Pi.ri;
fI0 = 0;

Xi = [I0_i x0_i u0 fI0];

x0_s = Ps.c*0 - 0;
I0_s = u0 / Ps.ri;
fI0 = 0;

Xs = [I0_s x0_s u0 fI0];


% Initial condition - steady state
fun = @(x) steady_equations(x, p, omega0, 0);

x_guess = [1; u0];

sol = fsolve(fun, x_guess);

V0 = sol(1);
u0 = sol(2);

Q0 = p.A * V0;

fprintf('Steady V0 = %.4f m/s\n', V0);
fprintf('Steady u0 = %.4f\n', u0);
fprintf('Steady Q0 = %.4f m3/s\n', Q0);

u = u0;

% Steady-state profile
% From the continuity and momentum equations, where d/dt is zero
% Continuity: dV/dx = 0 -> V0 = constant
% Momentum: g*(dH/dx) + f/(2*D)*V0^2 = 0 -> integration -> profile
H = HL - (f/(2*D*g))*V0^2 .* x;
V = ones(Nx,1) * V0;

Hb = H(end); % Penstock-turbine boundary
hf = HL - Hb;
hturb = Hb - Ht;

% Check of the initial state
fprintf('Hb (penstock end head) = %.8f m\n', Hb);
fprintf('hf (pipe friction loss)= %.8f m\n', hf);
fprintf('hturb (turb drop)      = %.8f m\n', hturb);
fprintf('HL-Ht                  = %.8f m\n', HL - Ht);
fprintf('hf+hturb               = %.8f m\n', hf + hturb);

hturb_from_Q = (Q0/(k*u0))^2;

fprintf('Hb-Ht           = %.6f m\n', Hb - Ht);
fprintf('(Q0/(p.k*u0))^2  = %.6f m\n', hturb_from_Q);

% Allocation
t_hist = (0:Nt-1)*dt;
omega_hist = zeros(1, Nt);
u_hist = zeros(1,Nt);
omega_ref_hist = zeros(1,Nt);
P_el_hist = zeros(1,Nt);

omega = omega0;
omega_hist(1) = omega;
omega_ref_hist(1) = omega;
u_hist(1)         = u0;


% Animation
save_every = 5;

Kmax = floor(Nt / save_every);

H_hist = zeros(Nx, Kmax);
V_hist = zeros(Nx, Kmax);
t_hist_anim = zeros(1, Kmax);

k_anim = 0;

% Mode switching times
t_sync = 60;
t_disc = 120;
P_el = 0;
u_end_isl = u0;
u_end_sync = u0;
P_el_end_sync = 0;

% Speed reference
ramp_time = 10;   % s
ramp_rate = (n_nom - n0) / ramp_time; 


% Electrical power reference
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

t_island_start = 0; % Time when the plant start-up begins

% Simulation loop
for n = 2:Nt

    t = t_hist(n);
    
    n_ref = min(n0 + ramp_rate*t, n_ref_final);
    omega_ref = 2*pi*n_ref/60;

    if t < t_island_start

        u = u0;

        [H, V, Q_T, H_b] = MOC_step(H, V, u, p);

        Mz = 0;
        [omega, P_T] = mechanics_step_island(omega, Q_T, H_b, Mz, p);

    elseif t >= t_island_start && t < t_sync
        % Island operation
        [u_cmd, Xi] = PID_parallel(omega_ref, omega, dt, Xi, Pi);

        % Servomechanism
        du = (u_cmd - u)/Tw;
        u = u + dt*du;
        
        % Hydraulics from MOC
        [H, V, Q_T, H_b] = MOC_step(H, V, u, p);
        
        % Check
        if t >= t_sync-0.05 && t < t_sync-0.05+dt
        fprintf(['t=%.3f, u=%.6f\n' ...
             'H(0)=%.6f, H(L)=%.6f, dH=%.6f\n' ...
             'H_b=%.6f, Q_T=%.6f\n' ...
             '(Q_T/(k*u))^2=%.6f, H_b-Ht=%.6f\n'], ...
        t, u, H(1), H(end), H(1)-H(end), H_b, Q_T, (Q_T/(p.k*u))^2, H_b-p.Ht);
        end
       % Turbine mechanics
        Mz = 0; % Island operation - load is zero
        [omega, P_T] = mechanics_step_island(omega, Q_T, H_b, Mz, p);

        P_el = 0;
        u_end_isl = u; % Save state for transition to the next operating mode

    elseif t>=t_sync && t<t_disc
        % Synchronous operation
        % Controller initialization

        Pref_pu = P_ref(n)/Pn_el; % Normalization
        Pel_pu  = P_el/Pn_el;

        if  abs(t - t_sync) < dt
            omega_sync = omega;
            % w0 = P_ref(n);
            % y0 = P_el;

            w0 = Pref_pu;
            y0 = Pel_pu;

            x0_s = Ps.c*w0 - y0;

            P0 = Ps.r0*(Ps.b*w0 - y0);
            D0 = 0;
            I0_s = (u_end_isl - P0 - Ps.rd*D0)/Ps.ri;

            Xs = [I0_s x0_s u_end_isl 0];
            
        end
        % Rotational speed remains synchronous
        omega = omega_sync;

        % Controller

        [u_cmd, Xs] = PID_parallel(Pref_pu, Pel_pu, dt, Xs, Ps);

        % Servomechanism
        du = (u_cmd - u)/Tw;
        u = u + dt*du;

        % Hydraulics
        [H, V, Q_T, H_b] = MOC_step(H, V, u, p);

        [P_el, P_T] = mechanics_step_sync(Q_T, H_b, p);
        
        % Save state for transition to the next operating mode
        P_el_end_sync = P_el;
        u_end_sync = u;

    else
        % Disconnection test
        % Controller initialization
        if abs(t -t_disc) < dt
            w0 = omega_nom;
            y0 = omega;

            x0_i = Pi.c*w0 - y0;

            P0 = Pi.r0*(Pi.b*w0 - y0);
            D0 = 0;
            I0_i = (u_end_sync - P0 - Pi.rd*D0)/Pi.ri;

            Xi = [I0_i x0_i u_end_sync 0];
        end
        
        % Speed reference - maintain nominal value
        n_ref = n_nom;
        omega_ref = (2*pi*n_ref)/60;

        [u_cmd, Xi] = PID_parallel(omega_ref, omega, dt, Xi, Pi);

        % Servomechanism
        du = (u_cmd - u)/Tw;
        u = u + dt*du;

        [H, V, Q_T, H_b] = MOC_step(H, V, u, p);
        
        % Power drop to zero
        Mz = 0;
        P_el = Mz/omega_nom;

        [omega_new, P_T] = mechanics_step_island(omega, Q_T, H_b, Mz, p);
        omega = omega_new;
        
    end

    omega_hist(n)     = omega;
    u_hist(n)         = u;
    omega_ref_hist(n) = omega_ref;
    P_el_hist(n)      = P_el;

    if mod(n, save_every) == 0
        k_anim = k_anim + 1;
        H_hist(:,k_anim) = H;
        V_hist(:,k_anim) = V;
        t_hist_anim(k_anim) = t;
    end
end

% Check
omega_ss_sim = mean(omega_hist(end-200:end));
n_ss_sim = omega_ss_sim * 60 / (2*pi);
n_syn = omega_nom * 60 / (2*pi);

fprintf('Simulated steady speed   n = %.3f rpm\n', n_ss_sim);
fprintf('Synchronous speed        n = %.3f rpm\n', n_syn);
fprintf('Speed deviation          = %.3e %%\n', (n_ss_sim - n_syn)/n_syn*100);



% Plots

% Full power plant cycle
figure
set(gcf,'Units','centimeters')
set(gcf,'Position',[0 0 12.5 18])

ax = subplot(3,1,1);
ax.FontSize = 10;
ax.TickLabelInterpreter = 'latex';
plot(t_hist, omega_hist*60/(2*pi), 'LineWidth', 1.5)
hold on
plot(t_hist, omega_ref_hist*60/(2*pi),'--', 'Color',[0.5 0.5 0.5],'LineWidth',1)
ylabel('n [rpm]','Interpreter','latex','FontSize',10)
legend({'$n$','$n_{ref}$'}, 'Interpreter','latex','FontSize',10,'Location', 'southeast')
title('HPP full cycle','Interpreter','latex','FontSize',10)
grid on

ax = subplot(3,1,2);
ax.FontSize = 10;
ax.TickLabelInterpreter = 'latex';
plot(t_hist, u_hist, 'LineWidth', 1.5)
ylabel('u [-]','Interpreter','latex','FontSize',10)
grid on
 
ax = subplot(3,1,3);
ax.FontSize = 10;
ax.TickLabelInterpreter = 'latex';
plot(t_hist, P_el_hist/1e6, 'LineWidth', 1.5)
hold on


mask = (t_hist >= t_sync) & (t_hist < t_disc);
plot(t_hist(mask), P_ref(mask)/1e6, '--', 'Color',[0.5 0.5 0.5],'LineWidth',1)

ylabel('$P_{el}$ [MW]','Interpreter','latex','FontSize',10)
xlabel('t [s]','Interpreter','latex','FontSize',10)
legend({'$P_{el}$','$P_{ref}$'}, 'Interpreter','latex','FontSize',10)
grid on


%exportgraphics(gcf, 'HP_plant_full_cycle.pdf','ContentType','vector');





% Plot H at different penstock locations over time, divided by operating mode

idx = round(linspace(1,Nx,10));
labels = compose('x = %.1f m', x(idx));

idx_island = t_hist_anim < t_sync;
idx_sync   = t_hist_anim >= t_sync & t_hist_anim < t_disc;
idx_disc   = t_hist_anim >= t_disc;

% Island operation
figure
set(gcf,'Units','centimeters');
set(gcf,'Position',[0 0 12.5 8.5]);
hold on

for k = 1:length(idx)
    plot(t_hist_anim(idx_island), H_hist(idx(k),idx_island),'LineWidth',1)
end

lgd = legend(labels, 'Interpreter','latex', ...
       'Location','southeast', ...
       'FontSize',10);
lgd.ItemTokenSize = [10 6];

xlabel('t [s]','Interpreter','latex','FontSize',10)
ylabel('H [m]','Interpreter','latex','FontSize',10)
title('Hydraulic head evolution in island mode','Interpreter','latex','FontSize',10)
grid on

ax_main = gca;
ax_main.FontSize = 10;
ax_main.TickLabelInterpreter = 'latex';

axis tight
box on
ax_main.LooseInset = ax_main.TightInset;

t = t_hist_anim(idx_island);
t1 = 32;
t2 = 40;

mask = (t > t1) & (t < t2);

ax_inset = axes('Position',[0.26 0.18 0.5 0.42]);
box on
hold on

for k = 1:length(idx)
    y = H_hist(idx(k), idx_island);
    plot(t(mask), y(mask),'LineWidth',1)
end

axis tight
grid on

ax_inset.FontSize = 9;
ax_inset.TickLabelInterpreter = 'latex';

yl = ylim;
y1 = floor(yl(1));
y2 = ceil(yl(2));

ylim([y1 y2])
yticks(y1:y2)

% Export
%exportgraphics(gcf,'HP_H_island.pdf','ContentType','vector','BackgroundColor','none');


% Synchronous operation
figure
set(gcf,'Units','centimeters');
set(gcf,'Position',[0 0 12.5 8.5]);
hold on

for k = 1:length(idx)
    plot(t_hist_anim(idx_sync), H_hist(idx(k),idx_sync),'LineWidth',1)
end

lgd = legend(labels, 'Interpreter','latex', ...
       'Location','southeast', ...
       'FontSize',10);
lgd.ItemTokenSize = [10 6];

xlabel('t [s]','Interpreter','latex','FontSize',10)
ylabel('H [m]','Interpreter','latex','FontSize',10)
title('Hydraulic head evolution in synchronous mode','Interpreter','latex','FontSize',10)
grid on

ax = gca;
ax.FontSize = 10;
ax.TickLabelInterpreter = 'latex';

axis tight
box on
ax.LooseInset = ax.TightInset;

% Export
%exportgraphics(gcf,'HP_H_sync.pdf','ContentType','vector','BackgroundColor','none');


% Disconnection test
figure
set(gcf,'Units','centimeters');
set(gcf,'Position',[0 0 12.5 8.5]);
hold on

for k = 1:length(idx)
    plot(t_hist_anim(idx_disc), H_hist(idx(k),idx_disc),'LineWidth',1)
end

lgd = legend(labels, 'Interpreter','latex', ...
       'Location','northeast', ...
       'FontSize',10);
lgd.ItemTokenSize = [10 6];

xlabel('t [s]','Interpreter','latex','FontSize',10)
ylabel('H [m]','Interpreter','latex','FontSize',10)
title('Hydraulic head evolution during disconnection test','Interpreter','latex','FontSize',10)
grid on

ax_main = gca;
ax_main.FontSize = 10;
ax_main.TickLabelInterpreter = 'latex';

axis tight
box on
ax_main.LooseInset = ax_main.TightInset;

% Time window for selection
t = t_hist_anim(idx_disc);
t1 = 154;
t2 = 162;
mask = (t > t1) & (t < t2);

ax_inset = axes('Position',[0.2 0.42 0.55 0.45]);
box on
hold on

for k = 1:length(idx)
    y = H_hist(idx(k), idx_disc);
    plot(t(mask), y(mask),'LineWidth',1)
end

axis tight
grid on

ax_inset.FontSize = 9;
ax_inset.TickLabelInterpreter = 'latex';

yl = ylim;
y1 = floor(yl(1));
y2 = ceil(yl(2));

ylim([y1 y2])
yticks(y1:y2)

% Export
%exportgraphics(gcf,'HP_H_disconnection.pdf','ContentType','vector','BackgroundColor','none');



% Animation

fig_anim = figure('Name','MOC wave animation');

for k = 1:k_anim

    if ~isvalid(fig_anim)
        break
    end

    figure(fig_anim)

    subplot(2,1,1)
    plot(x, H_hist(:,k), 'LineWidth', 2)
    ylabel('H [m]')
    ylim([100 210])
    grid on

    subplot(2,1,2)
    plot(x, V_hist(:,k), 'LineWidth', 2)
    xlabel('x [m]')
    ylabel('V [m/s]')
    ylim([0 3])
    grid on

    sgtitle(sprintf('t = %.3f s', t_hist_anim(k)))

    pause(0.0001)
end