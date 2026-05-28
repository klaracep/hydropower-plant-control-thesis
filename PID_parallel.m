function [u_sat, X] = PID_parallel(w, y, dt, X, P)
% PID controller function
% Tustin method

r0 = P.r0;
ri = P.ri;
rd = P.rd;

% Weights of the derivative term
b  = P.b;
c  = P.c;

N  = P.N; % Filter
K  = P.K; % Anti-windup coefficient

u_max = P.u_max;
u_min = P.u_min;
dut_max = P.dut_max;

% States
I = X(1);
x = X(2);
u_prev = X(3);
fI_prev = X(4);

% Derivative filter (Tustin)
a = dt*N/2;
x = ((1-a)*x + dt*N*(c*w - y))/(1+a);

D = N*(-x + c*w - y);

% Proportional term
Pterm = r0*(b*w - y);

% Controller output
u = Pterm + ri*I + rd*D;

% Rate limiter
du = u - u_prev;
du = max(min(du, dut_max*dt), -dut_max*dt);
u_rl = u_prev + du;

% Saturation
u_sat = min(max(u_rl, u_min), u_max);

% Integrator
fI = (w - y) + K*(u_sat - u);
I = I + dt/2*(fI + fI_prev);


X = [I x u_sat fI];

end