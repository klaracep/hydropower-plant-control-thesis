function [Pi, Ps] = params_pid()

% PID island
Pi.r0 = 3.43;
Pi.ri = 1.43;
Pi.rd = 0.391;
Pi.b = 0.5;
Pi.c = 0;
Pi.N = 1/0.04;
Pi.u_max = 1;
Pi.u_min = 0.002; % Small positive opening for numerical stability
Pi.dut_max = 0.2;
Pi.K = 10;  % anti wind-up


% PID sync
Ps.r0 =  0.849;
Ps.ri =  0.98;
Ps.rd =  0.042;
Ps.b = 0.5;
Ps.c = 0;
Ps.N = 1/0.04;
Ps.u_max = 1;
Ps.u_min = 0.002;
Ps.dut_max = 0.2;
Ps.K = 10;

end
