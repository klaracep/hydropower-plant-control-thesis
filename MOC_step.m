function [Hnew, Vnew, Q_T, H_b] = MOC_step(H, V, u, p)

% Penstock hydraulics using the Method of Characteristics
% Left boundary - reservoir, right boundary - turbine

Nx = p.Nx;
g  = p.g;
a  = p.a;
f  = p.f;
D  = p.D;
A  = p.A;
k  = p.k;
HL = p.HL;
Ht = p.Ht;
dt = p.dt;

Hnew = H;
Vnew = V;

% Left boundary - reservoir
Hnew(1) = HL;
Cm = V(2) - g/a*H(2) ...
     - f*dt/(2*D)*V(2)*abs(V(2));
Vnew(1) = Cm + g/a*HL;

% Right boundary - turbine
Cp = V(end-1) + g/a*H(end-1) ...
     - f*dt/(2*D)*V(end-1)*abs(V(end-1));

alpha = (g*A^2)/(a*k^2*u^2);
beta  = 1;
gamma = (g/a)*Ht - Cp;

disc = beta^2 - 4*alpha*gamma;

if disc < 0
    fprintf('disc<0 at u=%.5f, Cp=%.5f, alpha=%.3e, gamma=%.5f\n', u, Cp, alpha, gamma);
end

if disc <= 0
    Vb = 0;
    Hnew(end) = Ht;
else
    Vb = (-beta + sqrt(disc))/(2*alpha);
    Hnew(end) = Ht + (A*Vb/(k*u))^2;
end

Vnew(end) = Vb;
Hnew(end) = Ht + (A*Vb/(k*u))^2;

% Output
Q_T = A*Vb;
H_b = Hnew(end);

% Internal nodes
for i = 2:Nx-1
    Cp = V(i-1) + g/a*H(i-1) ...
         - f*dt/(2*D)*V(i-1)*abs(V(i-1));

    Cm = V(i+1) - g/a*H(i+1) ...
         - f*dt/(2*D)*V(i+1)*abs(V(i+1));

    Vnew(i) = 0.5*(Cp + Cm);
    Hnew(i) = (a/(2*g))*(Cp - Cm);
end

end