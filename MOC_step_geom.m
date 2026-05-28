function [Hnew, Vnew, Q_T, H_b] = MOC_step_geom(H, V, u, p)
% MOC modified for a four-segment penstock
% S1 - vertical section (D = constant)
% S2 - elbow
% S3 - conical section
% S4 - horizontal section (D = constant)
% S2 is represented as an additional loss between S1 and S3

Nx = p.Nx;
g  = p.g;
a  = p.a;
f  = p.f;
dt = p.dt;
dx = p.dx;

HL = p.HL;
Ht = p.Ht;
k  = p.k;

D_start = p.D;
D_end = p.D_end;

eta = p.eta;

% Segment boundaries

i_vert_end = round(p.L_vert/dx) + 1;
i_cone_end = i_vert_end + round(p.L_cone/dx);

% Conical section geometry (S3)

Nc = i_cone_end - i_vert_end + 1;

D3 = linspace(D_start,D_end,Nc);

b = (D_end-D_start)/((Nc-1)*dx);   % dD/dx

% Initialization

Hnew = H;
Vnew = V;

% SEGMENT 1 (vertical section)

for i = 2:i_vert_end-1

    Cp = V(i-1) + g/a*H(i-1) ...
        - f*dt/(2*D_start)*V(i-1)*abs(V(i-1));

    Cm = V(i+1) - g/a*H(i+1) ...
        - f*dt/(2*D_start)*V(i+1)*abs(V(i+1));

    Vnew(i) = 0.5*(Cp + Cm);
    Hnew(i) = (a/(2*g))*(Cp - Cm);

end

% INTERFACE 1 (vertical section -> elbow -> conical section)

D = D3(1);

C1p = V(i_vert_end-1) + g/a*H(i_vert_end-1) ...
      - f*dt/(2*D_start)*V(i_vert_end-1)*abs(V(i_vert_end-1));

C2m = V(i_vert_end+1) - g/a*H(i_vert_end+1) ...
      - f*dt/(2*D)*V(i_vert_end+1)*abs(V(i_vert_end+1)) ...
      + 2*a*V(i_vert_end+1)/D*b*dt;

Vnew(i_vert_end) = 0.5*(C2m + C1p);
Hnew(i_vert_end) = (a/(2*g))*(C1p - C2m);

% S2 - elbow -> additional loss

Hnew(i_vert_end) = Hnew(i_vert_end) ...
    - eta * Vnew(i_vert_end)^2/(2*g);

% SEGMENT 3 (conical section)

for j = 2:Nc-1

    i = i_vert_end + j - 1;

    Dm = D3(j-1);
    Dp = D3(j+1);

    Cp = V(i-1) + g/a*H(i-1) ...
        - f*dt/(2*Dm)*V(i-1)*abs(V(i-1)) ...
        - 2*a*V(i-1)/Dm*b*dt;

    Cm = V(i+1) - g/a*H(i+1) ...
        - f*dt/(2*Dp)*V(i+1)*abs(V(i+1)) ...
        + 2*a*V(i+1)/Dp*b*dt;

    Vnew(i) = 0.5*(Cp + Cm);
    Hnew(i) = (a/(2*g))*(Cp - Cm);

end

% INTERFACE 3 (conical section -> horizontal section)

D = D3(end);

C1p = V(i_cone_end-1) + g/a*H(i_cone_end-1) ...
      - f*dt/(2*D)*V(i_cone_end-1)*abs(V(i_cone_end-1)) ...
      - 2*a*V(i_cone_end-1)/D*b*dt;

C2m = V(i_cone_end+1) - g/a*H(i_cone_end+1) ...
      - f*dt/(2*D_end)*V(i_cone_end+1)*abs(V(i_cone_end+1));

Vnew(i_cone_end) = 0.5*(C2m + C1p);
Hnew(i_cone_end) = (a/(2*g))*(C1p - C2m);

% SEGMENT 4 (horizontal section)

for i = i_cone_end+1:Nx-1

    Cp = V(i-1) + g/a*H(i-1) ...
        - f*dt/(2*D_end)*V(i-1)*abs(V(i-1));

    Cm = V(i+1) - g/a*H(i+1) ...
        - f*dt/(2*D_end)*V(i+1)*abs(V(i+1));

    Vnew(i) = 0.5*(Cp + Cm);
    Hnew(i) = (a/(2*g))*(Cp - Cm);

end

% Left boundary - reservoir

Hnew(1) = HL;

Cm = V(2) - g/a*H(2) ...
    - f*dt/(2*D_start)*V(2)*abs(V(2));

Vnew(1) = Cm + g/a*HL;

% Right boundary - turbine

Cp = V(end-1) + g/a*H(end-1) ...
    - f*dt/(2*D_end)*V(end-1)*abs(V(end-1));

A = pi*D_end^2/4;

alpha = (g*A^2)/(a*k^2*u^2);
beta  = 1;
gamma = (g/a)*Ht - Cp;

disc = beta^2 - 4*alpha*gamma;

if disc <= 0

    Vb = 0;
    Hnew(end) = Ht;

else

    Vb = (-beta + sqrt(disc))/(2*alpha);
    Hnew(end) = Ht + (A*Vb/(k*u))^2;

end

Vnew(end) = Vb;

Q_T = A*Vb;
H_b = Hnew(end);

end