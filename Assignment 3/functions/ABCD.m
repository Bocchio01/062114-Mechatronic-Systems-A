function [A, B, C, D] = ABCD(x, u)

% Load system parameters
load("parameters.mat"); %#ok<LOAD>

L = @(T) L0 + beta1 * T + beta2 * T^2;
dLdT = @(T) beta1 + 2*beta2 * T;
ddLddT = @(T) 2*beta2; %#ok<NASGU>


% A coefficients
a11 = 0;
a12 = 1;
a13 = 0;
a14 = 0;

a21 = 1 / m * (-k1 -3*k3*x(1)^2);
a22 = 1 / m * (-c);
a23 = 1 / m * (alpha0);
a24 = 0;

a31 = 0;
a32 = 0;
a33 = -R / L(x(4));
a34 = (-R*x(3) + u) * ( -1 / L(x(4))^2 ) * dLdT(x(4));

a41 = 0;
a42 = 0;
a43 = 1 / C_T * (2*R*x(3));
a44 = 1 / C_T * (-h);

% B coefficients
b1 = 0;
b2 = 0;
b3 = 1 / L(x(4));
b4 = 0;

% C coefficients
c1 = 1;
c2 = 0;
c3 = 0;
c4 = 0;

% D coefficients
d1 = 0;


% State matrices
A = [
    a11 a12 a13 a14;
    a21 a22 a23 a24;
    a31 a32 a33 a34;
    a41 a42 a43 a44
    ];

B = [
    b1;
    b2;
    b3;
    b4
    ];

C = [c1 c2 c3 c4];

D = d1;

end