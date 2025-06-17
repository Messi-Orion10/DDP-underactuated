clear, clc, close all

%% parameters
N = 400;
H = 300;
iter_init = 1000;
iter = 10;
param.T = 0.01;
param.N = N;
param.H = H;
param.deltaJ_min = 1e-13;
param.J_div = 1e12;

% robot
g = 9.81;
l1 = 0.5;
l2 = 0.5;
m1 = 2;
m2 = 2;
I1 = 1/12*m1*l1^2;
I2 = 1/12*m2*l2^2;
d1 = l1/2;
d2 = l2/2;
param.a1 = I1 + m1*d1^2 + I2 + m2*(l1^2 + d2^2);
param.a2 = m2*l1*d2;
param.a3 = I2 + m2*d2^2;
param.a4 = g*(m1*d1 + m2*l1);
param.a5 = g*m2*d2;
param.b1 = 0.1;
param.b2 = 0.1;
param.act = 2; % pendubot = 1, acrobot = 2

% references
param.x0 = [0; 0; 0; 0]; % down-down
param.x_star = [pi; 0; 0; 0]; % up-up
param.lim = [-7, 7];
param.constr = 0; % 0 = null, saturation = 1, squashing = 2, QP = 3

% weights
P0 = diag([0; 0; 10; 10]);
P1 = diag([600; 400; 20; 20]);
P2 = diag([300; 200; 10; 10]);
P = NaN(4, 4, N+1);
P(:, :, 1:201) = repmat(P0, 1, 1, 201);
P(:, :, 202:301) = repmat(P1, 1, 1, 100);
P(:, :, 302:701) = repmat(P2, 1, 1, 400);
param.P = P;
param.R = 1;
param.mu = 0;
param.alpha0 = 0.01;
param.beta = 0;

%% initial control sequence
x0 = param.x0;

U_null = zeros(1, H);
U_null = [U_null, NaN];

tic
param.P = P(:, :, 1:H+1);
[U_init, J_init] = ddp_mpc(x0, U_null, iter_init, 1, param);
toc

%% mpc iterations
U = NaN(1, N+1);
X = NaN(4, N+1);
X(:, 1) = x0;

wb = waitbar(0, 'Simulating...');
tic
for t = 1:N
    param.P = P(:, :, t:t+H);
    x_t = X(:, t);

    [U_t, J] = ddp_mpc(x_t, U_init, iter, 0, param);
    u_t = U_t(1, 1);

    F = dyn(x_t, u_t, param);
    U(:, t) = u_t;
    X(:, t+1) = F;

    U_init = [U_init(1, 2:end-1), 0, NaN];

    waitbar(t/N, wb, sprintf('Simulating... [t = %d, cost = %.2f]', t, J));
end
toc
close(wb);

%% graphs
graphs;