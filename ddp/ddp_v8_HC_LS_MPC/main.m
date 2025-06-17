clear, clc, close all

%% parameters
N = 500;
H = 300;
iter_init = 1000;
iter = 10;
T = 0.01;
param.T = T;
param.H = H;
param.deltaJ_min = 10;
param.alpha_min = 1e-3;
param.div = 1e12;
param.verbose = 0;

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
param.lim = [-3.5, 3.5];
param.constr = 3; % 0 = null, saturation = 1, squashing = 2, QP = 3

% weights
P_light = diag([10; 10; 10; 10]);
% P_heavy = diag([100000; 100000; 1000; 1000]);
P_heavy = diag([600; 400; 30; 30]);
P = NaN(4, 4, N+H+1);
P(:, :, 1:H) = repmat(P_light, 1, 1, H);
P(:, :, H+1:end) = repmat(P_heavy, 1, 1, N+1);
P(:, :, H+2:end) = repmat(P_heavy, 1, 1, N);
param.R = 10;
param.mu = 0;
param.alpha0 = 1;
param.beta = 0.5;

%% initial control sequence
tic
U = NaN(1, N+1);
X = NaN(4, N+1);

x_1 = param.x0;
X(:, 1) = x_1;

U_null = zeros(1, H);
U_null = [U_null, NaN];

P_1 = P(:, :, 1:H+1);
[U_init, J_init] = ddp_mpc(x_1, U_null, P_1, iter_init, param);

u_1 = U_init(1, 1);
U(:, 1) = u_1;

%% mpc iterations
X(:, 1) = x_1;

wb = waitbar(0, 'Simulating...');
for t = 1:N
    P_t = P(:, :, t:t+H);
    x_t = X(:, t);

    [U_t, J] = ddp_mpc(x_t, U_init, P_t, iter, param);
    u_t = U_t(1, 1);

    d_t = 0;
    F = dyn_dist(x_t, u_t, d_t, param);

    U(:, t) = u_t;
    X(:, t+1) = F;

    U_init = [U_t(1, 2:end-1), U_t(1, end-1), NaN];

    waitbar(t/N, wb, sprintf('Simulating... [t = %d, cost-to-go = %.2f]', t, J));
end
toc
close(wb);

%% graphs
graphs_rec;