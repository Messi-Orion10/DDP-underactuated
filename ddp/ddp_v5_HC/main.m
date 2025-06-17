clear, clc, close all

%% parameters
N = 300;
param.N = N;
param.T = 0.01;
iter = 1000;
deltaJ_min = 1e-13;
J_div = 1e12;

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
param.constr = 3; % 0 = null, saturation = 1, squashing = 2, QP = 3

% weights
P0 = diag([0; 0; 10; 10]);
Pf = diag([600; 400; 20; 20]);
Tf = 100;
P = NaN(4, 4, N+1);
P(:, :, 1:N+1-Tf) = repmat(P0, 1, 1, N+1-Tf);
P(:, :, N+1-Tf+1:end) = repmat(Pf, 1, 1, Tf);
param.P = P;
param.R = 1;
param.mu = 0;
param.alpha0 = 1;

%% iterations
x0 = param.x0;
alpha0 = param.alpha0;
beta = param.beta;

U = NaN(1, N+1, iter);
U_init = zeros(1, N);
U(1, 1:end-1, 1) = U_init;

X = NaN(4, N+1, iter);
X(:, 1, 1) = x0;

for i = 1:N
    x_i = X(:, i, 1);
    u_i = U(1, i, 1);

    F = dyn(x_i, u_i, param);

    X(:, i+1, 1) = F;
end

X_1 = X(:, :, 1);
U_1 = U(1, :, 1);

tic
J = NaN(iter, 1);
J(1) = cost_tot(X_1, U_1, param);
fprintf('Finding a solution...\nIteration 1 (FE) [Cost: %.2f]\n', J(1));

for j = 1:iter-1
    X_j = X(:, :, j);
    U_j = U(1, :, j);
    alpha = alpha0/(1+beta*(j-1));

    [k, K] = bwd(X_j, U_j, param);
    [X_jj, U_jj] = fwd(X_j, U_j, k, K, alpha, param);
    J(j+1) = cost_tot(X_jj, U_jj, param);

    conv = abs(J(j+1) - J(j)) < deltaJ_min*J(j);
    div = J(j+1) > J_div;
    if(conv || div)
        J = J(1:j);
        X = X(:, :, 1:j);
        U = U(1, :, 1:j);
        fprintf('Early stop!\n');
        break;
    end

    fprintf('Iteration %d [Cost: %.2f]\n', j+1, J(j+1));

    X(:, :, j+1) = X_jj;
    U(1, :, j+1) = U_jj;
end
toc

%% graphs
graphs;