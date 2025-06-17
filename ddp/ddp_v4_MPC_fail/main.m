clear, clc, close all

%% parameters
N = 400;
H = 300; % prediction horizon
param.N = N;
param.H = H;
param.T = 0.01;
param.iter = 1000;
param.deltaJ_min = 1e-1;
param.J_div = 1e10;

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
param.act = 1; % pendubot = 1, acrobot = 2

% references
x0 = [0; 0; 0; 0]; % down-down
param.x_star = [pi; 0; 0; 0]; % up-up
param.lim = [-1, 1];

% weights
P0 = diag([0; 0; 10; 10]);
P1 = diag([600; 400; 10; 10]);
T1 = 100;
P2 = diag([20; 20; 10; 10]);
P = repmat(P0, 1, 1, N+H);
P(:, :, H+2-T1:H+1) = repmat(P1, 1, 1, T1);
P(:, :, H+2:end) = repmat(P2, 1, 1, N-1); % TEST
param.P = P;
param.R = 0.01;
param.mu = 1;
param.alpha = 0.1;

%% iterations
U = NaN(1, N+1);
X = NaN(4, N+1);
X(:, 1) = x0;
U_init = zeros(1, H);

wb = waitbar(0, 'Simulating...');

tic
for t = 1:N
    param.P = P(:, :, t:t+H);
    param.P(:, :, H)
    x_t = X(:, t);

    [U_t, J] = ddp_mpc(x_t, param);
    u_t = U_t(1, 1);

    F = dyn(x_t, u_t, param);
    U(:, t) = u_t;
    X(:, t+1) = F;

    waitbar(t/N, wb, sprintf('Simulating... [t = %d, cost = %.2f]', t, J));
end
toc

close(wb);

%% graphs
graphs;


% %% test
% X_TEST = NaN(4, H+1);
% X_TEST(:, 1) = x0;
% for i = 1:H
%     X_TEST(:, i+1) = dyn(X_TEST(:, i), U_t(:, i), param);
% end
% 
% 
% t = 0:0.01:H*T;
% mem = H;
% 
% X(1:2, :) = wrapToPi(X(1:2, :));
% 
% colors = {[0 0.4470 0.7410], [0.8500 0.3250 0.0980], [0.4660 0.6740 0.1880], [38, 38, 38]/255, [0.9290 0.6940 0.1250]};
% 
% fig3 = figure(3);
% set(gcf, 'Color', 'w');
% set(fig3, 'Position', [40, 100, 750, 600]);
% hold on;
% grid on;
% axis equal;
% axis([-1.2 1.2 -1.2 1.2]);
% 
% link1 = plot(NaN, NaN, 'LineWidth', 3, 'Color', colors{1});
% link2 = plot(NaN, NaN, 'LineWidth', 3, 'Color', colors{2});
% joint1 = plot(0, 0, 'o', 'MarkerSize', 10, 'MarkerFaceColor', colors{1});
% joint2 = plot(NaN, NaN, 'o', 'MarkerSize', 10, 'MarkerFaceColor', colors{2});
% ee_traj = plot(NaN, NaN, ':', 'LineWidth', 1.5, 'color', colors{4});
% 
% elapsed_time = text(-1.1, -1.1, 'Time: 0.00 s', 'FontSize', 12, 'Color', colors{4});
% 
% for i = 1:3:H+1
%     set(link1, 'XData', [0, l1*sin(X(1,i))], 'YData', [0, -l1*cos(X(1,i))]);
%     set(link2, 'XData', [l1*sin(X(1,i)), l1*sin(X(1,i)) + l2*sin(X(1, i) + X(2,i))], ...
%         'YData', [-l1*cos(X(1,i)), -l1*cos(X(1,i)) - l2*cos(X(1, i) + X(2,i))]);
%     set(joint2, 'XData', l1*sin(X(1,i)), 'YData', -l1*cos(X(1,i)));
%     set(ee_traj, 'XData', l1*sin(X(1,max([1, i-mem]):i)) + l2*sin(X(1, max([1, i-mem]):i) + X(2,max([1, i-mem]):i)), ...
%         'YData', -l1*cos(X(1,max([1, i-mem]):i)) - l2*cos(X(1, max([1, i-mem]):i) + X(2,max([1, i-mem]):i)));
%     set(elapsed_time, 'String', sprintf('Time: %.2f s', t(i)));
% 
%     drawnow;
% end