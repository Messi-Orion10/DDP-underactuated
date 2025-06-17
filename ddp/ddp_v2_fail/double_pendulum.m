clear, clc, close all
%% Parameters
T = 0.01; 

I1 = 1;
I2 = I1;
m1 = 0.5;
m2 = 0.75;
l1 = 0.4;
d1 = l1/2;
l2 = 0.6;
d2 = l2/2;
g = 9.81;
b1 = 0.2;
b2 = b1;
N = 1000;
P = diag([300;4;1;1]);
Pf = diag([300;10;1;1]);
R = diag([0.1; 0.2]);
J = 2;

a1 = I1 + m1*d1^2 + I2 + m2*(l1^2 + d2^2);
a2 = m2*l1*d2;
a3 = I2 + m2*d2^2;
a4 = g*(m1*d1 + m2*l1);
a5 = g*m2*d2;

x_star = [pi; 0; 0; 0]; % up-up (swing-up problem)

mu = 0;

param.T = T;
param.a1 = a1;
param.a2 = a2;
param.a3 = a3;
param.a4 = a4;
param.a5 = a5;
param.b1 = b1;
param.b2 = b2;
param.P = P;
param.Pf = Pf;
param.R = R;
param.N = N;
param.x_star = x_star;
param.mu = mu;

t = 0:0.01:N*T;
mem = N;
%% Dynamics 
x0 = [-pi/2; 0; 0; 0]; %down-down (initial condition)
U_init = zeros(2,N);
for i = 1:N
    U_init(1,i) = 0;
    U_init(2,i) = 0;
end

U = NaN(2,N+1,J);
U(:,1:end-1, 1) = U_init;
X = NaN(4,N+1,J);
X(:, 1, 1) = x0;

for i = 1:N
    X(:,i+1, 1) = F(X(:,i, 1), U(:,i, 1), param);
end

for j = 1:J-1
    [k, K] = bwd(X(:, :, j), U(:, :, j), param);
    [X(:, :, j+1), U(:, :, j+1)] = fwd(X(:, :, j), U(:, :, j), k, K, param);
    fprintf('|');
end

%% Plots
X_plot = X(:, :, J);
X_plot(1:2, :) = wrapToPi(X_plot(1:2, :));

% State trajectories
colors = {[0 0.4470 0.7410], [0.8500 0.3250 0.0980], [38, 38, 38]/255};
%colors = {[0, 112, 192]/255, [38, 38, 38]/255, [38, 38, 38]/255, [38, 38, 38]/255};

fig1 = figure(1);
set(gcf, 'Color', 'w');
set(fig1, 'Position', [750, 100, 750, 600]);
hold on;
grid on;
xlabel('Time [s]');
ylabel('Trajectories');
title('State trajectories');

state = gobjects(N, 1);
for j=1:2
    state(j) = plot(t, X_plot(j,:), 'LineWidth', 2, 'Color', colors{j});
end

for j=1:2
    state(j+2) = plot(t, X_plot(j+2,:), ':', 'LineWidth', 2, 'Color', colors{j});
end

%% Animation
fig2 = figure(2);
set(gcf, 'Color', 'w');
set(fig2, 'Position', [100, 100, 750, 600]);
hold on;
grid on;
axis equal;
axis([-1.2 1.2 -1.2 1.2]);

% Plot links and joints
link1 = plot(NaN, NaN, 'LineWidth', 2, 'Color', colors{1});
link2 = plot(NaN, NaN, 'LineWidth', 2, 'Color', colors{2});
joint1 = plot(NaN, NaN, 'o', 'MarkerSize', 10, 'MarkerFaceColor', colors{1});
joint2 = plot(NaN, NaN, 'o', 'MarkerSize', 10, 'MarkerFaceColor', colors{2});
end_effector = plot(NaN, NaN, ':', 'LineWidth', 2, 'color', colors{3});

elapsed_time = text(-1.1, -1.1, 'Time: 0.00 s', 'FontSize', 12, 'Color', 'k');

% Animation loop
for i = 1:3:N+1
    set(link1, 'XData', [0, l1*sin(X_plot(1,i))], 'YData', [0, -l1*cos(X_plot(1,i))]);
    set(link2, 'XData', [l1*sin(X_plot(1,i)), l1*sin(X_plot(1,i)) + l2*sin(X_plot(1, i) + X_plot(2,i))], ...
        'YData', [-l1*cos(X_plot(1,i)), -l1*cos(X_plot(1,i)) - l2*cos(X_plot(1, i) + X_plot(2,i))]);
    set(joint2, 'XData', l1*sin(X_plot(1,i)), 'YData', -l1*cos(X_plot(1,i)));
    set(end_effector, 'XData', l1*sin(X_plot(1,max([1, i-mem]):i)) + l2*sin(X_plot(1, max([1, i-mem]):i) + X_plot(2,max([1, i-mem]):i)), ...
        'YData', -l1*cos(X_plot(1,max([1, i-mem]):i)) - l2*cos(X_plot(1, max([1, i-mem]):i) + X_plot(2,max([1, i-mem]):i)));
    set(elapsed_time, 'String', sprintf('Time: %.2f s', t(i)));
    
    drawnow;
end
