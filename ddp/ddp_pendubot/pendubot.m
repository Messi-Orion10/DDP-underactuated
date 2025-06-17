clear, clc, close all
%% Parameters
T = 0.01; 

I1 = 0.9;
I2 = I1;
m1 = 1;
m2 = 1;
l1 = 0.5;
d1 = 0.5;
l2 = 0.5;
d2 = 0.5;
g = 9.81;
b1 = 0.1;
b2 = b1;
N = 1000;
Tf = 200; % <= N+1, >= 1
P = diag([0;0;10;10]);
Pf = diag([500;500;120;120]);
R = 10;
Jter = 20; % n iter = Jter-1

a1 = I1 + m1*d1^2 + I2 + m2*(l1^2 + d2^2);
a2 = m2*l1*d2;
a3 = I2 + m2*d2^2;
a4 = g*(m1*d1 + m2*l1);
a5 = g*m2*d2;

x0 = [0; 0; 0; 0]; %down-down (initial condition)
x_star = [pi; 0; 0; 0]; % up-up (swing-up problem)

mu = 0.000;

param.T = T;
param.a1 = a1;
param.a2 = a2;
param.a3 = a3;
param.a4 = a4;
param.a5 = a5;
param.b1 = b1;
param.b2 = b2;
param.Tf = Tf;
param.P = P;
param.Pf = Pf;
param.R = R;
param.N = N;
param.x_star = x_star;
param.mu = mu;

t = 0:0.01:N*T;
mem = N;
%% Dynamics 
U_init = zeros(1,N);
for i = 1:N
    U_init(i) = 0;
end

U = NaN(1,N+1,Jter);
U(1,1:end-1, 1) = U_init;
X = NaN(4,N+1,Jter);
X(:, 1, 1) = x0;

for i = 1:N
    X(:, i+1, 1) = F(X(:, i, 1), U(1, i, 1), param);
end

for j = 1:Jter-1
    [k, K] = bwd(X(:, :, j), U(1, :, j), param);
    [X(:, :, j+1), U(1, :, j+1)] = fwd(X(:, :, j), U(1, :, j), k, K, param);
    cost = J(X(:, :, j+1), U(1, :, j+1), param);
    fprintf('Iteration %d (Cost: %.2f)\n', j, cost);
end

%% Plots
X_plot = X(:, :, Jter);
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
