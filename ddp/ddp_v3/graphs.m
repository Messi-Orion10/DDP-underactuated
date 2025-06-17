T = param.T;
N = param.N;
lim = param.lim;

t = 0:0.01:N*T;
mem = N;

X_plot = X(:, :, end);
X_plot(1:2, :) = wrapToPi(X_plot(1:2, :));

U_plot = U(:, :, end);

colors = {[0 0.4470 0.7410], [0.8500 0.3250 0.0980], [0.4660 0.6740 0.1880], [38, 38, 38]/255, [0.9290 0.6940 0.1250]};

%% cost per iteration
fig1 = figure(1);
set(gcf, 'Color', 'w');
set(fig1, 'Position', [750, 100, 750, 600]);
hold on;
grid on;
xlabel('Iteration', 'Interpreter', 'latex');
ylabel('Cost', 'Interpreter', 'latex');
title('Cost per iteration', 'Interpreter', 'latex');
cost_iter = plot(1:size(J, 1), J, 'LineWidth', 2, 'Color', colors{5});

%% trajectories
fig2 = figure(2);
set(gcf, 'Color', 'w');
set(fig2, 'Position', [750, 100, 750, 600]);

subplot(2, 1, 1);
hold on;
grid on;
xlabel('Time [s]', 'Interpreter', 'latex');
ylabel('Angle [rad] / Angular velocity [rad/s]', 'Interpreter', 'latex');
title('State trajectories', 'Interpreter', 'latex');
legend('show', 'Location', 'best', 'Interpreter', 'latex');

states = gobjects(N, 1);
for j = 1:2
    states(j) = plot(t, X_plot(j,:), 'LineWidth', 2, 'Color', colors{j}, 'DisplayName', sprintf('$x_%d$', j));
end
for j = 1:2
    states(j + 2) = plot(t, X_plot(j + 2,:), ':', 'LineWidth', 2, 'Color', colors{j}, 'DisplayName', sprintf('$x_%d$', j + 2));
end

subplot(2, 1, 2);
hold on;
grid on;
xlabel('Time [s]', 'Interpreter', 'latex');
ylabel('Torque [Nm]', 'Interpreter', 'latex');
title('Control sequence', 'Interpreter', 'latex');
legend('show', 'Location', 'best', 'Interpreter', 'latex');

control = plot(t, U_plot(1,:), 'LineWidth', 2, 'Color', colors{3}, 'DisplayName', '$u$');
lim_lo = yline(lim(1), ':', 'LineWidth', 2, 'DisplayName', '$\underline{u}$', 'Color', colors{1});
lim_up = yline(lim(2), ':', 'LineWidth', 2, 'DisplayName', '$\overline{u}$', 'Color', colors{2});

%% animation
fig3 = figure(3);
set(gcf, 'Color', 'w');
set(fig3, 'Position', [40, 100, 750, 600]);
hold on;
grid on;
axis equal;
axis([-1.2 1.2 -1.2 1.2]);

link1 = plot(NaN, NaN, 'LineWidth', 3, 'Color', colors{1});
link2 = plot(NaN, NaN, 'LineWidth', 3, 'Color', colors{2});
joint1 = plot(0, 0, 'o', 'MarkerSize', 10, 'MarkerFaceColor', colors{1});
joint2 = plot(NaN, NaN, 'o', 'MarkerSize', 10, 'MarkerFaceColor', colors{2});
ee_traj = plot(NaN, NaN, ':', 'LineWidth', 1.5, 'color', colors{4});

elapsed_time = text(-1.1, -1.1, 'Time: 0.00 s', 'FontSize', 12, 'Color', colors{4});

for i = 1:3:N+1
    set(link1, 'XData', [0, l1*sin(X_plot(1,i))], 'YData', [0, -l1*cos(X_plot(1,i))]);
    set(link2, 'XData', [l1*sin(X_plot(1,i)), l1*sin(X_plot(1,i)) + l2*sin(X_plot(1, i) + X_plot(2,i))], ...
        'YData', [-l1*cos(X_plot(1,i)), -l1*cos(X_plot(1,i)) - l2*cos(X_plot(1, i) + X_plot(2,i))]);
    set(joint2, 'XData', l1*sin(X_plot(1,i)), 'YData', -l1*cos(X_plot(1,i)));
    set(ee_traj, 'XData', l1*sin(X_plot(1,max([1, i-mem]):i)) + l2*sin(X_plot(1, max([1, i-mem]):i) + X_plot(2,max([1, i-mem]):i)), ...
        'YData', -l1*cos(X_plot(1,max([1, i-mem]):i)) - l2*cos(X_plot(1, max([1, i-mem]):i) + X_plot(2,max([1, i-mem]):i)));
    set(elapsed_time, 'String', sprintf('Time: %.2f s', t(i)));
    
    drawnow;
end
