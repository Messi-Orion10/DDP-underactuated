function X = test(x0, U, param)
    H = param.H;
    l1 = 0.5;
    l2 = 0.5;
    T = param.T;
    N = param.N;

    %% rollout
    x_i = x0;
    X(:, 1) = x_i;
    for i = 1:H
        u_i = U(1, i);
        x_ii = dyn(x_i, u_i, param);
    
        X(:, i+1) = x_ii;
        x_i = x_ii;
    end
    
    %%
    t = 0:0.01:N*T;
    mem = N;
    
    X(1:2, :) = wrapToPi(X(1:2, :));
    
    colors = {[0 0.4470 0.7410], [0.8500 0.3250 0.0980], [0.4660 0.6740 0.1880], [38, 38, 38]/255, [0.9290 0.6940 0.1250]};
    
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
    
    for i = 1:3:H+1
        set(link1, 'XData', [0, l1*sin(X(1,i))], 'YData', [0, -l1*cos(X(1,i))]);
        set(link2, 'XData', [l1*sin(X(1,i)), l1*sin(X(1,i)) + l2*sin(X(1, i) + X(2,i))], ...
            'YData', [-l1*cos(X(1,i)), -l1*cos(X(1,i)) - l2*cos(X(1, i) + X(2,i))]);
        set(joint2, 'XData', l1*sin(X(1,i)), 'YData', -l1*cos(X(1,i)));
        set(ee_traj, 'XData', l1*sin(X(1,max([1, i-mem]):i)) + l2*sin(X(1, max([1, i-mem]):i) + X(2,max([1, i-mem]):i)), ...
            'YData', -l1*cos(X(1,max([1, i-mem]):i)) - l2*cos(X(1, max([1, i-mem]):i) + X(2,max([1, i-mem]):i)));
        set(elapsed_time, 'String', sprintf('Time: %.2f s', t(i)));
        
        drawnow;
    end

    close(gcf);
end