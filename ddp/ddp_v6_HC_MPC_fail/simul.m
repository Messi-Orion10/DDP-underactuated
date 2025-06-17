function X = simul(x0, U, param)
    N = param.N;

    X = NaN(4, N+1);
    X(:, 1) = x0;

    wb = waitbar(0, 'Simulating...');

    for t = 1:N
        x_t = X(:, t);
    
        u_t = U(1, t);
        F = dyn(x_t, u_t, param);
        X(:, t+1) = F;
    
        waitbar(t/N, wb, sprintf('Simulating... [t = %d]', t));
    end

    close(wb);
end

