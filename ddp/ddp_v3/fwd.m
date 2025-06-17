function [X_hat, U_hat] = fwd(X, U, k, K, alpha, param)
    N = param.N;

    X_hat = NaN(4, N+1);
    X_hat(:, 1) = X(:, 1);

    U_hat = NaN(1, N+1);
    
    for i = 1:N
        u_i = U(1, i);
        x_i = X(:, i);
        x_hat_i = X_hat(:, i);
        K_i = K(1, :, i);
        k_i = k(1, i);

        u_hat_i = u_i + alpha*k_i + K_i*(x_hat_i - x_i);
        x_hat_ii = dyn(x_hat_i, u_hat_i, param);

        X_hat(:, i+1) = x_hat_ii;
        U_hat(:, i) = u_hat_i;
    end
end