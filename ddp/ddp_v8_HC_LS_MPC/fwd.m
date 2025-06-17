function [X_hat, U_hat] = fwd(X, U, k, K, alpha, param)
    H = param.H;
    constr = param.constr;
    div = param.div;

    X_hat = NaN(4, H+1);
    X_hat(:, 1) = X(:, 1);

    U_hat = NaN(1, H+1);
    
    for i = 1:H
        u_i = U(1, i);
        x_i = X(:, i);
        x_hat_i = X_hat(:, i);
        K_i = K(1, :, i);
        k_i = k(1, i);

        if(abs(x_hat_i) > div)
            break;
        end

        u_hat_i = u_i + alpha*k_i + K_i*(x_hat_i - x_i);
        if constr ~= 0
            u_hat_i = bound(u_hat_i, param);
        end
        x_hat_ii = dyn(x_hat_i, u_hat_i, param);

        X_hat(:, i+1) = x_hat_ii;
        U_hat(:, i) = u_hat_i;
    end
end