function [U_opt, J_opt] = ddp_mpc(x, U_init, P, iter, param)
    H = param.H;
    deltaJ_min = param.deltaJ_min;
    alpha0 = param.alpha0;
    beta = param.beta;
    alpha_min = param.alpha_min;
    verbose = param.verbose;

    U(1, :, 1) = U_init;
    
    X = NaN(4, H+1, iter);
    X(:, 1, 1) = x;
    
    for i = 1:H
        x_i = X(:, i, 1);
        u_i = U(1, i, 1);
    
        F = dyn(x_i, u_i, param);
    
        X(:, i+1, 1) = F;
    end
    
    X_1 = X(:, :, 1);
    U_1 = U(1, :, 1);
    % disp(X_1(:, 1)) %
    % disp(U_1(1, 1)) %
    
    J = NaN(iter, 1);
    J(1) = cost(X_1, U_1, P, param);
    if verbose
        fprintf('Finding a solution...\nIteration 1 [Cost: %.2f]\n', J(1));
    end
    
    conv = 0;
    for j = 1:iter-1
        X_j = X(:, :, j);
        U_j = U(1, :, j);
        
        [k, K] = bwd(X_j, U_j, P, param);

        alpha = alpha0;
        [X_jj, U_jj] = fwd(X_j, U_j, k, K, alpha, param);
        J(j+1) = cost(X_jj, U_jj, P, param);
        while J(j+1) > J(j) - deltaJ_min || any(any(isnan(X_jj)))
            alpha = alpha*beta;
            [X_jj, U_jj] = fwd(X_j, U_j, k, K, alpha, param);
            J(j+1) = cost(X_jj, U_jj, P, param);
            if(alpha < alpha_min)
                conv = 1;
                break;
            end
        end

        if(conv)
            J = J(1:j);
            X = X(:, :, 1:j);
            U = U(1, :, 1:j);
            % disp(X(:, 2, end)); %
            % disp(U(1, 2, end)); %
            break;
        end
    
        if verbose
            fprintf('Iteration %d [Cost: %.2f]\n', j+1, J(j+1));
        end
    
        X(:, :, j+1) = X_jj;
        U(1, :, j+1) = U_jj;
    end
    
    J_opt = J(end);
    U_opt = U(1, :, end);
end

