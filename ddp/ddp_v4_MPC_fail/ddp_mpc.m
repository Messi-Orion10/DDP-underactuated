function [U_opt, J_opt] = ddp_mpc(x, param)
    H = param.H;
    J_div = param.J_div;
    deltaJ_min = param.deltaJ_min;
    iter = param.iter;

    U = NaN(1, H+1, iter);
    U_init = zeros(1, H);
    U(1, 1:end-1, 1) = U_init;
    
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
    
    J = NaN(iter, 1);
    J(1) = cost_tot(X_1, U_1, param);
    %fprintf('Finding a solution...\nIteration 1 (FE) [Cost: %.2f]\n', J(1));
    
    for j = 1:iter-1
        X_j = X(:, :, j);
        U_j = U(1, :, j);
    
        [k, K] = bwd(X_j, U_j, param);
        [X_jj, U_jj] = fwd(X_j, U_j, k, K, param);
        J(j+1) = cost_tot(X_jj, U_jj, param);
    
        conv = abs(J(j+1) - J(j)) < deltaJ_min*J(j);
        div = J(j+1) > J_div;
        if(conv || div)
            J = J(1:j);
            X = X(:, :, 1:j);
            U = U(1, :, 1:j);
            %fprintf('Early stop!\n');
            break;
        end
    
        %fprintf('Iteration %d [Cost: %.2f]\n', j+1, J(j+1));
    
        X(:, :, j+1) = X_jj;
        U(1, :, j+1) = U_jj;
    end
    
    J_opt = J(end);
    U_opt = U(1, :, end);
end