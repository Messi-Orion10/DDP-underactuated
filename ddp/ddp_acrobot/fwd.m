function [X_hat, U_hat] = fwd(X, U, k, K, param)
    N = param.N;

    X_hat = NaN(4,N+1);
    X_hat(:, 1) = X(:,1); %x0

    U_hat = NaN(1,N+1);
    
    for i = 1:N
        U_hat(1,i) = U(1,i) + k(1,i) + K(1,:,i)*(X_hat(:,i)-X(:,i));
        X_hat(:,i+1) = F(X_hat(:,i), U_hat(1,i), param);
    end
end