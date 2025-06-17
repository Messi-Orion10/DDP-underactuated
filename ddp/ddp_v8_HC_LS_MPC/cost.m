function J = cost(X, U, P, param)
    H = param.H;
    R = param.R;
    x_star = param.x_star;

    U(isnan(U)) = 0;
    J = 0;
    for i = 1:(H+1)
        x_i = X(:, i);
        u_i = U(1, i);
        P_i = P(:, :, i);

        J = J + 1/2*(x_i-x_star).'*P_i*(x_i-x_star) + 1/2*u_i.'*R*u_i;
    end
end