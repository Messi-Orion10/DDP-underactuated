function L = cost(x, u, i, param)
    P = param.P;
    R = param.R;
    x_star = param.x_star;

    L = 1/2*(x-x_star).'*P(:, :, i)*(x-x_star) + 1/2*u.'*R*u;
end