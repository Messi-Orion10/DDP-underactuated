function L = L(x, u, param)
    P = param.P;
    R = param.R;
    x_star = param.x_star;

    L = 1/2*(x-x_star).'*P*(x-x_star) + 1/2*u.'*R*u;
end