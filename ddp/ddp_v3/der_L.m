function [Lx, Lu, Lxx, Luu] = der_L(x, u, i, param)
    P = param.P;
    R = param.R;
    x_star = param.x_star;

    P_i = P(:, :, i);

    Lx = (x - x_star).'*P_i;
    Lu = u*R;
    Lxx = P(:, :, i);
    Luu = R;
end