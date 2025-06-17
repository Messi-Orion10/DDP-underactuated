function [Lx, Lu, Lxx, Luu] = der_L(x, u, P, param)
    R = param.R;
    x_star = param.x_star;

    Lx = (x - x_star).'*P;
    Lu = u*R;
    Lxx = P;
    Luu = R;
end