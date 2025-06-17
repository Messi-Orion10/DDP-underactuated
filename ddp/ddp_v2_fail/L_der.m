function [Lx, Lu, Lxx, Luu] = L_der(x, u, param)
    x_star = param.x_star;
    P = param.P;
    R = param.R;
    x1 = x(1);
    x2 = x(2);
    x3 = x(3);
    x4 = x(4);
    x1_star = x_star(1);
    x2_star = x_star(2);
    x3_star = x_star(3);
    x4_star = x_star(4);
    P1 = P(1, 1);
    P2 = P(2, 2);
    P3 = P(3, 3);
    P4 = P(4, 4);

    Lx = [P1*sin(x1-x1_star), P2*sin(x2-x2_star), P3*(x3-x3_star), P4*(x4-x4_star)];
    Lxx = diag([P1*cos(x1-x1_star), P2*cos(x2-x2_star), P3, P4]);
    Lu = u.'*R;
    Luu = R;
end