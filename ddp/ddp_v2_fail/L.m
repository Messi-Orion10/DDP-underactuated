function L = L(x, u, param)
    P = param.P;
    R = param.R;
    P1 = P(1, 1);
    P2 = P(2, 2);
    P34 = P(3:4, 3:4);
    x_star = param.x_star;
    x1 = x(1);
    x2 = x(2);
    x1_star = x_star(1);
    x2_star = x_star(2);
    x34 = [x(3); x(4)];
    x34_star = [x_star(3); x_star(4)];

    L = 1/2*(x34-x34_star).'*P34*(x34-x34_star) + 1/2*u.'*R*u + P1*(1-cos(x1-x1_star)) + P2*(1-cos(x2-x2_star));
end