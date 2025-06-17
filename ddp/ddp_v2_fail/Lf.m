function Lf = Lf(x, param)
    Pf = param.Pf;
    Pf1 = Pf(1, 1);
    Pf2 = Pf(2, 2);
    Pf34 = Pf(3:4, 3:4);
    x_star = param.x_star;
    x1 = x(1);
    x2 = x(2);
    x1_star = x_star(1);
    x2_star = x_star(2);
    x34 = [x(3); x(4)];
    x34_star = [x_star(3); x_star(4)];

    Lf = 1/2*(x34-x34_star).'*Pf34*(x34-x34_star) + Pf1*(1-cos(x1-x1_star)) + Pf2*(1-cos(x2-x2_star));
end