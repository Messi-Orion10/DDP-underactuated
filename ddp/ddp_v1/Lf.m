function Lf = Lf(x, param)
    Pf = param.Pf;
    x_star = param.x_star;

    Lf = 1/2*(x-x_star).'*Pf*(x-x_star);
end