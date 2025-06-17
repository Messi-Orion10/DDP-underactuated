function F = dyn_dist(x, u, d, param)
    T = param.T;
    act = param.act;
    a1 = param.a1;
    a2 = param.a2;
    a3 = param.a3;
    a4 = param.a4;
    a5 = param.a5;
    b1 = param.b1;
    b2 = param.b2;

    x1 = x(1);
    x2 = x(2);
    x3 = x(3);
    x4 = x(4);

    M = [a1+2*a2*cos(x2) a3+a2*cos(x2); 
        a3+a2*cos(x2) a3];
    b = diag([b1;b2]);
    c = a2*sin(x2)*[x4+2*x3; x3^2];
    e = [a4*sin(x1)+a5*sin(x1+x2); 
        a5*sin(x1+x2)];

    if(act == 2)
        u_bar = [0; u];
    else
        u_bar = [u; 0];
    end

    f12 = [x3; x4];
    f34 = M\(-b*f12 - c - e + u_bar + d);
    f = [f12; f34];

    F = x + T*f;
end