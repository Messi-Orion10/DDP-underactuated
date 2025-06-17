function u_bound = bound(u, param)
    lim = param.lim;
    constr = param.constr;
    
    lb = lim(1);
    ub = lim(2);

    u_bound = u;
    if constr == 1 || constr == 3
        u_bound = min(ub, max(lb, u));
    else 
        if constr == 2
            u_bound = (ub - lb)/2*tanh(2*u/(ub - lb)) + (ub + lb)/2;
        end
    end
end