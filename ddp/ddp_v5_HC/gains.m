function [k, K] = gains(u, Qu, Quu, Qux, param)
    lim = param.lim;
    mu = param.mu;
    constr = param.constr;

    Quu_tilde = Quu + mu;

    k = -Qu/Quu_tilde;
    K = -Qux/Quu_tilde;
    if constr == 3
        lb = lim(1);
        ub = lim(2);

        [k, clamped] = boxQP(Quu_tilde, Qu, lb - u, ub - u);
        
        if clamped
            K = zeros(1, 4);
        end
    end
end


        

