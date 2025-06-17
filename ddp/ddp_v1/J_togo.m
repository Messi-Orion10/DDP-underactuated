function J_togo = J_togo(x, u, i, param)
    J_togo = 0;
    N = param.N;

    for j=i:N-1
        J_togo = J_togo + L(x(:,j), u(:,j), param);
    end
    
    Lf = L(x(:,N), [0;0], param);
    J_togo = J_togo + Lf;
end