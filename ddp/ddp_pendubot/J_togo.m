function J_togo = J_togo(x, u, i, param)
    J_togo = 0;
    N = param.N;

    for j=i:N-1
        J_togo = J_togo + L(x(:,j), u(1,j), param);
    end
    
    J_togo = J_togo + Lf(x(:,N), param);
end