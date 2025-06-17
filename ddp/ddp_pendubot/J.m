function J = J(x, u, param)
    J = 0;
    N = param.N;

    for i=1:N-1
        J = J + L(x(:,i), u(1,i), param);
    end
    
    J = J + Lf(x(:,N), param);
end