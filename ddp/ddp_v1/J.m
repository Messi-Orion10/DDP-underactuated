function J = J(x, u, param)
    J = 0;
    N = param.N;

    for i=1:N-1
        J = J + L(x(:,i), u(:,i), param);
    end
    
    Lf = L(x(:,N), [0;0], param);
    J = J + Lf;
end