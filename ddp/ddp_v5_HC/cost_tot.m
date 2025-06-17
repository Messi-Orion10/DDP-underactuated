function J = cost_tot(X, U, param)
    N = param.N;

    J = 0;
    for i = 1:N
        J = J + cost(X(:, i), U(1, i), i, param);
    end
    J = J + cost(X(:, N+1), 0, N+1, param);
end