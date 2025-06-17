function J = cost_tot(X, U, param)
    H = param.H;

    J = 0;
    for i = 1:H
        J = J + cost(X(:, i), U(1, i), i, param);
    end
    J = J + cost(X(:, H+1), 0, H+1, param);
end