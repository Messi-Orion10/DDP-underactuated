function [k, K] = bwd(X, U, param)
    H = param.H;
    mu = param.mu;

    Vx = NaN(1, 4, H+1);
    Vxx = NaN(4, 4, H+1);
    k = NaN(1, H);
    K = NaN(1,4, H);

    x_fin = X(:, H+1);
    [Lx, ~, Lxx, ~] = der_L(x_fin, NaN, H+1, param);
    Vx(1, :, H+1) = Lx;
    Vxx(:, :, H+1) = Lxx;

    % Backward pass
    for i = H:-1:1
        x_i = X(:, i);
        u_i = U(1, i);
        Vx_ii = Vx(1, :, i+1);
        Vxx_ii = Vxx(:, :, i+1);

        [Fx, Fu, ~, ~, ~] = der_F(x_i, u_i, param);
        [Lx, Lu, Lxx, Luu] = der_L(x_i, u_i, i, param);
        
        Qx = Lx + Vx_ii*Fx;
        Qu = Lu + Vx_ii*Fu;
        Qxx = Lxx + Fx.'*Vxx_ii*Fx;
        Quu = Luu + Fu.'*Vxx_ii*Fu;
        Quu_tilde = Quu + mu;
        Qux = Fu.'*Vxx_ii*Fx;

        assert(abs(det(Quu_tilde)) > 1e-9, 'Quu_tilde is nearly singular');

        k_i = -squeeze(Quu_tilde)\(Qu.');
        K_i = -squeeze(Quu_tilde)\(squeeze(Qux));
        Vx_i = Qx - (K_i.'*Quu*k_i).';
        Vxx_i = Qxx - K_i.'*Quu*K_i;

        k(1, i) = k_i;
        K(1, :, i) = K_i;
        Vx(1, :, i) = Vx_i;
        Vxx(:, :, i) = Vxx_i;
    end
end