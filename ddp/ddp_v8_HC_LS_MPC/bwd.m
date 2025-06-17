function [k, K] = bwd(X, U, P, param)
    H = param.H;

    Vx = NaN(1, 4, H+1);
    Vxx = NaN(4, 4, H+1);
    k = NaN(1, H);
    K = NaN(1,4, H);

    x_end = X(:, H+1);
    u_end = NaN;
    P_end = P(:, :, H+1);
    [Lx, ~, Lxx, ~] = der_L(x_end, u_end, P_end, param);
    Vx(1, :, H+1) = Lx;
    Vxx(:, :, H+1) = Lxx;

    for i = H:-1:1
        x_i = X(:, i);
        u_i = U(1, i);
        P_i = P(:, :, i);
        Vx_ii = Vx(1, :, i+1);
        Vxx_ii = Vxx(:, :, i+1);

        [Fx, Fu, ~, ~, ~] = der_F(x_i, u_i, param);
        [Lx, Lu, Lxx, Luu] = der_L(x_i, u_i, P_i, param);
        
        Qx = Lx + Vx_ii*Fx;
        Qu = Lu + Vx_ii*Fu;
        Qxx = Lxx + Fx.'*Vxx_ii*Fx;
        Quu = Luu + Fu.'*Vxx_ii*Fu;
        Qux = Fu.'*Vxx_ii*Fx;

        [k_i, K_i] = gains(u_i, Qu, Quu, Qux, param);
        
        Vx_i = Qx - (K_i.'*Quu*k_i).';
        Vxx_i = Qxx - K_i.'*Quu*K_i;

        k(1, i) = k_i;
        K(1, :, i) = K_i;
        Vx(1, :, i) = Vx_i;
        Vxx(:, :, i) = Vxx_i;
    end
end