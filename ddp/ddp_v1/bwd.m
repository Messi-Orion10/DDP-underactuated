function [k, K] = bwd(X, U, param)
    N = param.N;
    P = param.P;
    R = param.R;
    Pf = param.Pf;
    x_star = param.x_star;
    mu = param.mu;

    Qx = NaN(1,4,N);
    Qu = NaN(1,2,N);
    Qxx = NaN(4,4,N);
    Quu = NaN(2,2,N);
    Quu_tilde = NaN(2,2,N);
    Qux = NaN(2,4,N);
    Vx = NaN(1,4,N+1);
    Vxx = NaN(4,4,N+1);

    k = NaN(2,N);
    K = NaN(2,4,N);

    Vx(1, :, N+1) = (X(:,N+1)-x_star).'*Pf;
    Vxx(:, :, N+1) = Pf;

    for i = N:-1:1
        [Fx, Fu, Fxx, Fux, Fxu] = F_der(X(:,i), U(:,i), param);  
        % for j = 1:4
        %     Vx_prime_Fxx = squeeze(Vx(1, j, i+1)*Fxx(j, :, :));
        %     Vx_prime_Fux = squeeze(Vx(1, j, i+1)*Fux(j, :, :));
        %     Vx_prime_Fxu = squeeze(Vx(1, j, i+1)*Fxu(j, :, :));
        % end
        Qx(1,:,i) = (X(:,i)-x_star).'*P + Vx(1, :, i+1)*Fx;
        Qu(1,:,i) = U(:,i).'*R + Vx(1, :, i+1)*Fu;
        Qxx(:,:,i) = P + Fx.'*Vxx(:, :, i+1)*Fx;% + Vx_prime_Fxx;
        Quu(:,:,i) = R + Fu.'*Vxx(:, :, i+1)*Fu;
        Quu_tilde(:,:,i) = Quu(:,:,i) + mu*eye(2);
        Qux(:,:,i) = Fu.'*Vxx(:, :, i+1)*Fx;% + Vx_prime_Fux;
        Qxu(:,:,i) = Fx.'*Vxx(:, :, i+1)*Fu;% + Vx_prime_Fxu;

        k(:,i) = -squeeze(Quu_tilde(:,:,i))\(Qu(1,:,i).');
        K(:,:,i) = -squeeze(Quu_tilde(:,:,i))\(squeeze(Qux(:,:,i)));

        Vx(1,:,i) = Qx(1,:,i) - (K(:,:,i).'*Quu(:,:,i)*k(:,i)).';
        Vxx(:,:,i) = Qxx(:,:,i) - K(:,:,i).'*Quu(:,:,i)*K(:,:,i);

        %Vxx(:,:,i) = Vxx(:,:,i) + mu*eye(4);
    end
end