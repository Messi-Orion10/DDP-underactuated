% function [k, K] = bwd(X, U, param)
%     N = param.N;
%     P = param.P;
%     R = param.R;
%     Pf = param.Pf;
%     x_star = param.x_star;
%     mu = param.mu;
% 
%     Qx = NaN(1,4,N);
%     Qu = NaN(1,1,N);
%     Qxx = NaN(4,4,N);
%     Quu = NaN(1,1,N);
%     Quu_tilde = NaN(1,1,N);
%     Qux = NaN(1,4,N);
%     Vx = NaN(1,4,N+1);
%     Vxx = NaN(4,4,N+1);
% 
%     k = NaN(1,N);
%     K = NaN(1,4,N);
% 
%     Vx(1, :, N+1) = (X(:,N+1)-x_star).'*Pf;
%     Vxx(:, :, N+1) = Pf;
% 
%     for i = N:-1:1
%         [Fx, Fu, Fxx, Fux, Fxu] = F_der(X(:,i), U(1,i), param);  
%         % for j = 1:4
%         %     Vx_prime_Fxx = squeeze(Vx(1, j, i+1)*Fxx(j, :, :));
%         %     Vx_prime_Fux = squeeze(Vx(1, j, i+1)*Fux(j, :, :));
%         %     Vx_prime_Fxu = squeeze(Vx(1, j, i+1)*Fxu(j, :, :));
%         % end
%         Qx(1,:,i) = (X(:,i)-x_star).'*P + Vx(1, :, i+1)*Fx;
%         Qu(1,1,i) = U(1,i).'*R + Vx(1, :, i+1)*Fu;
%         Qxx(:,:,i) = P + Fx.'*Vxx(:, :, i+1)*Fx;% + Vx_prime_Fxx;
%         Quu(1,1,i) = R + Fu.'*Vxx(:, :, i+1)*Fu;
%         Quu_tilde(1,1,i) = Quu(1,1,i) + mu;
%         Qux(1,:,i) = Fu.'*Vxx(:, :, i+1)*Fx;% + Vx_prime_Fux;
%         Qxu(:,1,i) = Fx.'*Vxx(:, :, i+1)*Fu;% + Vx_prime_Fxu;
% 
%         k(1,i) = -squeeze(Quu_tilde(1,1,i))\(Qu(1,1,i).');
%         K(1,:,i) = -squeeze(Quu_tilde(1,1,i))\(squeeze(Qux(1,:,i)));
% 
%         Vx(1,:,i) = Qx(1,:,i) - (K(1,:,i).'*Quu(1,1,i)*k(1,i)).';
%         Vxx(:,:,i) = Qxx(:,:,i) - K(1,:,i).'*Quu(1,1,i)*K(1,:,i);
% 
%         %Vxx(:,:,i) = Vxx(:,:,i) + mu*eye(4);
%     end
% end


function [k, K] = bwd(X, U, param)
    % Extract parameters
    N = param.N;
    P = param.P;
    R = param.R;
    Tf = param.Tf;
    Pf = param.Pf;
    x_star = param.x_star;
    mu = param.mu;

    % Initialize matrices
    Qx = NaN(1,4,N);
    Qu = NaN(1,1,N);
    Qxx = NaN(4,4,N);
    Quu = NaN(1,1,N);
    Quu_tilde = NaN(1,1,N);
    Qux = NaN(1,4,N);
    Qxu = NaN(4,1,N);  % Added Qxu initialization
    Vx = NaN(1,4,N+1);
    Vxx = NaN(4,4,N+1);

    k = NaN(1,N);
    K = NaN(1,4,N);

    P_seq = NaN(4, 4, N+1);
    P_seq(:, :, 1:N+1-Tf) = repmat(P, 1, 1, N+1-Tf);
    P_seq(:, :, N+1-Tf+1:end) = repmat(Pf, 1, 1, Tf);

    % Terminal cost
    Vx(1, :, N+1) = (X(:,N+1) - x_star).' * P_seq(:, :, N+1);
    Vxx(:, :, N+1) = P_seq(:, :, N+1);

    % Backward pass
    for i = N:-1:1
        % Compute derivatives
        [Fx, Fu, Fxx, Fux, Fxu] = F_der(X(:,i), U(1,i), param);  
        
        % Value function derivatives
        Qx(1,:,i) = (X(:,i) - x_star).' * P_seq(:, :, i) + Vx(1, :, i+1) * Fx;
        Qu(1,1,i) = U(1,i).' * R + Vx(1, :, i+1) * Fu;
        Qxx(:,:,i) = P_seq(:, :, i) + Fx.' * Vxx(:, :, i+1) * Fx;  % + Vx_prime_Fxx;
        Quu(1,1,i) = R + Fu.' * Vxx(:, :, i+1) * Fu;
        Quu_tilde(1,1,i) = Quu(1,1,i) + mu;
        Qux(1,:,i) = Fu.' * Vxx(:, :, i+1) * Fx;  % + Vx_prime_Fux;
        Qxu(:,1,i) = Fx.' * Vxx(:, :, i+1) * Fu;  % + Vx_prime_Fxu;

        % Ensure Quu_tilde is invertible
        assert(abs(det(Quu_tilde(1,1,i))) > 1e-9, 'Quu_tilde is nearly singular');

        % Control updates
        k(1,i) = -squeeze(Quu_tilde(1,1,i)) \ (Qu(1,1,i).');
        K(1,:,i) = -squeeze(Quu_tilde(1,1,i)) \ (squeeze(Qux(1,:,i)));

        % Update value function
        Vx(1,:,i) = Qx(1,:,i) - (K(1,:,i).' * Quu(1,1,i) * k(1,i)).';
        Vxx(:,:,i) = Qxx(:,:,i) - K(1,:,i).' * Quu(1,1,i) * K(1,:,i);
    end
end