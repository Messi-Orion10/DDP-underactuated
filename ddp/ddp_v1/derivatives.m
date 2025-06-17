clc, clear, close all
%% derivatives
syms a [5,1];
syms x [4,1];
syms u [2,1];
syms r [2,1];
syms p [4,1];
syms T b1 b2;
M = [a1+2*a2*cos(x2) a3+a2*cos(x2); 
        a3+a2*cos(x2) a3];
b = diag([b1;b2]);
c = a2*sin(x2)*[x4+2*x3; x3^2];
e = [a4*sin(x1)+a5*sin(x1+x2); 
    a5*sin(x1+x2)];
f = [x3;x4;0;0];
f(3:4) = M\(-b*[x3;x4]-c-e+u);
F = x + T*f;
for i = 1:4
    Fx(:,i) = diff(F,x(i));
end
Fx = simplify(Fx);
for i = 1:2
    Fu(:,i) = diff(F,u(i));
end
Fu = simplify(Fu);
for i = 1:4
    Fxx(:,:,i) = diff(Fx,x(i));
    Fux(:, :, i) = diff(Fu, x(i));
end
Fux = simplify(Fux);
Fxx = simplify(Fxx);
for i = 1:2
    Fuu(:, :, i) = diff(Fu, u(i));
    Fxu(:, :, i) = diff(Fx, u(i));
end
Fuu = simplify(Fuu);
Fxu = simplify(Fxu);

%% Apri un file di testo per scrittura
fileID = fopen('tensor_derivate.txt', 'w');
tens = Fxu;
% Scrivi il tensore nel file
for i = 1:size(tens, 1)
    for j = 1:size(tens, 2)
        for k = 1:size(tens,3)
            fprintf(fileID, 'Fxu(%d,%d,%d) = %s;\n', i, j, k, char(tens(i,j,k)));
        end
    end
end

% Chiudi il file
fclose(fileID);
