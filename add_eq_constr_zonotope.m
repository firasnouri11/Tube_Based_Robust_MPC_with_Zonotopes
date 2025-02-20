function [Feq,geq] = add_eq_constr_zonotope(x_init, A,B,K,N,Xc,Xc_robust,W,Z)
% 
A_cl = A - B*K;
n_w = size(W.G,2);
%Construct equality constraints matrix that satisfies Feq * x = geq
nx = size(A,1);
Xeq = []; Ueq = [];
geq = zeros(nx*N,1);
for i=1:N
    Xeq = blkdiag(Xeq,-A);
    Ueq = blkdiag(Ueq,-B);
end
Xeq = [Xeq zeros(size(Xeq,1),nx)];
row = 1;
col = size(A,2) +1;
for j=1:size(Xeq,1)
    for k= (size(A,2)+1):size(Xeq,2)
        if (row == j && col == k)
            Xeq(j,k) = 1;
            row = row + 1;
            col = col +1;
        end
    end
end
Feq = [Xeq Ueq];

F_X_predicted = eye(nx);
F_G_x = -Xc_robust.G;
F = [F_X_predicted zeros(size(F_X_predicted,1),(nx*N+N+nx+n_w)) F_G_x];
Feq = [Feq zeros(size(Feq,1), size(F,2)-size(Feq,2)); F];
g = [Xc.c];
%{
for i=1:N
    g = [g; Xc.c];
end
%}
geq = [geq;g];

% construct [Gx_tilde*Phi_x G_epsilon*Phi_epsilon] - G_x*Gamma_1 = 0
n_x_tilde = size(Xc_robust.G,2);
n_epsilon = size(Z.G,2); 
F_G_x_tilde = [zeros(n_x_tilde*nx, nx*(N+1)+N) blkdiag(Xc_robust.G(1:nx,1), Xc_robust.G(1:nx,2))];
F_G_x_tilde = [F_G_x_tilde zeros(size(F_G_x_tilde,1),nx+n_w) -blkdiag(Xc_robust.G,Xc_robust.G)];
Feq = [Feq zeros(size(Feq,1),nx+nx); F_G_x_tilde];
geq = [geq; zeros(size(F_G_x_tilde,1),1)];
F_Gamma_1 = -Xc.G;
F_G_epsilon = Z.G(1:nx,1);
for i=1:(n_epsilon-1)
    F_Gamma_1 = blkdiag(F_Gamma_1,-Xc.G);
    F_G_epsilon = blkdiag(F_G_epsilon, Z.G(1:nx,i+1));
end
F_Gamma_1 = [F_Gamma_1 zeros(n_epsilon*nx,nx*(n_x_tilde+n_epsilon)+nx+nx) F_G_epsilon];
Feq = blkdiag(Feq, F_Gamma_1);
geq = [geq; zeros(size(F_Gamma_1,1),1)];

% Construct c_x - (c_x_tilde + c_epsilon) = G_x * Beta_1
F_Beta_1 = [zeros(nx,nx*(N+1)+N+nx+n_w+nx+2*nx*(n_x_tilde+n_epsilon)) Xc.G zeros(nx, nx+n_epsilon)];
g_Beta_1 = Xc.c - (Xc_robust.c + Z.c);
Feq = [Feq; F_Beta_1];
geq = [geq; g_Beta_1];

% Construct A_cl * G_epsilon * Phi_epsilon - G_epsilon * Gamma_2 = 0
A_cl_G_epsilon = A_cl * Z.G;
F_A_cl_G_epsilon = zeros(nx*n_epsilon, n_epsilon);
F_Gamma_2 = zeros(nx*n_epsilon, n_epsilon^nx);
for i = 1:n_epsilon
    F_A_cl_G_epsilon((2*i-1):2*i, i) = A_cl_G_epsilon(:,i);
    F_Gamma_2((2*i-1):2*i, ((i-1)*n_epsilon+1):(i*n_epsilon)) = -Z.G;
end
F_18_e = [F_A_cl_G_epsilon F_Gamma_2];
Feq = [Feq zeros(size(Feq,1),n_epsilon*n_epsilon); zeros(size(F_18_e,1),size(Feq,2)-n_epsilon) F_18_e];
geq = [geq; zeros(size(F_18_e,1),1)];

% (18f) Construct G_w = G_epsilon * Gamma_3
F_Gamma_3 = [Z.G];
for i = 2:n_w
    F_Gamma_3 = blkdiag(Z.G,Z.G);
end
F_Gamma_3 = [zeros(size(F_Gamma_3,1),size(F_Gamma_2,2)) F_Gamma_3]; 
Feq = blkdiag(Feq,F_Gamma_3);
geq = [geq; reshape(W.G,[],1)];

% (18g) Construct (I-A_cl)*c_epsilon -c_w = G_epsilon * Beta_2
I_A_cl = (eye(nx) - A_cl)*Z.c - W.c;
F_Beta_2 = [zeros(size(Z.G,1),n_w*n_epsilon) Z.G];
Feq = blkdiag(Feq,F_Beta_2);
geq = [geq; I_A_cl];

Feq = [Feq zeros(size(Feq,1),n_epsilon)];

end
