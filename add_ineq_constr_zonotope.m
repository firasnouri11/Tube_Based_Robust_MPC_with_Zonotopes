function [F_ineq,g_ineq] = add_ineq_constr_zonotope(x_init, Z, Xc_robust,Uc_robust,n_x,n_w, N)
% Construct inequality constraints matrix that satisfies Fineq * x <= gineq
%  
n_epsilon = size(Z.G,2);
X_init = x_init + Z;% 
% construct -Phi_x <= Ksi_x <= Phi_x into F_ineq * x_ineq = g_ineq
F_ineq = [-eye(n_x) zeros(n_x, n_w) -eye(n_x); -eye(n_x) zeros(n_x, n_w) eye(n_x)];
F_ineq = [zeros(2*n_x, n_x*(N+1)+N) F_ineq];

g_ineq = zeros(size(F_ineq,1),1);


% Construct |Gamma_1|*1 + |Beta_1| <= 1

% Introduce slack variables -Lambda_1 <= Gamma_1 <= Lambda_1
% 2 inequalities for every element of Gamma_1 : -Lambda_1 -Gamma_1 <= 0 
% and Gamma_1 - Lambda_1 <= 0
F_Gamma_1 = [-eye(n_x*(n_x+n_epsilon)); eye(n_x*(n_x+n_epsilon))];
F_Lambda_1 = [-eye(n_x*(n_x+n_epsilon)); -eye(n_x*(n_x+n_epsilon))];
F_Gamma_Lambda = [zeros(size(F_Gamma_1,1),n_x*(N+1)+N+n_x+n_w+n_x) F_Gamma_1 F_Lambda_1 zeros(size(F_Gamma_1,1),2*n_x)];
g_Gamma_Lambda = zeros(size(F_Gamma_Lambda,1),1);

% Introduce slack variables -Miu_1 <= Beta_1 <= Miu_1
% 2 inequalities for every element of Beta_1 : -Miu_1 -Beta_1 <= 0 
% and Beta_1 - Miu_1 <= 0
F_Beta_1 = [-eye(n_x); eye(n_x)];
F_Miu_1 = [-eye(n_x); -eye(n_x)];
F_Beta_Miu = [zeros(size(F_Beta_1,1),n_x*(N+1)+N+n_x+n_w+n_x+2*n_x*(n_x+n_epsilon)) F_Beta_1 F_Miu_1];
g_Beta_Miu = [zeros(size(F_Beta_Miu,1),1)];

% Construct inequality Lambda_1*1 + Miu_1 <= 1
F_row_Lambda_1 = [repmat([1 0], 1, n_x+n_epsilon) 0; 0 repmat([1 0], 1, n_x+n_epsilon)];
F_row_Miu_1 = [eye(n_x)];
F_row_Lambda_Miu = [zeros(size(F_row_Lambda_1,1),n_x*(N+1)+N+n_x+n_w+n_x+n_x*(n_x+n_epsilon)) F_row_Lambda_1 zeros(size(F_row_Lambda_1,1),n_x-1) F_row_Miu_1];
g_row_Lambda_Miu = ones(n_x,1);

F_ineq = [F_ineq zeros(size(F_ineq,1),2*(n_x*(n_x+n_epsilon)+n_x)); F_Gamma_Lambda; F_Beta_Miu; F_row_Lambda_Miu];
g_ineq = [g_ineq; g_Gamma_Lambda; g_Beta_Miu; g_row_Lambda_Miu];


% (18h) Construct |Gamma_2|*1 + |Gamma_3|*1 + |Beta_2| <= 1

% Introduce slack variables -Lambda_2 <= Gamma_2 <= Lambda_2
% 2 inequalities for every element of Gamma_2 : -Lambda_2 -Gamma_2 <= 0 
% and Gamma_2 - Lambda_2 <= 0
F_Gamma_2 = [-eye(n_epsilon^2); eye(n_epsilon^2)];
F_Lambda_2 = [-eye(n_epsilon^2); -eye(n_epsilon^2)];
%F_Gamma_Lambda_2 = [zeros(size(F_Gamma_2,1),n_x) F_Gamma_2 F_Lambda_2];
F_Gamma_Lambda_2 = [zeros(size(F_Gamma_2,1),n_epsilon) F_Gamma_2 F_Lambda_2];
F_ineq = blkdiag(F_ineq,F_Gamma_Lambda_2);
g_Gamma_Lambda_2 = zeros(size(F_Gamma_Lambda_2,1),1);
g_ineq = [g_ineq; g_Gamma_Lambda_2];

% Introduce slack variables -Lambda_3 <= Gamma_3 <= Lambda_3
% 2 inequalities for every element of Gamma_3 : -Lambda_3 -Gamma_3 <= 0 
% and Gamma_3 - Lambda_3 <= 0
F_Gamma_3 = [-eye(n_epsilon*n_w); eye(n_epsilon*n_w)];
F_Lambda_3 = [-eye(n_epsilon*n_w); -eye(n_epsilon*n_w)];
F_Gamma_Lambda_3 = [F_Gamma_3 F_Lambda_3];
F_ineq = blkdiag(F_ineq,F_Gamma_Lambda_3);
g_Gamma_Lambda_3 = zeros(size(F_Gamma_Lambda_3,1),1);
g_ineq = [g_ineq; g_Gamma_Lambda_3];

% Introduce slack variables -Miu_2 <= Beta_2 <= Miu_2
% 2 inequalities for every element of Beta_2 : -Miu_2 -Beta_2 <= 0 
% and Beta_2 - Miu_2 <= 0
F_Beta_2 = [-eye(n_epsilon); eye(n_epsilon)];
F_Miu_2 = [-eye(n_epsilon); -eye(n_epsilon)];
F_Beta_Miu_2 = [F_Beta_2 F_Miu_2];
F_ineq = blkdiag(F_ineq, F_Beta_Miu_2);
g_Beta_Miu_2 = [zeros(size(F_Beta_Miu_2,1),1)];
g_ineq = [g_ineq; g_Beta_Miu_2];

% Construct inequality Lambda_2*1 + Lambda_3*1 + Miu_2 - Phi_epsilon*1 <= 0
F_row_Lambda_2 = repmat(eye(n_epsilon), 1, n_epsilon);
F_row_Lambda_3 = repmat(eye(n_epsilon), 1, n_w);
F_row_Miu_2 = [eye(n_epsilon)];
F_row_Phi_epsilon = -eye(n_epsilon);
F_row_Lambda_Miu_2_3 = [zeros(size(F_row_Lambda_2,1),n_x*(N+1)+N+n_x+n_w+n_x+2*n_x*(n_x+n_epsilon)+2*n_x) F_row_Phi_epsilon zeros(size(F_row_Lambda_2,1),n_epsilon^2) F_row_Lambda_2 zeros(size(F_row_Lambda_2,1),n_epsilon*n_w) F_row_Lambda_3 zeros(size(F_row_Lambda_2,1),n_epsilon) F_row_Miu_2];
g_row_Lambda_Miu_2_3 = zeros(n_epsilon,1);

F_ineq = [F_ineq; F_row_Lambda_Miu_2_3];
g_ineq = [g_ineq; g_row_Lambda_Miu_2_3];
%{
% Construct -Phi_epsilon <= Ksi_epsilon <= Phi_epsilon 
% 2 inequalities for every element of Ksi_epsilon : -Phi_epsilon -Ksi_epsilon <= 0
% and Ksi_epsilon - Phi_epsilon <= 0
F_Ksi_epsilon = [-eye(n_epsilon); eye(n_epsilon)];
F_Phi_epsilon = [-eye(n_epsilon); -eye(n_epsilon)];
F_Phi_Ksi_epsilon = [zeros(size(F_Phi_epsilon,1),n_x*(N+1)+N+n_x+n_w+n_x+2*(n_x*(n_x+n_epsilon)+n_x)) F_Phi_epsilon zeros(size(F_Ksi_epsilon,1),2*(n_epsilon^2+n_epsilon*n_w+n_epsilon)) F_Ksi_epsilon];
g_Phi_Ksi_epsilon = zeros(size(F_Phi_Ksi_epsilon,1),1);
F_ineq = [F_ineq zeros(size(F_ineq,1),n_epsilon); F_Phi_Ksi_epsilon];
g_ineq = [g_ineq; g_Phi_Ksi_epsilon];
%}
end