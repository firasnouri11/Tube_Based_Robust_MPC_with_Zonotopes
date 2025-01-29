function [H_new,f_new] = construct_cost_function(H,f,lambda,n_x,n_w, n_epsilon)

% build the H matrix and f vector to accomodate the scaling matrix
% variables
H_scaling_matrix = -1 * blkdiag(eye(n_x), lambda * eye(n_w)); % H matrix for quadprog;
f_scaling_matrix = zeros(size(H_scaling_matrix,1),1);
H_new = blkdiag(H,H_scaling_matrix);
f_new = [f;f_scaling_matrix];

% add zeros for other variables that are not optimized
H_Ksi_x = zeros(n_x, n_x);
H_Gamma_1 = zeros(n_x * (n_x + n_epsilon), n_x * (n_x + n_epsilon));
H_Lambda_1 = zeros(n_x * (n_x + n_epsilon), n_x * (n_x + n_epsilon));
H_Beta_1 = zeros(n_x, n_x);
H_Miu_1 = zeros(n_x, n_x);
H_Phi_epsilon = zeros(n_epsilon, n_epsilon);
H_Gamma_2 = zeros(n_epsilon^2, n_epsilon^2);
H_Lambda_2 = zeros(n_epsilon^2, n_epsilon^2);
H_Gamma_3 = zeros(n_epsilon * n_w, n_epsilon * n_w);
H_Lambda_3 = zeros(n_epsilon * n_w, n_epsilon * n_w);
H_Beta_2 = zeros(n_epsilon, n_epsilon);
H_Miu_2 = zeros(n_epsilon, n_epsilon);
H_Ksi_epsilon = zeros(n_epsilon, n_epsilon);
H_c_epsilon = zeros(n_x,n_x);

f_Ksi_x = zeros(n_x, 1);
f_Gamma_1 = zeros(n_x * (n_x + n_epsilon), 1);
f_Lambda_1 = zeros(n_x * (n_x + n_epsilon), 1);
f_Beta_1 = zeros(n_x, 1);
f_Miu_1 = zeros(n_x, 1);
f_Phi_epsilon = zeros(n_epsilon, 1);
f_Gamma_2 = zeros(n_epsilon^2, 1);
f_Lambda_2 = zeros(n_epsilon^2, 1);
f_Gamma_3 = zeros(n_epsilon * n_w, 1);
f_Lambda_3 = zeros(n_epsilon * n_w, 1);
f_Beta_2 = zeros(n_epsilon,1);
f_Miu_2 = zeros(n_epsilon,1);
f_Ksi_epsilon = zeros(n_epsilon,1);
f_c_epsilon = zeros(n_x,1);

% construct all matrices
H_new = blkdiag(H_new, H_Ksi_x, H_Gamma_1, H_Lambda_1, H_Beta_1, H_Miu_1, H_Phi_epsilon, H_Gamma_2, H_Lambda_2, H_Gamma_3, H_Lambda_3, H_Beta_2, H_Miu_2 );
f_new = [f_new; f_Ksi_x; f_Gamma_1; f_Lambda_1; f_Beta_1; f_Miu_1; f_Phi_epsilon; f_Gamma_2; f_Lambda_2; f_Gamma_3; f_Lambda_3; f_Beta_2; f_Miu_2 ];
end
