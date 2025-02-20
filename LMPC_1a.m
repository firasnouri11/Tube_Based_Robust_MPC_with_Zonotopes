% MATLAB program for Linear MPC: Single-input system (alternate code for LMPC_1.m with constraints defined using lb and ub)
clear all;
close all

% System parameters and simulation parameters
A = [1 1; 0 1];
B = [0; 1];

NT=100;N=5;n=2;m=1; 
Q=eye(n); QN=Q; R=0*eye(m); [K,P] = dlqr(A,B,Q,R);
A_cl = A - B .* K; 
x_min=[-1;-5];x_max=[60;5]; 
umin=-1;umax=1;
x_min = repmat(x_min,1,N+1); x_max = repmat(x_max,1,N+1);
Xc = zonotope([29.5;0], [30.5 0;0 5]);
Uc = zonotope(polytope([-1 1]));

% define start point
x0=[53;-5]; 
x=zeros(n,NT+1); x(:,1)=x0;
Xk=zeros(n*(N+1),1); Xk(1:n,1)=x0;
u=zeros(m,NT);
Uk=zeros(m*N,1);
zk=[Xk;Uk]; 

x_ref = zeros(n,(N+NT+1));

% define disturbance set W and calculate bounding set Epsilon
w = zeros(n,NT); 
W = zonotope([0; 0], [0.5 0;0 0.5]); % convex set of disturbance
n_w = size(W.G,2);
lambda = 10;
s = 5;

[Z, Xc_robust, Uc_robust] = compute_disturbance_invariance_set(A,B,K,W,s,Xc,Uc);

n_epsilon = size(Z.G,2);

% constructing H and f to comply with 0.5x.*H*x + f.*x cost function of quadprog 
for i=1:N+1
    AX((i-1)*n+1:i*n,:)=A^(i-1);
end
for i=1:N+1
  for j=1:N
      if i>j
          BU((i-1)*n+1:i*n,(j-1)*m+1:j*m)=A^(i-j-1)*B;
      else
          BU((i-1)*n+1:i*n,(j-1)*m+1:j*m)=zeros(n,m);
      end    
  end
end

QX=Q;RU=R;
for i=1:N-1
  QX=blkdiag(QX,Q); RU=blkdiag(RU,R);
end
QX=blkdiag(QX,P);
Weight_matrix=blkdiag(QX,RU);

Phi_x = zeros(n,NT);
Phi_epsilon = zeros(n_epsilon,NT);
computation_time = zeros(1,NT);

for k=1:NT
   % define current state x
   xk=x(:,k);
   
   [Feq,geq] = add_eq_constr_zonotope(xk,A,B,K,N,Xc,Xc_robust,W,Z);
   [Fineq,gineq] = add_ineq_constr_zonotope(xk, Z, Xc_robust,Uc_robust,n,n_w, N);
   
   scaled_x_ref = -1 * [reshape(x_ref(:,k:(k+N)),[],1); zeros(N,1)];
   f = Weight_matrix* scaled_x_ref;
   f((n*N+1):(n*N+2)) = [0;0];
   
   % only use this H matrix & f vector if we only want set containment (No MPC)
   Weight_matrix = zeros(size(Weight_matrix));
   f = zeros(size(f,1),1);
   [H, f] = construct_cost_function(Weight_matrix,f,lambda,n,n_w, n_epsilon);
   
   zk = [xk; zeros(size(H,2)-size(xk,1),1)];
   
   options = optimoptions('quadprog', 'Algorithm', 'active-set', 'MaxIter', 2000, 'Display', 'iter');
   tic
   [z1,fval]=quadprog(H,f,Fineq,gineq,Feq,geq,[],[],zk,options);
   computation_time(k) = toc;
   
   k_w = n*(N+1)+N+n; %+n_w+n+2*n*(n+n_epsilon)+2*n;
   k_epsilon = n*(N+1)+N+n+n_w+n+2*n*(n+n_epsilon)+2*n;
   Phi_x(:,k) = z1(n*(N+1)+N+1:n*(N+1)+N+n);
   Phi_epsilon(:,k) = z1((k_epsilon+1):(k_epsilon+n_epsilon));
   Phi_w = diag(z1((k_w+1):(k_w+n_w)));

   zk = z1;
end

average_Phi_x = diag(mean(Phi_x,2));
average_Phi_epsilon = diag(mean(Phi_epsilon,2));
Z_scaled = zonotope(Z.c,Z.G*average_Phi_epsilon);

figure(1)
plot(Xc, [1 2],'k','FaceColor',[0.9290 0.6940 0.1250]);
hold on
plot(Xc_robust,[1 2],'k','FaceColor',[0 0.4470 0.7410])
plot(minkDiff(Xc,Z_scaled), [1 2],'k','FaceColor',[0.4940 0.1840 0.5560])
plot(zonotope(Xc_robust.c,Xc_robust.G*average_Phi_x), [1 2],'k','FaceColor',[0.8500 0.3250 0.0980]); % plotting scaled Xc set
plot(interval([1.04;-4],[1.2;-3]), [1 2], 'EdgeColor',[0.5 0.5 0.5], 'LineWidth',1);
xlabel('$\textbf{x}_{1}$','Interpreter','latex','FontSize',14);ylabel('$\textbf{x}_{2}$','Interpreter','latex','FontSize',14);
legend('$X$','$\hat{X}$','$X \ominus \tilde{\varepsilon}(\Phi_{\varepsilon})$','$\tilde{X}(\Phi_x)$','Interpreter','latex','FontSize',14);
title('$\tilde{X}(\Phi_x) \subseteq X \ominus \tilde{\varepsilon}(\Phi_{\varepsilon}) \subseteq \hat{X} = X \ominus \varepsilon$','Interpreter','latex', 'FontSize',16)
ylim([x_min(2)-2 x_max(2)+2])
xlim([2*x_min(1) x_max(1)+2])

figure(2)
plot(zonotope(Xc_robust.c,Xc_robust.G*average_Phi_x),[1 2],'k','FaceColor',[0 0.4470 0.7410]);
xlabel('$\textbf{x}_{1}$','Interpreter','latex','FontSize',14);ylabel('$\textbf{x}_{2}$','Interpreter','latex','FontSize',14);
legend('$\tilde{X}(\Phi_x)$','Interpreter','latex','FontSize',14);
xlim([2*x_min(1) x_max(1)+2])
ylim([-2.1 2.1])

figure(3)
plot(Z_scaled,[1 2],'k','FaceColor',[0.8500 0.3250 0.0980]); % plotting scaled Xc set
hold on
plot(Z,[1 2],'k','FaceColor',[0 0.4470 0.7410]);
xlabel('$\textbf{x}_{1}$','Interpreter','latex','FontSize',14);ylabel('$\textbf{x}_{2}$','Interpreter','latex','FontSize',14);
legend('$\tilde{\varepsilon}(\Phi_{\varepsilon})$','$\varepsilon$','Interpreter','latex','FontSize',14);
grid on

figure(4)
plot(Xc, [1 2],'k','FaceColor',[0.9290 0.6940 0.1250]);
hold on
plot(Xc_robust,[1 2],'k','FaceColor',[0 0.4470 0.7410])
plot(minkDiff(Xc,Z_scaled), [1 2],'k','FaceColor',[0.4940 0.1840 0.5560])
xlabel('$\textbf{x}_{1}$','Interpreter','latex','FontSize',14);ylabel('$\textbf{x}_{2}$','Interpreter','latex','FontSize',14);
legend('$X$','$\hat{X}$','$X \ominus \tilde{\varepsilon}(\Phi_{\varepsilon})$','Interpreter','latex','FontSize',14);
title('$\tilde{X}(\Phi_x) \subseteq X \ominus \tilde{\varepsilon}(\Phi_{\varepsilon}) \subseteq \hat{X} = X \ominus \varepsilon$','Interpreter','latex', 'FontSize',16)
xlim([1.04 1.2])
ylim([-4 -3])
