% MATLAB program for Linear MPC: Single-input system (alternate code for LMPC_1.m with constraints defined using lb and ub)
clear all;
close all

% System parameters and simulation parameters
A = [1 1; 0 1];
B = [0; 1];

NT=100;N=15;n=2;m=1; 
Q=eye(n); QN=Q; R=0*eye(m); [K,P] = dlqr(A,B,Q,R); % at first we don't need K if we don't have any noise
A_cl = A - B .* K; 
x_min=[-1;-5];x_max=[60;5]; 
umin=-1;umax=1;
x_min = repmat(x_min,1,N+1); x_max = repmat(x_max,1,N+1);
Xc = zonotope([29.5;0], [30.5 0;0 5]);
Uc = zonotope(polytope([-1 1]));

% define start point
x0=[13;-3]; 
x=zeros(n,NT+1); x(:,1)=x0;
Xk=zeros(n*(N+1),1); Xk(1:n,1)=x0;
u=zeros(m,NT);
Uk=zeros(m*N,1);
zk=[Xk;Uk]; 

% define reference trajectory to be tracked
x_ref = zeros(n,(N+NT+1)); x_ref(1,10:40) = 40; x_ref(1,65:85) = 10; x_ref(2,10:80) = 0;
u_ref = zeros(1,(N+NT));u_ref(:,10:80) = pinv(B)*(eye(n)-A)*x_ref(:,10:80);

% define disturbance set W and calculate bounding set Epsilon
w = zeros(n,NT); 
W = zonotope([0; 0], [0.1 0;0  0.1]); % convex set of disturbance
n_w = size(W.G,2);
lambda = 10;
%W = zonotope([0; 0], [ 0; 1]);
[Z, Xc_robust, Uc_robust] = compute_disturbance_invariance_set(A,B,K,W,N,Xc,Uc);
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
QX=blkdiag(QX,QN);
Weight_matrix=blkdiag(QX,RU);
%{
lb=[reshape(x_min,1,[]) umin*ones(1,N*m)]; % set lower bound for state x and input u
ub=[reshape(x_max,1,[]) umax*ones(1,N*m)];
%}
x_predicted = zeros(n*(N+1),NT+1);
u_predicted = zeros(m*N,NT);

% simulating system with MPC
for k=1:NT
   % define current state x
   xk=x(:,k);
   
   [Feq,geq] = add_eq_constr_zonotope(xk,A,B,K,N,Xc,Xc_robust,W,Z);
   [Fineq,gineq] = add_ineq_constr_zonotope(xk, Z, Xc_robust,Uc_robust,n,n_w, N);

   %{ 
   this was part of fmincon solver, no need to uncomment
   z_ref_iter = [reshape(x_ref(:,k:(N+k)),[],1); u_ref(k:(N+k-1))'];
   % define cost function as a reference tracking problem (input is not weighted) 
   fun = @(z)(z_ref_iter - z)'*H*(z_ref_iter -z);
   
   F=[];g=[];Feq=[eye((N+1)*n) -BU];geq=AX*xk;
   %z1=fmincon(fun,zk,F,g,Feq,geq,lb,ub);
   %}
   
   scaled_x_ref = -1 * [reshape(x_ref(:,k:(k+N)),[],1); zeros(N,1)];
   f = Weight_matrix* scaled_x_ref;
   f((n*N+1):(n*N+2)) = [0;0];
   [H, f] = construct_cost_function(Weight_matrix,f,lambda,n,n_w, n_epsilon);
   zk = [xk; zeros(size(H,2)-size(xk,1),1)];
   
   options = optimoptions('quadprog', 'Algorithm', 'active-set',  'Display', 'iter');
   %options = optimoptions('quadprog', 'Display', 'iter');
   tStart = cputime;
   z1=quadprog(H,f,Fineq,gineq,Feq,geq,[],[],zk,options);
   tEnd = cputime - tStart

   x_predicted(:,k) = z1(1:n*(N+1),1);
   u_predicted(:,k) = z1((N+1)*n+1:(N+1)*n+N,1);
   u(1,k) = u_predicted(1,k) + K * (-xk + x_predicted(1:n,k)); % optimal input with stabilizing feedback control law K
   %u(:,k) = max(min(u(:,k),umax), umin); % check plausibility
   
   % create random noise based on scaled disturbance set W_scaled
   k_epsilon = n*(N+1)+N+n; %+n_w+n+2*n*(n+n_epsilon)+2*n;
   Phi_w = diag(z1((k_epsilon+1):(k_epsilon+n_w)));
   W_scaled = zonotope(W.c,W.G*Phi_w);
   w(:,k) = randPoint(W_scaled);

   % propagate to the system
   x(:,k+1)= A * x(:,k) + B * u(:,k) + [1;1] .* w(:,k);
   zk = z1;
end


% plotting response
figure(1)
time = (0:NT);
subplot(2,1,1)
plot(time,x(1,:),'r.-','LineWidth',.7) 
hold on
plot(time,x(2,:),'k.-','LineWidth',.7)
hold on
stairs(time,x_ref(1,1:(NT+1)), 'b--', 'LineWidth',.7)
plot(time,x_predicted(1,:), 'm.-', 'LineWidth', .7)
legend('$x_1$','$x_2$','$x_{reference}$','$x_{predicted/error}$','Interpreter','latex');
%axis ([0 50-10 10])
xlabel('$k$','Interpreter','latex');ylabel('$\textbf{x}_{k}$','Interpreter','latex');
grid on
ax = gca;
set(gca,'xtick',[0:5:NT])
%set(gca,'ytick',(-2*max(x_min(:))):10:(1.5*max(x(1,:))))
ax.GridAlpha = 1
ax.GridLineStyle = ':'
subplot(2,1,2)
stairs(time(1:end-1),u,'r.-','LineWidth',.7)
hold on
stairs(time(1:NT),u_predicted(1,:), 'b.-', 'LineWidth', .7)
%axis([0 50 -10 0])
legend('$u$','$u_{predicted}$','Interpreter','latex');
xlabel('$k$','Interpreter','latex');ylabel('${u}_{k}$','Interpreter','latex');
grid on
ax = gca;
ylim([2*umin 2*umax])
%set(gca,'xtick',[0:5:NT])
%set(gca,'ytick',(-2):.5:(2))
ax.GridAlpha = 1
ax.GridLineStyle = ':'
print -dsvg lmpc1
figure(2)
plot(zonotope(Xc.c,Xc.G*diag(zk((2*(N+1)+N+1):(2*(N+1)+N+n))))) % plotting scaled Xc set
hold on
plot(Xc)
plot(x(1,:),x(2,:))
plot(x(1,1),x(2,1), 'r')
plot(x(1,NT),x(2,NT), 'x')


function [Error_Zonotope, Xc_robust, Uc_robust] = compute_disturbance_invariance_set(A,B,K,W,N,Xc,Uc)
    % calculate the bounding set Epsilon as an outer approximation of the
    % MRPI set satisfying minkowski sum of A_cl*Epsilon and disturbance W as a subset of Epsilon 
    A_cl = A - B .* K;
    Epsilon = W;
    for j = 1:N
        Epsilon = Epsilon + A_cl^j * W;
    end
    Error_Zonotope = Epsilon;
    Xc_robust = minkDiff(Xc,Epsilon);
    Uc_robust = minkDiff(Uc,K * Epsilon);
    %{
    figure; hold on
    plot(e,[1 2], 'g');
    plot(Xc,[1 2], 'b');
    %hold on
    plot(Xc_robust,[1 2], 'r'); %plot(Uc_robust,[1 2]);
    %}
end





