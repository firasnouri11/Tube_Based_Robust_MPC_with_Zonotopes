% This code demonstrates a simple Tube-Based MPC to control a system
% under disturbances while keeping the actual system within a "tube"
% around the nominal trajectory.


addpath(genpath('C:\Users\firas\CORA'));


%Initialize the system states and define disturbance zonotope.
%Compute the nominal control input using a nominal MPC.
%Calculate the error between actual and nominal states.
%Adjust the control input using feedback.
%Propagate the disturbance zonotope through the system dynamics.
%Apply constraints robustly using the propagated zonotope.
%Update the nominal and actual states.
%Repeat the process for each time step.



%xn_nom : his is the ideal state trajectory as it assumes there are no disturbances or uncertainties affecting the system.
%the nominal system evolves like : xnom(t+1) = A*x_nom(t) + B*u_nom(t) 
% it is a reference to guide the actual system.

%xactual: is the actual state of the system with the effects of
%disturbances and errros  . it evolves : xact(t+1) = A*x_act(t) + B*u_act(t) + w(t) 



% System Dynamics: x(t+1) = A*x(t) + B*u(t) + w(t)
A = [1 1; 0 1];  % State transition matrix
B = [0.5; 1];    % Control input matrix

x0 = [2; 0];     % Initial state (both nominal and actual start here), hshould be nominal 



% Define Disturbance Zonotope 
disturbance_center = 0;  
disturbance_generators = 0.1;  % Disturbance range [-0.1, 0.1]
disturbance_zonotope = [disturbance_center, disturbance_generators]; % Simple zonotope





% Cost Function Parameters
Q = eye(2);  % Penalizes state deviations
R = 1;       % Penalizes control effort



% Constraints
u_max = 1;          % Maximum allowable control input
x_max = [10; 10];   % Maximum allowable state values


K = -dlqr(A, B, Q, R); % dlqr ?
%compute a state feedback gain matrix  K 
T = 15;  % Total simulation time steps


x_nom = x0;           % Nominal state trajectory
x_actual = x0;        % Actual state trajectory (subject to disturbance)
x_nom_traj = x0;      % Store nominal states
x_actual_traj = x0;   % Store actual states
u_traj = [];          % Store control inputs


     for t = 1:T

    u_nom = -0.5 * x_nom(1);  % Use a nominal system (ignores disturbances) to calculate an optimal control u_nom(t)
    u_nom = max(min(u_nom, u_max), -u_max);  % Ensure within [-u_max, u_max]

    e = x_actual - x_nom;  % Calculate error between actual and nominal states

    u_actual = u_nom + K * e;  % Apply error correction to nominal control 
    %K is the feedback gain to bring the actual state back toward the nominal trajectory.

    u_actual = max(min(u_actual, u_max), -u_max);
    
 
        x_nom_next = A * x_nom + B * u_nom;

w_zono = disturbance_zonotope * [-1, 1];  % Zonotope disturbance
        x_actual_next = A * x_actual + B * u_actual + w_zono;

        %insure both states and control input stay within the liimits

x_actual_next = max(min(x_actual_next, x_max), -x_max);   




 x_nom_traj = [x_nom_traj, x_nom_next];  % Store nominal trajectory
    x_actual_traj = [x_actual_traj, x_actual_next];  % Store actual trajectory
    u_traj = [u_traj, u_actual];  % Store control inputs

 x_nom = x_nom_next;
    x_actual = x_actual_next;
end

      
