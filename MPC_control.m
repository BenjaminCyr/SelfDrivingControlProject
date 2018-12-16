close all; clear all;

%% parameters
ub_df = 0.5;
lb_df = -0.5;

ub_Fx = 2500;
lb_Fx = -5000;

dt=0.01;

options = optimoptions('quadprog', 'Display', 'none');

%% load reference trajectory
load('Trajectory.mat')

T = traj_total(:,1);
%time span
tspan = [0:dt:T(end)]';

% First Turn
% start_point = 1000;
% end_point = 3000;

%Final Turns
% start_point = 10000;
% end_point = 13500;

%Middle Turn
% start_point = 8500;
% end_point = start_point+1000;


start_point = 1;
end_point = length(tspan);

%Interpolate;
Z_ref=interp1(T, traj_total(:,2:9),tspan);
Y_ref = Z_ref(:,1:6);
U_ref = Z_ref(:, 7:8);

x_c = @(t) Z_ref(round(t/dt+1), 1);
u_c = @(t) Z_ref(round(t/dt+1), 2);
y_c = @(t) Z_ref(round(t/dt+1), 3);
v_c = @(t) Z_ref(round(t/dt+1), 4);
psi_c = @(t) Z_ref(round(t/dt+1), 5);
r_c = @(t) Z_ref(round(t/dt+1), 6);

state_c = @(t) [x_c(t); u_c(t); y_c(t); v_c(t); psi_c(t); r_c(t);];

Fx_c = @(t) Z_ref(round(t/dt+1), 7);
df_c = @(t) Z_ref(round(t/dt+1), 8);

input_c = @(t) [Fx_c(t) df_c(t)];

A_c = @(t) calculate_A(state_c(t), input_c(t));
B_c = @(t) calculate_B(state_c(t), input_c(t));

A=@(i) eye(6) + dt*A_c((i-1)*dt);

B=@(i) dt*B_c((i-1)*dt);

%% 2.2 Number of decision variables for colocation method
npred=10;

%% 2.5 simulate controller working from initial condition [0.25;-0.25;-0.1]
%use ode45 to between inputs
% we begin by defining the cost function
Q = diag([100 1 100 1 1 1]);
R = diag([.0000001 1]);

H = zeros(6*(npred+1) + 2*npred);
c = zeros(6*(npred+1) + 2*npred,1);

H((6*npred+1):(6*npred+6),(6*(npred)+1):(6*npred+6)) = Q;
for i = 1:npred
    H((6*(i-1)+1):(6*(i-1)+6),(6*(i-1)+1):(6*(i-1)+6)) = Q;
    H((2*(i-1)+6*(npred+1)+1):(2*i+6*(npred+1)), (2*(i-1)+6*(npred+1)+1):(2*i+6*(npred+1))) = R;
end




load('TestTrack.mat')

plot(TestTrack.bl(1,:), TestTrack.bl(2,:), 'b')
hold on
plot(TestTrack.br(1,:), TestTrack.br(2,:), 'b')
hold on
plot(TestTrack.cline(1,:), TestTrack.cline(2,:), 'y')
hold on
plot(Z_ref(:,1),Z_ref(:,3), 'g');
hold on


Y = zeros(length(tspan), 6);
U = zeros(length(tspan), 2);
Y(start_point,:) = Y_ref(start_point,:);

last_U = [U_ref(start_point,2) U_ref(start_point,1)];

tic;
for i = start_point:end_point %length(tspan)-1 
   error = Y(i, :)-Y_ref(i,:);
   fprintf('Error %d:  [%f %f %f %f %f %f]\n', i, error);
   window = min(npred, length(tspan)-i);
   if window > 1
       [Aeq, beq] = eq_cons(i,A,B,window, error);
       [Lb, Ub] = bound_cons(i, U_ref, [lb_Fx ub_Fx; lb_df ub_df], window);
       
       H_section = H((npred-window)*6+1 : 6*npred+6+2*window, (npred-window)*6+1 : 6*npred+6+2*window);
       optimal_z = quadprog( H_section, c(1:8*window+6), [], [], Aeq, beq, Lb, Ub, [], options);
   end
   
   y_mpc = optimal_z(1:6*(window+1));
   u_mpc = optimal_z(6*(window+1)+1:6*(window+1)+2);
   current_U = [U_ref(i,2)+u_mpc(2) U_ref(i, 1)+u_mpc(1)];
   u = [last_U; current_U];
   %use ode45 to simulate nonlinear system, f, forward 1 timestep
   ytemp = forwardIntegrateControlInput(u, Y(i,:));
%    [~, ytemp] = ode45(@(t, x) bike_odefun(x, current_U), [0 dt], Y(i,:));
   Y(i+1,:) = ytemp(end,:);
%    
%    plot(Y(i:i+1, 1),Y(i:i+1, 3), 'r')
%    hold on
   
   U(i,:) = current_U;
   last_U = current_U;
end
toc;

plot(Y(start_point:end_point,1),Y(start_point:end_point, 3), 'r');

ROB599_ControlsProject_part1_input = U(start_point:end_point,:);
START = Y_ref(start_point,:);
save('ROB599_ControlsProjectTeam35.mat', 'ROB599_ControlsProject_part1_input', 'START', 'traj_total');


function [Aeq,beq]=eq_cons(idx,A,B,npred,x0)
    %build matrix for A_i*x_i+B_i*u_i-x_{i+1}=0
    %in the form Aeq*z=beq
    %initial_idx specifies the time index of initial condition from the reference trajectory 
    %A and B are function handles above
    
    n = size(A(1), 2); %size of state vector
    m = size(B(1), 2); % size of input vector
    
    %Aeq*[x_0, y_0, psi_0, ... ,u_0, delta_0, ...] =beq
    Aeq = zeros(n*(npred+1), n*(npred+1)+m*npred);
    beq = zeros(n*(npred+1), 1);
    
    Aeq( 1:n, 1:n ) = eye( n ); 
    beq( 1:n ) = x0';
    for i = 1:npred-1
        Aeq((i*n+1):((i+1)*n), (i*n+1):((i+1)*n)) = eye(n);
        Aeq((i*n+1):((i+1)*n), ((i-1)*n+1):(i*n)) = -A(idx+i-1);
        Aeq((i*n+1):((i+1)*n), ((i-1)*m+(npred+1)*n+1):(i*m+(npred+1)*n)) = -B(idx+i-1);
    end
end

function [Lb,Ub]=bound_cons(idx, U_ref, input_range, npred)
    %initial_idx is the index along uref the initial condition is at
    Lb = -Inf(6*(npred+1)+2*npred, 1);
    Ub = Inf(6*(npred+1)+2*npred, 1);
    
    for i = 1:npred
        Lb(2*i+6*(npred+1)-1:2*i+6*(npred+1)) = input_range(:,1) - U_ref(idx+i-1,:)';
        Ub(2*i+6*(npred+1)-1:2*i+6*(npred+1)) = input_range(:,2) - U_ref(idx+i-1,:)';
    end
end