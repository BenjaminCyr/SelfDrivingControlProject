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
tspan = 0:dt:T(end);

%Interpolate;
Z_ref=interp1(T, traj_total(:,2:9),tspan, 'spline')';
Y_ref = Z_ref(1:6,:);
U_ref = Z_ref(7:8,:);

x_c = @(t) Z_ref(1, round(t/dt+1));
u_c = @(t) Z_ref(2, round(t/dt+1));
y_c = @(t) Z_ref(3, round(t/dt+1));
v_c = @(t) Z_ref(4, round(t/dt+1));
psi_c = @(t) Z_ref(5, round(t/dt+1));
r_c = @(t) Z_ref(6, round(t/dt+1));

state_c = @(t) [x_c(t); u_c(t); y_c(t); v_c(t); psi_c(t); r_c(t);];

Fx_c = @(t) Z_ref(7, round(t/dt+1));
df_c = @(t) Z_ref(8, round(t/dt+1));

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
Q = diag([1 1 1 1 1 1]);
R = diag([.0000000001 1]);

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
plot(Z_ref(1,:),Z_ref(3,:), 'g');
hold on


Y = zeros(6,length(tspan));
U = zeros(2,length(tspan));
Y(:,1) = Y_ref(:,1);


for i = 1:length(tspan)-1 
   disp(i)
   error = Y(:,i)-Y_ref(:,i)
   window = min(npred, length(tspan)-i);
   if window > 1
       [Aeq, beq] = eq_cons(i,A,B,window, error);
       [Lb, Ub] = bound_cons(i, U_ref, [lb_Fx ub_Fx; lb_df ub_df], window);

       H_section = H((npred-window)*6+1 : 6*npred+6+2*window, (npred-window)*6+1 : 6*npred+6+2*window);
       optimal_z = quadprog( H_section, c(1:8*window+6), [], [], Aeq, beq, Lb, Ub, [], options);
   end
   
   y_mpc = optimal_z(1:6*(window+1));
   u_mpc = optimal_z(6*(window+1)+1:6*(window+1)+2);
   u = U_ref(:,i) + u_mpc;
   %use ode45 to simulate nonlinear system, f, forward 1 timestep
   [~,ytemp]=ode45(@(t,x) bike_odefun(x,u, false),[0 dt],Y(:,i));
   Y(:,i+1) = ytemp(end,:)';
   U(:,i) = u;
end

plot(Y(1,1:i),Y(3,1:i), 'r');

U = [U(2,:)' U(1,:)'];
save('Input.mat', 'U');


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
    beq( 1:n ) = x0;
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
        Lb(2*i+6*(npred+1)-1:2*i+6*(npred+1)) = input_range(:,1) - U_ref(:,idx+i-1);
        Ub(2*i+6*(npred+1)-1:2*i+6*(npred+1)) = input_range(:,2) - U_ref(:,idx+i-1);
    end
end

% function x_close = closest_point(x, traj, index)
%     distances = vecnorm(x - traj(index));
%     distance2 = vecnorm(x - traj(index+1));
%     while distance2 > distance
%         distance = distance2;
%         index = index + 1;
%         distance2 = vecnorm(x - traj(index+1));
%     end
%     x_close = traj(:,index);
% end

% function x_close = closest_point(x, traj)
%     distances = sqrt((x(1)-traj(1,:)).^2 + (x(3)-traj(3,:)).^2);
% 
%     [~, ind] = min(distances);
%     x_close = traj(:,ind);
% end