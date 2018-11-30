clear all; close all;

load('TestTrack.mat')

cline = TestTrack.cline;
% T = 0.5
% dt = 0.01;
nsteps = 10+1;

min_T = 0.01*nsteps;
max_T = 10;

highest_x = 1500;
lowest_x = 0;
% highest_x = 1;
% lowest_x = -100;

highest_y = 850;
lowest_y = -200;
% highest_y = 1000;
% lowest_y = -200;

ub_delta = 0.5;
lb_delta = -0.5;

ub_Fx = 2500;
lb_Fx = -5000;

% weights = [1 0.1 1 0 0 0 0 0 nsteps];
weights = [1 .001 .1];

% 
% %remember the format for z is as follows:
% %z=[x0 y0 th0 x1 y1 th1 ... xn yn thn u0 d0 ... u(n-1) d(n-1) T]';
%     
% %1.1

ub = [repmat([highest_x Inf highest_y Inf Inf Inf], 1, nsteps) repmat([ub_Fx ub_delta], 1, nsteps-1)];

lb = [repmat([lowest_y 0 lowest_y -Inf -Inf -Inf], 1, nsteps) repmat([lb_Fx lb_delta], 1, nsteps-1)];
 
% %1.4
% %%%%%%%%%%%%%%% no need to change these lines  %%%%%%%%%%%%%%%%%%%%%%
options = optimoptions('fmincon', 'ConstraintTolerance', 1e-6, 'SpecifyConstraintGradient',true,...
                       'SpecifyObjectiveGradient',true,...
                       'MaxFunctionEvaluations', 1000, 'MaxIterations', 3000, 'Display','off',...
                       'CheckGradients',false);
%                          'Display','iter')
% 

% u0 = [ones(1, nsteps-1); zeros(1, nsteps-1)];
% u0_t=@(t) [interp1(0:dt:(nsteps-2)*dt,u0(1,:),t,'previous','extrap');...
%         interp1(0:dt:(nsteps-2)*dt,u0(2,:),t,'previous','extrap')];
    
% x0 = ode1(@(t,x) bike_odefun(x, u0_t(t)), 0:dt:T, initial_state);
T = 2;

is = [287 5 -176 0 2 0];
% initial_state = [0 1 0 0 0 0];
fs = [0 0 0 0 0 0];

plot(TestTrack.bl(1,:), TestTrack.bl(2,:), 'b')
hold on
plot(TestTrack.br(1,:), TestTrack.br(2,:), 'b')
hold on
plot(TestTrack.cline(1,:), TestTrack.cline(2,:), 'y')
hold on

final_point = [cline(1, end); cline(2, end) + 20];
cline = [cline final_point];

num_iter = 100;
% U_total = zeros(2, num_iter*length(0:0.01:T));
% U_total = [];
traj_total  = [];
index = 1;
i = 1;

while is(3) < cline(2,end) 
    [close_point, index] = closest_point([is(1) is(3)], cline, index);
    fs = [close_point(1) 0 close_point(2) 0 0 0];
    
%     disp([is(1) is(3)]);
%     disp([fs(1), fs(3)]);
    
    x0 = linspace(is(1), fs(1), nsteps);
    y0 = linspace(is(3), fs(3), nsteps);
    
%     vec = [is(1)-fs(1) fs(3)-is(3)]
    df_guess = 0;
    U0 = repmat([100 df_guess], 1, nsteps-1);

    z0 = zeros(nsteps*8-2, 1);
    z0(1:nsteps*6) = repmat(is, 1, nsteps);
    z0(1:6:nsteps*6) = x0;
    z0(2:6:nsteps*6) = ones(1,nsteps);
    z0(3:6:nsteps*6) = y0;
%     z0(nsteps*6+1:nsteps*8-2) = U0;

    nc=@(z) nonlcon(z, nsteps, is, T);
    cf=@(z) costfun(z, nsteps, fs, weights);

    A = [];
    b = [];

    fprintf('%d:\t', i);
    [z,fval,~,output] = fmincon(cf,z0,A,b,[],[],lb',ub',nc,options);
    fprintf('fval = %e\tFeasibility = %e\n', fval, output.constrviolation);
    

    Y0=reshape(z(1:6*nsteps),6,nsteps)';
    U=reshape(z(6*nsteps+1:8*nsteps-2),2,nsteps-1);

    x_opt = z(1:6:(nsteps)*6);
    u_opt = z(2:6:(nsteps)*6);
    y_opt = z(3:6:(nsteps)*6);
    v_opt = z(4:6:(nsteps)*6);
    psi_opt = z(5:6:(nsteps)*6);
    r_opt = z(6:6:(nsteps)*6);

    Fx = z((nsteps)*6+1:2:(nsteps)*8-2);
    delta = z((nsteps)*6+2:2:(nsteps)*8-2);
    dt = T/(nsteps - 1);

    % Y0=reshape(z(1:6*nsteps),6,nsteps)';
    % U=reshape(z(6*nsteps+1:end-1),2,nsteps-1);



    u=@(t) [interp1(0:dt:(nsteps-2)*dt,U(2,:),t,'previous','extrap');...
            interp1(0:dt:(nsteps-2)*dt,U(1,:),t,'previous','extrap')];
    
    U_i = u(0:0.01:T);
%     U_total(:, (i - 1)*length(U_i)+1:i*length(U_i)) = U_i;
    traj_total = [traj_total; Y0(1:end-1, :) U']; 
%     U_total = [U_total U_i];
%     [Y1, T1] = forwardIntegrateControlInput(U_i', is);
    
    plot(Y0(:,1),Y0(:,3), 'r')
    hold on
    
    is = Y0(end,:);
    i = i + 1;
end


save('Trajectory.mat', 'traj_total');
% is = [287 5 -176 0 2 0];    
% [Y1, T1] = forwardIntegrateControlInput(U_total', is);


xlabel('x');
ylabel('y');
legend('left', 'right', 'center', 'fmincon trajectory','ode45 trajectory using x0 = [0;0;0]');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [pt, current_index] = closest_point(pt_in, boundary, current_index)
    distance = sqrt((pt_in(1) - boundary(1,current_index))^2 + (pt_in(2) - boundary(2,current_index))^2);
    while distance < 20 && current_index < length(boundary)       
       current_index = current_index+1;
       distance = sqrt((pt_in(1) - boundary(1,current_index))^2 + (pt_in(2) - boundary(2,current_index))^2);
    end
    pt = boundary(:,current_index);
end