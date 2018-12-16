clear all; close all;

load('TestTrack.mat')
cline = TestTrack.cline;
leftline = TestTrack.bl;
rightline = TestTrack.br;
plot(TestTrack.bl(1,:), TestTrack.bl(2,:), 'k')
hold on
plot(TestTrack.br(1,:), TestTrack.br(2,:), 'k')
hold on
plot(TestTrack.cline(1,:), TestTrack.cline(2,:), 'y')
hold on



nsteps = 10+1;

min_T = 0.01*nsteps;
max_T = 10;

highest_x = 1500;
lowest_x = 0;

highest_y = 850;
lowest_y = -200;

ub_delta = 0.5;
lb_delta = -0.5;

ub_Fx = 2500;
lb_Fx = -5000;

weights = [1 .00001 1 5];


ub = [repmat([highest_x Inf highest_y Inf Inf Inf], 1, nsteps) repmat([ub_Fx ub_delta], 1, nsteps-1) max_T];

lb = [repmat([lowest_x -Inf lowest_y -Inf -Inf -Inf], 1, nsteps) repmat([lb_Fx lb_delta], 1, nsteps-1) min_T];
 

% %%%%%%%%%%%%%%% no need to change these lines  %%%%%%%%%%%%%%%%%%%%%%
options = optimoptions('fmincon', 'ConstraintTolerance', 1e-6, 'SpecifyConstraintGradient',true,...
                       'SpecifyObjectiveGradient',true,...
                       'MaxFunctionEvaluations', 3000, 'MaxIterations', 10000, 'Display','none',...
                       'CheckGradients',false);
%                          'Display','iter')

is = [287 5 -176 0 2 0];
fs = [0 0 0 0 0 0];

% x = mean([cline(1,:); leftline(1,:)]);
% y = mean([cline(2,:); leftline(2,:)]);
x = mean([cline(1,:); rightline(1,:)]);
y = mean([cline(2,:); rightline(2,:)]);
cline = [x; y];


final_point = [cline(1, end)+10; cline(2, end) + 20];
cline = [cline final_point];

num_iter = 100;
traj_total  = [];
dt_total = [];
last_time = 0;
index = 1;
i = 1;

last_input = [0 0];

while is(3) < cline(2,end)
    [close_point, index] = closest_point([is(1) is(3)], cline, index);
    fs = [close_point(1) 0 close_point(2) 0 0 0];
    
    plot(is(1),is(3), 'ro')
    plot(fs(1),fs(3), 'go')
    
    x0 = linspace(is(1), fs(1), nsteps);
    y0 = linspace(is(3), fs(3), nsteps);

    z0 = zeros(nsteps*8-2+1, 1);
    z0(1:nsteps*6) = repmat(is, 1, nsteps);
    z0(1:6:nsteps*6) = x0;
    z0(3:6:nsteps*6) = y0;
    
    z0(end) = 2;

    nc=@(z) nonlcon(z, nsteps, is);
    cf=@(z) costfun(z, nsteps, fs, weights);

    fprintf('%d:\t', i);
    [z,fval,~,output] = fmincon(cf,z0,[],[],[],[],lb',ub',nc,options);
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
    T = z(end);
    dt = T/(nsteps - 1);

    dt_total = linspace(last_time, last_time+T, nsteps);
    traj_total = [traj_total;dt_total(1:end-1)' Y0(1:end-1, :) U']; 
    last_time = last_time+T;
    
    plot(Y0(:,1),Y0(:,3), 'g')
    hold on
    
    is = Y0(end,:);
    last_input = [Fx(end) delta(end)];
    i = i + 1;
end

traj_total(end,1)

% save('Trajectory.mat', 'traj_total');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [pt, current_index] = closest_point(pt_in, boundary, current_index)
    distance = sqrt((pt_in(1) - boundary(1,current_index))^2 + (pt_in(2) - boundary(2,current_index))^2);
    while distance < 20 && current_index < length(boundary)       
       current_index = current_index+1;
       distance = sqrt((pt_in(1) - boundary(1,current_index))^2 + (pt_in(2) - boundary(2,current_index))^2);
    end
    pt = boundary(:,current_index);
end