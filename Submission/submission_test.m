clear all; close all;

num_runs = 3;

load('ROB599_ControlsProject_part1_Team35.mat');
load('TestTrack.mat');

figure;
plot(TestTrack.bl(1,:), TestTrack.bl(2,:), 'k')
hold on
plot(TestTrack.br(1,:), TestTrack.br(2,:), 'k')
hold on
plot(TestTrack.cline(1,:), TestTrack.cline(2,:), 'y')
hold on
Y_part1 = forwardIntegrateControlInput2(ROB599_ControlsProject_part1_input);
info_part1 = getTrajectoryInfo(Y_part1, ROB599_ControlsProject_part1_input, [], TestTrack)

plot(Y_part1(:,1),Y_part1(:,3), 'b');

infos{1,num_runs} = [];
times{1, num_runs} = [];
Nobs = zeros(1, num_runs);
percents = zeros(1, num_runs);
for i = 1:num_runs
    figure;
    plot(TestTrack.bl(1,:), TestTrack.bl(2,:), 'k')
    hold on
    plot(TestTrack.br(1,:), TestTrack.br(2,:), 'k')
    hold on
    plot(TestTrack.cline(1,:), TestTrack.cline(2,:), 'y')
    hold on
    
    Nobs(i) = randi([10 25]);
    Xobs = generateRandomObstacles(Nobs(i), TestTrack);
    for j = 1:length(Xobs)
        obstacle = Xobs{j};
        plot(obstacle(:,1), obstacle(:,2), 'r');
    end

    fprintf("%d: %d Obstacles\n", i, Nobs(i));
    comp_time = tic;
    U = ROB599_ControlsProject_part2_Team35(TestTrack, Xobs);
    time = toc(comp_time)
    times{i} = time;
    
    Y_part2 = forwardIntegrateControlInput2(U);
    plot(Y_part2(:,1),Y_part2(:,3), 'g');
    
    info = getTrajectoryInfo(Y_part2, U, Xobs, TestTrack)
    infos{i} = info;
    percents(i) = info.percent_of_track_completed;
end

avg_percent = mean(percents)