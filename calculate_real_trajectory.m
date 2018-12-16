clear all; close all;


Nobs = 15;

load('ROB599_ControlsProjectTeam35.mat');
load('TestTrack.mat');

Y_part1 = forwardIntegrateControlInput(ROB599_ControlsProject_part1_input);
info_part1 = getTrajectoryInfo(Y_part1, ROB599_ControlsProject_part1_input, [], TestTrack)

plot(Y_part1(:,1),Y_part1(:,3), 'b');

Xobs = generateRandomObstacles(Nobs, TestTrack);

figure;
plot(TestTrack.bl(1,:), TestTrack.bl(2,:), 'k')
hold on
plot(TestTrack.br(1,:), TestTrack.br(2,:), 'k')
hold on
plot(TestTrack.cline(1,:), TestTrack.cline(2,:), 'y')
hold on
for i = 1:length(Xobs)
    obstacle = Xobs{i};
    plot(obstacle(:,1), obstacle(:,2), 'r');
end

tic;
U = ROB599_ControlsProject_part2_Team35_fmincon(TestTrack, Xobs);
toc;

Y_part2 = forwardIntegrateControlInput(U);
info_part2 = getTrajectoryInfo(Y_part2, U, Xobs, TestTrack)

plot(Y_part2(:,1),Y_part2(:,3), 'g');
