clear all; close all;

load('ROB599_ControlsProjectTeam35.mat');
load('TestTrack.mat');

Y = forwardIntegrateControlInput(ROB599_ControlsProject_part1_input);

figure;
plot(TestTrack.bl(1,:), TestTrack.bl(2,:), 'b')
hold on
plot(TestTrack.br(1,:), TestTrack.br(2,:), 'b')
hold on
plot(TestTrack.cline(1,:), TestTrack.cline(2,:), 'y')
hold on
plot(Y(:,1),Y(:,3), 'r');

info = getTrajectoryInfo(Y, ROB599_ControlsProject_part1_input, [], TestTrack)

% Y = ROB599_ControlsProject_part2_Team35(TestTrack, [])

% info = getTrajectoryInfo(Y, ROB599_ControlsProject_part1_input, [], TestTrack)