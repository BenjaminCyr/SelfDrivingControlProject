load('Input.mat');
load('TestTrack.mat');

Y = forwardIntegrateControlInput(U);

figure;
plot(TestTrack.bl(1,:), TestTrack.bl(2,:), 'b')
hold on
plot(TestTrack.br(1,:), TestTrack.br(2,:), 'b')
hold on
plot(TestTrack.cline(1,:), TestTrack.cline(2,:), 'y')
hold on
plot(Y(:,1),Y(:,3), 'r');

info = getTrajectoryInfo(Y, U, [], TestTrack)