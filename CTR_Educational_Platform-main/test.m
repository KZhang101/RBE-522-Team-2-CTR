%% Set up
clear; clc;

startPose = Pose(0, 0, 0, 0, 0, 0);
drive_bot = Drive(startPose);
% dist = -10;

% drive_bot.travel_for(dist,dist,dist,0,0,0)
drive_bot.set_current_pose_as_home();

%% Test linear 1
% drive_bot.travel_for(0,0,1,0,0,0)

%% Test linear 2
% drive_bot.travel_for(0,-1,0,0,0,0)
% 
%drive_bot.travel_for(20,0,0,0,0,0)
% %% Test linear 3
% drive_bot.travel_for(0,0,10,0,0,0)  
% 
% %% Test rotation 1
% drive_bot.travel_for(0,0,0,-5,0,0)
% 
% %% Test rotation 2
% drive_bot.travel_for(0,0,0,0,-10,0)
% 
% %% Test rotation 3
% drive_bot.travel_for(0,0,0,0,0,90)