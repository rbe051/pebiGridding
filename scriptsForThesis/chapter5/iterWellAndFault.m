clc; clear; close all


%% Create a well and fault Path

well  = {[0.2,0.2; 0.8,0.8]; ...
         [0.4,0.5; 0.4,0.8]};
fault = {[0.2,0.6; 0.8,0.4]};

%% Plot
figure(); hold on
plot(well(:,1), well(:,2))
plot(fault(:,1),fault(:,2))
axis equal
axis([0,1,0,1])

%% Split well and fault at intersection
