clear; close all
%addpath ../../voronoi2D/
%% Set fault
% 90 degrees
alpha = 90/2*pi/180;
a = 0.8 - 0.5;
b = a*tan(alpha);
f1 = {[0.5-a,0.5-b; 0.5+a,0.5+b],[0.5-a,0.5+b;0.5+a,0.5-b]};


% 45 degrees
alpha = 45/2*pi/180;
a = 0.8 - 0.5;
b = a*tan(alpha);
f2 = {[0.5-a,0.5-b; 0.5+a,0.5+b],[0.5-a,0.5+b;0.5+a,0.5-b]};

% 20 degrees
alpha = 20/2*pi/180;
a = 0.8 - 0.5;
b = a*tan(alpha);
f3 = {[0.5-a,0.5-b; 0.5+a,0.5+b],[0.5-a,0.5+b;0.5+a,0.5-b]};



%% Set grid parameters
rGs = 0.1;
fGf = 1.0;
pdims = [1,1];
%% Create grid
G1 = compositePebiGrid(rGs,pdims,'faultlines',f1,'faultgridfactor',fGf);
G2 = compositePebiGrid(rGs,pdims,'faultlines',f2,'faultgridfactor',fGf);
G3 = compositePebiGrid(rGs,pdims,'faultlines',f3,'faultgridfactor',fGf);


%% Plot
color = get(gca,'ColorOrder');
close all
figure(1);
plotGrid(G1,'facecolor','none');
plotLinePath(f1,'--','color',color(2,:),'linewidth',1);
axis equal off tight

figure(2)
plotGrid(G2,'facecolor','none');
plotLinePath(f2,'--','color',color(2,:),'linewidth',1);
axis equal off tight

figure(3);
plotGrid(G3,'facecolor','none')
plotLinePath(f3,'--','color',color(2,:),'linewidth',1);
axis equal off tight

%% Save
figure(1)
print('../../../../master/thesis/fig/ch06/showInterSectAt90','-depsc')
figure(2)
print('../../../../master/thesis/fig/ch06/showInterSectAtSharp','-depsc')
figure(3)
print('../../../../master/thesis/fig/ch06/showInterSectAtVerySharp','-depsc')
