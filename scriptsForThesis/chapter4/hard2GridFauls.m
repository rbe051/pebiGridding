clear; close all
%% Axis
ax = [0.2,0.8,0.3,0.8];
%% Sharp intersections
f = {[0.2,0.45;0.8,0.55],...
     [0.2,0.55;0.8,0.45]};
color = get(gca, 'colorOrder');
plotLinePath(f, 'color',color(2,:))
axis equal off
axis(ax)

print('../../../../master/thesis/fig/ch04/hardFaultGridSharp','-dpng')
%% multiple intersections
f = {[0.2,0.5;0.8,0.6],...
     [0.2,0.7;0.8,0.4],...
     [0.5,0.2;0.5,0.8]};
color = get(gca, 'colorOrder');
figure()
plotLinePath(f, 'color',color(2,:))
axis equal off
axis(ax)
print('../../../../master/thesis/fig/ch04/hardFaultGridMultInt','-dpng')

%% barley intersecting 
f = {[0.2,0.7;0.8,0.7],...
     [0.5,0.3;0.5,0.72]};
color = get(gca, 'colorOrder');
figure()
plotLinePath(f, 'color',color(2,:))
axis equal off
axis(ax)
print('../../../../master/thesis/fig/ch04/hardFaultGridBarInt','-dpng')