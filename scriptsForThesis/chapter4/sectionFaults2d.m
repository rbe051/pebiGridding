clc; clear all; close all


%% Create a fault
fault    = [0.2,0.8; 0.5,0.5; 0.8,0.5];
dFault   = 0.1;
circFact = 0.6;

%% Interpolate fault
linesDist = sqrt(sum(diff(fault,[],1).^2,2));
linesDist = [0; linesDist]; % add the starting point
cumDist   = cumsum(linesDist);
dt        = cumDist(end)/ceil(cumDist(end)/dFault);
newPtsT   = 0:dt:cumDist(end);
CC        = interp1(cumDist, fault, newPtsT);

%% Calculate circle radiuses


d  = sqrt(sum((CC(2:end,:)-CC(1:end-1,:)).^2, 2));
CR = circFact*[d(1); (d(1:end-1) + d(2:end))/2; d(end)];
            
            
%% Calculate midpoints and circle offsets
a = (d.^2 - CR(2:end).^2 + CR(1:end-1).^2)./(2*d);
h = sqrt(CR(1:end-1).^2 - a.^2);
    
% Calculate normal vectors
n1 = (CC(2:end,:)-CC(1:end-1,:))./repmat(d,1,2); %Unit vector
n2 = [-n1(:, 2), n1(:,1)];                       %Unit normal


%% Calculate intersections
lPts = CC(1:end-1,:) + bsxfun(@times, a, n1)  ...
                     + bsxfun(@times, h,  n2);
rPts = CC(1:end-1,:) + bsxfun(@times, a, n1)  ...
                     - bsxfun(@times, h,  n2);
pts = [lPts; rPts];

%% Plotting
figure(); hold on
theta= linspace(0,2*pi)';
for i = 1:size(CC,1)
  X = repmat(CC(i,:),100,1) + repmat(CR(i),100,2).*[cos(theta), sin(theta)];
  plot(X(:,1),X(:,2),'k')
end
plot(fault(:,1), fault(:,2));
plot(pts(:,1),pts(:,2),'.','markersize',15);
axis equal tight

%print('../../../../master/thesis/fig/ch04/fault2DexampleImp','-depsc')