clc; close all


% pts = [0.65,-0.2;-0.5,0.7;0,0; 0.3, 0.7; -0.2,1.25; 0.85,1.8; 1.1,0.5; 0.5,1.2];
% %pts = 10*randn(6,2);
% x = 1:3;
% y = 1:3;
% [X,Y] = meshgrid(x,y);
% int = 1:size(X,1);
% X(int,int) = X(int,int) + 1*randn(size(X(int,int)));
% Y(int,int) = Y(int,int) +1*randn(size(X(int,int)));
% %X(2,1) = 0;
% %pts = [X(:),Y(:)];
% G =delaunayTriangulation(pts);
% Gt = triangleGrid(G.Points, G.ConnectivityList)
% Gp = pebi(Gt);
% [CC,r] = circumcenter(G);
% 
% circles = {};
% n = 100;
% theta = linspace(0,2*pi,n)';
% for i = 1:length(CC)
%    circles{i} = repmat(CC(i,:),n,1) + r(i) * [cos(theta),sin(theta)];
% end
% 
% plotGrid(Gt, 'facecolor', 'none')
% hold on
% plot(pts(:,1), pts(:,2), '.', 'markersize', 30)
% axis equal off
% figure()
% hold on
% plotGrid(Gp,'facecolor','none')
% plot(pts(:,1), pts(:,2), '.', 'markersize', 30)
% axis equal off
% figure()
% plotGrid(Gt, 'facecolor', 'none')
% hold on
% plotGrid(Gp,'facecolor','none')
% plot(pts(:,1), pts(:,2), '.', 'markersize', 30)
% for i = 1:numel(circles)
%    y = circles{i}
%    %plot(y(:,1), y(:,2), 'color', [0.7,0.7,0.7])
% end
% 
% %plot(CC(:,1), CC(:,2),'o')
% %plot(Gp.nodes.coords(:,1), Gp.nodes.coords(:,2),'ro')
% axis equal off


% %% Locally delaunay
% pts = [1,0; 1.75, 0.5; 1, 1; 0.25, 0.5];%;0.05,0.8;0.95,1.2];
% 
% %%Locally delaunay
% lines = [1,2; 3, 1; 4, 1; 2,3; 3,4];
% triangle2 = [2,3,4,2];
% cc = [1,0.5];
% r = 0.5
% 
% %%not locally delaunay
% % cc = [1, 0.8125];
% % r = 0.8125;
% % lines = [1,2; 2, 4; 4, 1; 2,3; 3,4];
% % triangle2 = [2,3,4,2];
% 
% n=100
% theta = linspace(0,2*pi,n)';
% circle = repmat(cc,n,1) +r*[cos(theta), sin(theta)]; 
% 
% %col =  [0, 0.4470, 0.7410]
% figure
% hold on
% for i = 1:size(lines,1)
%     plot(pts(lines(i,:),1), pts(lines(i,:),2), 'color','k')
% end
% plot(pts(:,1), pts(:,2), '.', 'markersize', 30)
% plot(circle(:,1),circle(:,2), 'color', [0.7,0.7,0.7])
% 
% 
% axis equal off



% %% Sufficient condition well gridding.
% pts = [0,0.5;1,0.5;1+1/sqrt(2),0.5-1/sqrt(2)];
% 
% lines = [1,2;2,3];
% cc = [0.5,0.5;1,0.5];
% r = [0.5,1];
% n = 100;
% theta = (linspace(0,2*pi, n))';
% 
% circle1 = repmat(cc(1,:),n,1) + r(1)*[cos(theta), sin(theta)]; 
% circle2 = repmat(cc(2,:),n,1) + r(2)*[cos(theta), sin(theta)]; 
% 
% col =  [0, 0.4470, 0.7410]
% figure
% hold on
% for i = 1:size(lines,1)
%     plot(pts(lines(i,:),1), pts(lines(i,:),2), 'color','k')
% end
% plot(pts(:,1), pts(:,2), '.', 'markersize', 30, 'color', col)
% plot(circle1(:,1),circle1(:,2), 'color', [0.7,0.7,0.7])
% plot(circle2(:,1),circle2(:,2), 'color', [0.7,0.7,0.7])
% 
% axis equal off




