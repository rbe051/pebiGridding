clear; close all

%% create a well

w = {[0.2,0.5;0.8,0.5], ...
     [0.2,0.2;0.8,0.8]};
f = {[0.7,0.2;0.7,0.85]};
gs = 0.01;

d = {@(p) 0.005*ones(size(p,1),1)};
wCut = [0,0];
[wp,~, pp] = createWellGridPoints(w, gs,'protLayer',true,'protD',d,'wCut',wCut);

%% plot
figure();hold on

plot(wp(:,1), wp(:,2),'.','markersize',15);
plot(pp(:,1), pp(:,2),'.','markersize',15);


%% Create Grid
G = compositePebiGrid(gs,[1,1], 'wellLines',w,'faultLines',f,...
                      'protLayer',true,'protD',d)

%% Plot Grid

plotGrid(G)


%% Case 2 Wavy well

x = 1:0.01:9;
y = 5+4*sin(2*pi*x/10);
w = {[x',y']};

gs = 0.1;

d = {@(p) 0.05*rand(size(p,1),1)};


%% Create Grid
G = compositePebiGrid(gs,[10,10], 'wellLines',w,...
                      'protLayer',true,'protD',d)

%% Plot Grid

plotGrid(G)

