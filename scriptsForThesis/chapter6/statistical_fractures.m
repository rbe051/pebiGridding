close all; clear
%% load dataset

load('dataset/statistical_fractures.mat')

l = mat2cell(fl, ones(size(fl,1),1), size(fl,2));
l = l';
offset = 2;
l = cellfun(@(c) reshape(c', 2,[])' + offset,l,'un',false);

plotLinePath(l)
%% Create grid
eps = [3,5,1,1,1,2,2,1.5,1,1.5,1.5,1.5,1,1.5,1.5,2.5,1.5];
%refPts = [15.4,77.5; 13,69.6;23.9,96.7;23.8,95.9;23.7,95.2;23.7,94.4;...
%          23.5,92.7;23.3,90.8;7.1,85;11,85.6;11.1,86.5;7,64.2;26.95,75;...
%          24, 68;18.8,67.2;22.9,22.8;24.8,36.9] + offset;
refPts = [17.4,79.5; 15.0,71.6; 25.9,98.7; 25.8,97.9; ...
          25.7,97.2; 25.7,96.4; 25.5,94.7; 25.3,92.8; ...
          9.1,87.0;  13.0,87.6; 13.1,88.5; 9.0 ,66.2; ...
          28.95,77;  26.0,70.0; 20.8,69.2; 24.9,24.8; 26.8,38.9];
amp = [0.5,0.7,0.3,0.3,0.4,0.4,0.4,.6,.3,.45,.45,.2,.3,0.25,0.55,0.45,0.6];
faultRho =@(p) min(ones(size(p,1),1), ...
                   min(bsxfun(@times, amp,exp(bsxfun(@rdivide, pdist2(p,refPts),eps))),[],2));
%%
%close all
G = pebiGrid(15,[35,120],'faultLines',l,'faultGridFactor',1/40,...
             'circleFactor',0.62,'faultRefinement',true,'faultEps',5,...
             'faultRho', faultRho)
           

%% plot


figure(); hold on
plotGrid(G)
plotLinePath(l,'--')
axis equal

%% Save
save('statistical_fractures_edgeCentered')