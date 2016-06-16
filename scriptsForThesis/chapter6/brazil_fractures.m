close all; clear
%% load dataset

load('dataset/brazil_fractures.mat')


l = mat2cell(fl, ones(size(fl,1),1), size(fl,2));
l = l';
offset = 2;
l = cellfun(@(c) reshape(c', 2,[])' + offset,l,'un',false);

plotLinePath(l)
%% Create grid
protD = @(p) 0.1*ones(size(p,1),1);
G = pebiGrid(500,[1050,1050],'wellLines',l,'wellGridFactor',1/80,...
             'wellRefinement',true,'epsilon',50,'protLayer',true,...
             'protD', protD)
           

%% plot


figure(); hold on
plotGrid(G)
plotLinePath(l,'--')
axis equal

%% Save
