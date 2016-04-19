clc; clear; close all



fileName = {'unitCube100','unitCube500','unitCube1000'};

box = [-1,-1,-1;  ...
        1,-1,-1;  ...
        1, 1,-1;  ...
       -1, 1,-1;  ...
       -1,-1, 1;  ...
        1,-1, 1;  ...
        1, 1, 1;  ...
       -1, 1, 1];
edges = [1 2 3 4;5 6 7 8;3 4 8 7;1 2 6 5;2 3 7 6;1 4 8 5];

for i = 1:numel(fileName)
    load(fileName{i});
    figure();
    plot(f);
    ylim([min(f)/1.005, min(f)*1.05])
    xlabel('# of iterations')
    ylabel('F(x)')
    name = strcat(fileName{i},'F(x)');
    set(gca,'FontSize',15)
    saveas(gcf,name,'epsc')
    
    figure()
    semilogy(g)
    xlabel('# of iterations')
    ylabel('||\nabla F(x)||')
    name = strcat(fileName{i},'gradF(x)');
    set(gca,'FontSize',15)
    saveas(gcf,name,'epsc')

    % Plot grid
    Cbi = G.faces.neighbors(any(G.faces.neighbors==0,2),:);
    Cbi = Cbi(:);
    Cbi = Cbi(Cbi~=0);
    Cbi = unique(Cbi);
    
    Cb = false(G.cells.num,1);
    Ci = true(G.cells.num,1);
    Cb(Cbi) = true(size(Cbi));
    Ci(Cbi) = false(size(Cbi));

        
    G = computeGeometry(G);
    c = G.cells.centroids(:,1)>-0.5;

    Ci = Ci&c;
    Cb = Cb&c;
    figure()
    plotGrid(G,Ci,'facecolor','g');
    plotGrid(G,Cb);

    patch('vertices', box, 'faces', edges,'facealpha',0)
    view(-50,30)
    axis equal
    axis off
    
    name = strcat(fileName{i},'Grid');
    saveas(gcf, name, 'epsc')

end




