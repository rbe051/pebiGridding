% Copyright (C) 2016 Runar Lie Berge. See COPYRIGHT.TXT for details.
close all; clear
%% Examples
% example1: Single fault intersected by several wells 
% example2: Complex wells intersecting
ex = 3;


%% Generate grid
switch ex
    case 1
        wellLine = {[0.6,0.2;0.65,0.6],...        
                    [0.3,0.3;0.7,0.8],...
                    [0.6,0.2;0.85,0.4],...
                    [0.15,0.7;0.4,0.7]};
        %wellLine ={[0.3,0.3;0.7,0.8]}; 

        x     = linspace(0.2, 0.8, 10);
        y     = 0.8 - 0.5*x - 0.05* sin(6*pi*x);
        fault = {[x' , y']};
        pri   = [2,3,1,4,5];
        Gp = compositePebiGrid(1/24, [1, 1], ...
                               'wellLines', wellLine, 'wellGridFactor', 0.5^2, ...
                               'faultLines',fault, 'faultGridFactor', 1/sqrt(2),...
                               'circleFactor', 0.6,'mlqtMaxLevel', 2, ...
                               'mlqtLevelSteps',[0.06,0.025]', 'priOrder', pri);
        Gdist = pebiGrid(1/24, [1, 1], 'wellLines', wellLine, ...
                        'wellGridFactor', 0.5^2, 'wellRefinement', true, ...
                        'epsilon',1/12, ...
                        'faultlines', fault, 'circleFactor', .6,...
                        'faultGridFactor', 1/sqrt(2),'priOrder', pri);
    case 2
        fault    = {};
        x        = linspace(0.2,0.8);
        wellLine = {[0.5,0.2; 0.5,0.3;0.47,0.4;0.4,0.5; 0.33,0.6;0.26,0.7], ...
                    [0.5,0.3;0.53,0.4;0.58,0.5],...            
                    [0.5,0.45;0.5,0.55;0.45,0.65;0.4,0.75;0.38,0.85],...
                    [0.5,0.55;0.55,0.65;0.6,0.75;0.62,0.85]};

        Gp = compositePebiGrid(1/14, [1, 1], 'wellLines', wellLine, ...
                              'wellGridFactor', 0.5^2, ...
                              'mlqtMaxLevel', 2, 'mlqtLevelSteps',[0.09,0.04]');
        Gdist = pebiGrid(1/14, [1, 1], 'wellLines', wellLine, ...
                        'wellGridFactor', 0.5^2, 'wellRefinement',true, 'epsilon',1/7);
    case 3
        [nx,ny,nz] = deal(20,10,10);
        [Lx,Ly,Lz] = deal(400,200,50);

        gridSize   = norm([Lx,Ly])/norm([nx,ny]);

        fault      = {[40,60;360,180], [40,160;250,40]};
        wellLine   = {[40,100], [350, 150], [160,140]};

        mlqtMax      = 2;                 % Set number of reminement levels
        wellGridSize = 0.8/2^mlqtMax;     % Set relative well grid size
        mlqtSizes = 2.0*linspace(gridSize,gridSize*wellGridSize,mlqtMax+1)';
                                          % Set distance around wells to be
                                          % refined.

        wellEps = sqrt(gridSize)*8;       % Size around wells to be refined
                
        Gp = compositePebiGrid(gridSize, [Lx, Ly], 'wellLines', wellLine, ...
                               'wellGridFactor', wellGridSize,...
                               'mlqtMaxLevel', 2, 'faultLines', fault);

        Gdist = pebiGrid(gridSize, [Lx, Ly], 'wellLines', wellLine, ...
                        'wellRefinement', true, 'wellGridFactor', wellGridSize,...
                        'epsilon',50, 'faultLines', fault);
                    
        % Create 3D grid
        Gp3D    = makeLayeredGrid(Gp, nz);
        Gp3D.nodes.coords(:,3) = Gp3D.nodes.coords(:,3)*Lz/nz;
        Gp3D    = computeGeometry(Gp3D);
        Gdist3D = makeLayeredGrid(Gdist, 5);
        Gdist3D.nodes.coords(:,3) = Gdist3D.nodes.coords(:,3)*Lz/nz;
        Gdist3D = computeGeometry(Gdist3D);

    otherwise
        error('Unknown Example')
end
                  
                  
              
%% Plotting                       
orange = [1,138/255,0.1];      
figure()
hold on
plotGrid(Gp, 'faceColor', 'none')
axis equal tight off
hold on
for i = 1:numel(wellLine)
  line = wellLine{i};
  if size(line,1) == 1
      plot(line(1,1), line(1,2),'.r', 'markersize', 8);
  end
  plot(line(:, 1), line(:, 2),'r');
end
for i = 1:numel(fault)
  line = fault{i};
  plot(line(:, 1), line(:, 2),'color',orange);
end

figure()
hold on
plotGrid(Gdist, 'faceColor', 'none')
axis equal tight off
hold on
for i = 1:numel(wellLine)
  line = wellLine{i};
  if size(line,1) == 1
      plot(line(1,1), line(1,2),'.r', 'markersize', 8);
  end
  plot(line(:, 1), line(:, 2),'r');
end
for i = 1:numel(fault)
  line = fault{i};
  plot(line(:, 1), line(:, 2),'color', orange);
end

    
if ex == 3
    figure()
    hold on
    plotGrid(Gp3D)
    axis equal tight off
    hold on
    for i = 1:numel(wellLine)
      line = wellLine{i};
      if size(line,1) == 1
          plot3(line(1,1), line(1,2), -0.01*ones(size(line,1),1),'.r', 'markersize', 8);
      end
      plot3(line(:, 1), line(:, 2),-0.01*ones(size(line,1),1),'r');
    end
    for i = 1:numel(fault)
      line = fault{i};
      plot3(line(:, 1), line(:, 2),-0.01*ones(size(line,1),1),'color',orange);
    end

    figure()
    hold on
    plotGrid(Gdist3D)
    axis equal tight off
    hold on
    for i = 1:numel(wellLine)
      line = wellLine{i};
      if size(line,1) == 1
          plot3(line(1,1), line(1,2), -0.01*ones(size(line,1),1),'.r', 'markersize', 8);
      end
      plot(line(:, 1), line(:, 2),'r');
    end
    for i = 1:numel(fault)
      line = fault{i};
      plot3(line(:, 1), line(:, 2),-0.01*ones(size(line,1),1),'color',orange);
    end
    
end
