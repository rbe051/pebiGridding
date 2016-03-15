function [] = plotFault(Gp, varargin)

    faultNodes = find(Gp.faces.isFault);
    numFaultLine = length(faultNodes);

    faultsLine = {};
    j = 1;

    if numFaultLine < 1
        return
    end 

    line  = Gp.faces.nodes(Gp.faces.nodePos(faultNodes(1)) ...
                          :Gp.faces.nodePos(faultNodes(1)+1)-1);

    for i = 2:numFaultLine
        lineSegment = Gp.faces.nodes(Gp.faces.nodePos(faultNodes(i)) ...
                                    :Gp.faces.nodePos(faultNodes(i)+1)-1);
        if lineSegment(1) == line(end);
             line = [line; lineSegment(2)];
        else

            faultsLine{j} = line;
            line = lineSegment;
            j = j +1;
        end
    end
    faultsLine{j} = line;


    hold on
    for i = 1:numel(faultsLine)
        faultCoords = Gp.nodes.coords(faultsLine{i},:);
        plot(faultCoords(:,1),faultCoords(:,2), varargin{:});
    end
end
