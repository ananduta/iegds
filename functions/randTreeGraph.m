function Adj = randTreeGraph(noNodes,maxBranch_perNode)
% W. Ananduta
% 27/01/2022
% Generate a random tree graph given number of nodes and maximum branches per node

    %maxBranch_perNode = 3;
    
    Adj = zeros(noNodes);
    set_assignedNodes = [];
    
    % node 1 is always the root
    set_assignedNodes = [1];
    for i = 1:noNodes
        if find(set_assignedNodes==i)
            noBranch = randi(maxBranch_perNode);
            countBranch = 0;
            for j=1:noNodes
                if isempty(find(set_assignedNodes==j)) && countBranch~=noBranch
                    Adj(i,j) = 1;
                    countBranch = countBranch + 1;
                    set_assignedNodes = [set_assignedNodes; j];
                end
            end
        end
    end
    Adj = Adj + Adj';
end