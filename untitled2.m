%%
for i = 1:p.gn.n
    for jj = 1:p.gn.noN(i)
        j = p.gn.N{i}(jj);
        
        deltaFlowDir{i,j} = p.alphaFixed{i,j} - oSOCP.alpha_s{i,j};
        
        deltaFlow{i,j} = o{1}.phi{i,j} - oSOCP.phi{i,j};
    end
end