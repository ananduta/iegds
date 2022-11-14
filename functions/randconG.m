function A = randconG(n,sp)
% Generate a random connected graph
% W. Ananduta
% 19/07/2019

% Input: number of nodes, sparsity sp in (0,1)
% Output: Adjacency matrix of the graph


% number of edges
n_E = floor(sp*n*(n-1)/2);

% break conditions
if n_E < n-1
    disp('generation failed. Increase sparsity');
    return
end
if sp > 1
    disp('sparsity cannot be larger than 1');
    return
end
if sp < 0
    disp('sparsity cannot be smaller than 0');
    return
end

N = 1:n;
A = [zeros(n-1,1), eye(n-1);
     0, zeros(1,n-1)];
A = A + A';

c_ad=n_E-(n-1); % counter of additional edges

while c_ad > 0
    i = randi(n);
    j= randi(n);
    if i==j
        continue;
    end
    if A(i,j) == 1
        continue;
    end
    A(i,j) =1;
    A(j,i) =1;
    c_ad = c_ad-1;
end

end

