function B_n = idx_bus(N_b,n)
% determine index of bus where each prosumer is connected to
% W. Ananduta
% 02/10/2020

% Input: N_b, cell type, each cell represents a bus. It is the set of prosumers connected to the associated bus
% Output: B_n, vector with the size of the number of prosumers in the
% network

B_n = zeros(n,1);
for y = 1:length(N_b)
    for ii = 1:length(N_b{y})
        i = N_b{y}(ii);
        B_n(i) = y;
    end
end

end