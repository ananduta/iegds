function nab_phi_i = grad_sel_hsdm(y,np,i)
% Generate gradient of selection function
% W. Ananduta
% 15/07/2021

Ni=length(np.N{i});

g = ones(4+Ni,1);
g(1,1) = 0;
g(2,1) = 0;
g(4,1) = 0;

c = zeros(4+Ni,1);
c(3,1) = np.pmg_ref;

gh = kron(ones(np.h,1),g);
ch = kron(ones(np.h,1),c);

nab_phi_i = gh.*y - ch;




end