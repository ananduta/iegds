function np = alg_param_fbf(np)
% Step-size selection for FBF
% W. Ananduta
% 15/07/2021


% Lipschitz constant of pseudogradient
L_f = max([2*max(np.q_dg), 2*max(np.q_st),2*np.d_mg]);

% Lipschitz constant of C
for i=1:np.n
    noNi(i) = length(np.N{i});
end
L_c = max([2,2*max(noNi)]);

L_B = max(L_f, L_c);

ss_ref = 1/(L_B+0.1); %% COMPUTE BASED ON LIPSCHITZ CONSTANTS!

% rho
np.rho = ss_ref*ones(np.n,1);

% tau_tr
np.tau_tr = 0.4*ones(np.n);

% gamma
np.gamma = ss_ref*ones(np.n,1);

np.gamma_gu = ss_ref*ones(np.n,1);

% tau_mg
np.tau_mg = ss_ref*ones(np.n,1);

np.tau_gu = ss_ref*ones(np.n,1);

% dummy var
np.A = cell(np.n,1);
for i=1:np.n
    Ni=sum(np.Adj(i,:));
    np.A{i} = zeros((4+Ni)*np.h);
end
end
