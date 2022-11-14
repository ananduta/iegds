function np = alg_param_pfb(np)
% Step-size selection for pFB
% W. Ananduta
% 15/07/2021

c = 0.01;
eta = 1/max([2*max(np.q_dg), 2*max(np.q_st),2*np.d_mg]);
for i=1:np.n
    noNi(i) = length(np.N{i});
end
delta = 1/min([eta,  1/(2*max(noNi))])+0.1;

% rho
np.rho = (1/delta)*c*ones(np.n,1);

% tau_tr
np.tau_tr = 1/(1+2+delta)*c*ones(np.n);


% gamma
%np.gamma = ss_ref*zeros(np.n,1);
%np.tau_mg = zeros(np.n,1);
for i = 1:np.n
    % tau_mg
    np.tau_mg(i) = 1/(1+2*noNi(i)+delta)*c;
    np.gamma(i) = 1/(2*noNi(i)+delta)*c;
end

% dummy var
np.A = cell(np.n,1);
for i=1:np.n
    Ni=sum(np.Adj(i,:));
    np.A{i} = zeros((3+Ni)*np.h);
end
end
