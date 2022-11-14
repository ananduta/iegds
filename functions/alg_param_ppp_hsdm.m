function np = alg_param_ppp_hsdm(np)
% Step-size selection
% W. Ananduta
% 24/07/2019

np.a_bar = cell(np.n,1);
% alphas
for i=1:np.n
    Ni=sum(np.Adj(i,:));
    
    % alphas
    a_dg = 1;
    a_st = 1;
    a_mg = 2*(np.n*np.d_mg + 2);
    a_tr =2*Ni*ones(Ni,1);
    a_all = [a_dg;a_st;a_mg;a_tr];
    a_ik = diag(a_all);
    np.A{i}=kron(eye(np.h),a_ik);


    % delta
%    np.delta(i,1) = 0.8/(Ni);
    np.delta(i,1) = 1/(np.n);
end

% beta
b = 0.49;
np.beta_tr = b*ones(np.n);

% gamma
gam = 0.49;
np.gamma = gam*ones(np.n,1);



end
