function np = alg_param_37b1(np)
% Step-size selection
% W. Ananduta
% 28/09/2020

np.A = cell(np.n,1);
a_add = 1;
% alphas
for i=1:np.n
    Ni=sum(np.Adj(i,:));
    
    % alphas
    
    
    a_pi = 1 + a_add;
    a_st = 1 + a_add;
    a_mg = 3 + np.n*np.d_mg + a_add; 
    a_mgs = 3 + np.n*np.d_mgs + a_add; 
    a_tr =(a_add+Ni)*ones(Ni,1); 
    
    a_all = [a_pi;a_st;a_mg;a_mgs;a_tr];
    a_ik = diag(a_all);
    np.A{i}=kron(eye(np.h),a_ik);    % Structure/order different than what's written on the paper, due to dynamics of storage


    

end


% beta
b = 0.49;
np.beta_tr = b*ones(np.n);

for y=1:np.b
    
    a_the = a_add;
    a_v = a_add;
    a_tg = 2 + a_add;
    a_p = (1 + a_add)*ones(length(np.B{y}),1);
    a_q = (0 + a_add)*ones(length(np.B{y}),1);
    
    a_all = [a_the; a_v; a_tg; a_p; a_q];
    np.A_DSO{y} = kron(diag(a_all),eye(np.h));
    
    np.beta_pb(y) = 0.5/(length(np.B{y})+2*length(np.N_b{y})+a_add);
    

end

np.beta_tg = 1/(np.n+np.b+a_add);

np.gamma_mg = 1/(np.n+a_add);

end
