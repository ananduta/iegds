function np = alg_param_D(np)
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
    
    a_tr =(a_add+Ni)*ones(Ni,1); 
    
    a_all = [a_pi;a_st;a_mg;a_tr];
    a_ik = diag(a_all);
    np.A{i}=kron(eye(np.h),a_ik);    % Structure/order different than what's written on the paper, due to dynamics of storage


    

end


% beta
b = 0.49;
np.beta_tr = b*ones(np.n);

for y=1:np.b
    
    a_the =  1/(2*(a_add+sum(np.Bnet(y,:))+sum(np.Gnet(y,:)))); 
    a_v = 1/(2*(a_add+ sum(np.Bnet(y,:))-sum(np.Gnet(y,:))));    
    %a_the = 2 + a_add;
    %a_v = 2 + a_add;
    a_tg = 2 + a_add;
    a_p = (2 + a_add)*ones(length(np.B{y}),1);
    a_q = (1 + a_add)*ones(length(np.B{y}),1);
    
    a_all = [a_the; a_v; a_tg; a_p; a_q];
    np.A_DSO{y} = kron(diag(a_all),eye(np.h));
    
    np.beta_pb(y) = 0.2/(length(np.B{y})+2*length(np.N_b{y})+length(np.B{y})+a_add);
    
    np.gamma_mg(y) = 0.1/(length(np.N_b{y})+length(np.B{y})+a_add);
    
    np.beta_tg(y) = 0.1/(1+2*length(np.N_b{y})+length(np.B{y})+a_add);
    
    for zz = 1:length(np.B{y})
        z = np.B{y}(zz);
        np.beta_p(y,z) = 0.1/(2*abs(np.Bnet(y,z))+2*abs(np.Gnet(y,z))+2+a_add);
        np.beta_q(y,z) = 0.1/(2*abs(np.Bnet(y,z))+2*abs(np.Gnet(y,z))+2+a_add);
        
    end
end
np.beta_pc = 0.1;
np.beta_qc = 0.1;




np.delta_mg = 0.49;
np.delta_tg = 0.49;

end
