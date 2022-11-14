function nabx_L = gradx_L(np,u,sigma_mg,lambda_mg,sigma_gu,lambda_gu,mu_tr,t,i)

    nabJ = 2*np.Qh{i}*u + np.ch{i};
                
    pmg_i_t = np.Smg{i}*u;
    dgu_i_t = np.Sgu{i}*u;
    c1 = [];
    for h = 1:np.h
        %c_1 =  [0; 0; np.q_mg*(s.ph_mg(h,t)+pmg_i_t(h)); zeros(sum(np.Adj(i,:)),1)]; 
        c_1 =  [0; 0; np.d_mg*(sigma_mg(h,t)-pmg_i_t(h)); np.d_gu*(sigma_gu(h,t)-dgu_i_t(h)); zeros(sum(np.Adj(i,:)),1)]; 
        c1 = [c1;c_1];
    end
    

    c4 = [np.Smg{i};-np.Smg{i}]'*lambda_mg{i}(:,t);
    
    c6 = [np.Sgu{i};-np.Sgu{i}]'*lambda_gu{i}(:,t);

    c5 = zeros(np.h*(4+sum(np.Adj(i,:))),1);
    for jj=1:length(np.N{i})
        j = np.N{i}(jj);
        c5 = c5 + np.Str{i,j}'*mu_tr{i,j}(:,t);           
    end

    
    nabx_L = nabJ + c1 + c4 + c5 + c6;
    
    
end