function [s,o] = com_cost1(s,np)
    s.t = length(s.res1);
    
    for i = 1:np.n
        Qh = np.Qh{i};
        ch = np.ch{i};
        z = s.u{i}(:,s.t);
        
        f_mg = np.d_mg*(s.sigma_mg(:,s.t)')*np.Smg{i}*z;
        f_mgs = np.d_mgs*(s.sigma_mg(:,s.t)')*np.Smgs{i}*z;
        s.J(i,1) = z'*Qh*z + ch'*z + f_mg + f_mgs;
        s.J_mg(i,1) = f_mg;
        s.J_mgs(i,1) = f_mgs;
        
        % save output
        o.u{i}= z;
        o.J(i,1) = s.J(i,1);
        o.J_mg(i,1) = f_mg;
        o.J_mgs(i,1) = f_mgs;
    end
    
    for y=1:np.b
        o.u_DSO{y} = s.u_DSO{y}(:,s.t);
    end
    
    s.Jt = sum(s.J);
    s.Jt_mg= sum(s.J_mg);
    s.J_pd = np.d_mg*(s.sigma_mg(:,s.t)')*s.sum_p_pd;
    
    o.Jt = s.Jt;
    o.Jt_mg = s.Jt_mg;
    o.J_pd = s.J_pd;
end