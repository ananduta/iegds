function p = alg_param(p)
% Step-size selection
% W. Ananduta
% 24/07/2019


% alphas
for i=1:p.n
    %Ni=sum(p.en.Adj(i,:));
    
    % alphas
    a_all = ones(1,p.nu(i));
    a_all(3) = 2*(p.n*p.en.d_mg + 2);
    a_all(4) = 2*(p.n*p.gn.d_gu + 2);
    a_all(1, 5:4+p.en.noN(i)) = 2*p.en.noN(i)*ones(1,p.en.noN(i));
    a_all(1, p.nx(i)+ 1) = 3;
    a_all(1, p.nx(i)+ 3: p.nu(i)) = p.gn.r*p.gn.noN(i)*ones(1,  p.nu(i)-p.nx(i)- 2);
    p.Alpha{i} = kron(eye(p.h),diag(a_all));
    %a_dg = 1;
    %a_st = 1;
    %a_mg = 2*(np.n*np.d_mg + 2);
    %a_gu = 2*(np.n*np.d_gu + 2);
    %a_tr =2*Ni*ones(Ni,1);
    %a_all = [a_dg;a_st;a_mg;a_gu;a_tr];
    %a_ik = diag(a_all);
    %np.A{i}=kron(eye(np.h),a_ik);
    %p.alpha(i) = 1;%(p.nx(i)+p.n*max(p.en.d_mg,p.gn.d_gu));

    % delta
%    np.delta(i,1) = 0.8/(Ni);
    
end
p.delta = 1/(p.n+1);
    p.delta_gu = 1/(p.n+1);
% beta
b = 0.49;
p.beta_tr = b*ones(p.n);

p.beta_phi = b*ones(p.n);

p.beta_gf = b*0.02/p.gn.r*ones(p.n);

p.beta_psi = b*0.1*ones(p.n);

% gamma
gam = 0.1/(p.n);
p.gamma = gam*ones(p.n,1);
p.gamma_gu = gam*ones(p.n,1);
end
