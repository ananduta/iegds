function p = alg_param2(p)
% Step-size selection
% W. Ananduta
% 24/07/2019


offset = 1;
% alphas
for i=1:p.n
    %Ni=sum(p.en.Adj(i,:));
    
    % alphas
    a_all = p.gn.r*p.n*10*ones(1,p.nu(i));
    %a_all = ones(1,p.nu(i));
    % x
    a_all(3) = (p.n*p.en.d_mg + p.n + offset);
    a_all(4) = (p.n*p.gn.d_gu + p.n + offset);
    a_all(1, 5:4+p.en.noN(i)) = (2 + offset)*ones(1,p.en.noN(i));
    
    % y
    a_all(1, p.nx(i)+ 1) = 2*p.gn.noN(i) + offset;
    for jj = 1:p.gn.noN(i)
        j = p.gn.N{i}(jj);
        id_nei_s = 2 + (2+p.gn.r)*(jj-1) + 1;
        id_nei_e = 2+(2+p.gn.r)*jj;
        
        % define the parameters of the affine function at each region
        nf = @(y) y^2/p.gn.cf(i,j);
        pwaf = pwa_approx_nf(p.r,-p.gn.phi_max(i,j),p.gn.phi_max(i,j),nf);
        
        a_all(1, id_nei_s: id_nei_e) = [p.gn.r*p.n*10 p.gn.r*p.n*10 2*abs(pwaf.a')] + offset*ones(1,2+p.gn.r);
        
    end
    % aux penalty term
    a_all(1,p.nx(i)+p.ny(i)-p.nt(i)+1:p.nx(i)+p.ny(i)) = zeros(1,p.nt(i));
    
    % z
    for jj = 1:p.gn.noN(i)
        j = p.gn.N{i}(jj);
        id_delta_psi = p.nx(i)+p.ny(i)+(1+3*p.gn.r)*(jj-1) + 1;
        
        a_all(id_delta_psi) = 2 + offset;
        
        % define the parameters of the affine function at each region
        nf = @(y) y^2/p.gn.cf(i,j);
        pwaf = pwa_approx_nf(p.r,-p.gn.phi_max(i,j),p.gn.phi_max(i,j),nf);
        
        for m =1:p.gn.r
            id_delta_m = p.nx(i)+p.ny(i)+(1+3*p.gn.r)*(jj-1) + 1+3*(m-1) + 3;
            
            a_all(id_delta_m) = 2*abs(pwaf.b(m))*p.gn.r + offset;
        end
        
    end
    
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
b = 2+offset;
p.beta_tr = 1/b*ones(p.n);

p.beta_phi = 1/b*0.1*ones(p.n);

nf = @(y) y^2/max(max(p.gn.cf));
pwaf = pwa_approx_nf(p.gn.r,-max(max(p.gn.phi_max)),max(max(p.gn.phi_max)),nf);
c = max(max(abs(pwaf.a)),max(abs(pwaf.b)));
p.beta_gf = 0.1/(b*c*p.gn.r)*ones(p.n);

p.beta_psi = 1/b*ones(p.n)*0.2;

% gamma
gam = 0.1/(p.n);
p.gamma = gam*ones(p.n,1);
p.gamma_gu = gam*ones(p.n,1);
end
