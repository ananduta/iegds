function p = alg_param3(p)
% Step-size selection
% W. Ananduta
% 12/11/2021


offset = 1;
% alphas
for i=1:p.n
    %Ni=sum(p.en.Adj(i,:));
    
    % alphas
    a_all = 2*ones(1,p.nu(i));
    %a_all = ones(1,p.nu(i));
    % x
%     a_x = 2*ones(1,p.nx(i));
%     a_x(p.id_mg) = (p.n*p.en.d_mg + p.n + offset);
%     a_x(p.id_gu) = (p.n*p.gn.d_gu + p.n + offset);
%     %a_x(1, 5:4+p.en.noN(i)) = (2 + offset)*p.en.noN(i)*ones(1,p.en.noN(i));
%     
    normPF = 0;
    normPFt = 0;
%     
    for jj = 1:p.en.noN(i)
        j = p.en.N{i}(jj);
        
        normPF = normPF + norm([p.m.PF{i,j}{1};p.m.PF{j,i}{2}]);
        normPFt = normPFt + norm([p.m.PF{i,j}{1};p.m.PF{j,i}{2}]');
        
        %a_x(p.id_eg+jj) = normPFt+offset;
    end
%     a_x(p.id_th) = normPFt+offset;
%     a_x(p.id_v) = normPFt+offset;
    
    a_x = max([normPFt+offset,p.n*p.en.d_mg + p.n + offset,p.n*p.gn.d_gu + p.n + offset])*ones(1,p.nx(i));
    % y
    %a_all(1, p.nx(i)+ 1) = offset;
    normR = 0;
    normG = 0;
    normRt = 0;
    normGt = 0;
    for jj = 1:p.gn.noN(i)
        j = p.gn.N{i}(jj);
        
        normR = normR + norm([p.m.Hc{i,j}{1};p.m.Hc{j,i}{2}]);
        normRt = normRt + norm([p.m.Hc{i,j}{1};p.m.Hc{j,i}{2}]');
        normG = normG + norm([p.m.G{i,j}{1};p.m.G{j,i}{2}]);
        normGt = normGt + norm([p.m.G{i,j}{1};p.m.G{j,i}{2}]');
        
    end
    a_y = zeros(1,p.ny(i));
    a_y(1) = 2*(normRt+normGt)+offset;
    a_y(2) = 1+offset;
    a_y(1,3:p.ny(i)-p.nt(i)) = kron(ones(1,p.gn.noN(i)),[(2+offset)*p.gn.noN(i) normRt normRt*ones(1,p.gn.r)]);
      
        
    % z
    a_z =(max(normRt,normGt) + p.n*offset)*ones(1,p.nz(i));
%    a_all(1,p.nx(i)+p.ny(i)+1:p.nu(i)) = max(normRt,normGt) + p.n*offset;
    
    a_all = [a_x a_y a_z];
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
    
    p.beta_pl(i) = 1/(normPF+offset);
    p.beta_gf(i) = 1/(normR+offset);

    p.beta_psi(i) = 1/(normG+offset);
end

p.delta = 1/(p.n+1);
p.delta_gu = 1/(p.n+1);
% beta
b = 2+offset;


p.beta_phi = 1/b*ones(p.n);




% gamma
gam = 0.1/(p.n);
p.gamma = gam*ones(p.n,1);
p.gamma_gu = gam*ones(p.n,1);
end
