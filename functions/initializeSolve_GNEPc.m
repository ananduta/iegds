function [p,s] = initializeSolve_GNEPc(p,r)

    %% INITIALIZATION

    

    

    % Generate matrices for  constraints       
    p = build_mat_exP2P_pen(p);  
    
    % Assigning parameters of the algorithm
    p = alg_param3(p);     % CHECK THIS
    
    
    % Generate matrices for cost function
    p = build_mat_cost_pen(p);
    
    % Initialization of the variables
    k=1;
    s.sigma_mg = p.en.sumPd(1:p.h);
    s.sigma_gu = p.gn.sumGd(1:p.h);
    
    for i=1:p.n
            % decision variables
            %s = initialize_u_quadp(p,s,i);
           if r == 1
                p.init = 0;
                s.u{i}(:,1) = p.init*ones(size(p.m.Aineq{i},2),1);
           else
               p.init = rand;
               s.u{i}(:,1) = p.init*ones(size(p.m.Aineq{i},2),1);
          end
            % x
            s.p_di{i}(:,1) = p.m.Sdg{i}*s.u{i}(:,1);
            s.p_ch{i}(:,1) = p.m.Sch{i}*s.u{i}(:,1);
            s.p_dh{i}(:,1) = p.m.Sdh{i}*s.u{i}(:,1);
            s.p_mg{i}(:,1) = p.m.Smg{i}*s.u{i}(:,1);
            s.d_gu{i}(:,1) = p.m.Sgu{i}*s.u{i}(:,1);
            for jj = 1:length(p.en.N{i})
                j = p.en.N{i}(jj);
                s.p_l{i,j}(:,1) = p.m.Spl{i,j}*s.u{i}(:,1);
            end

            s.sigma_mg(:,1) = s.sigma_mg(:,1) + s.p_mg{i}(:,1);
            s.sigma_gu(:,1) = s.sigma_gu(:,1) + s.d_gu{i}(:,1);
            
            
            % y
            s.psi{i}(:,1) = p.m.Spsi{i}*s.u{i}(:,1);
            for jj = 1:length(p.gn.N{i})
                j = p.gn.N{i}(jj);
                s.phi{i,j}(:,1) = p.m.Sphi{i,j}*s.u{i}(:,1);
                s.pen{i,j}(:,1) = p.m.Spen{i,j}*s.u{i}(:,1);
            end
            
            
            
            
    end
        % Dual variables
        for i=1:p.n
            for jj=1:length(p.en.N{i})
                j = p.en.N{i}(jj);
                sl.c_pl{i,j}(:,1) = p.m.PF{i,j}{1}*s.u{i}(:,1) + p.m.PF{i,j}{2}*s.u{j}(:,1);
                
                s.mu_pl{i,j}(:,1) = p.init*ones(p.h,1);
            end
            for jj = 1:length(p.gn.N{i})
                j = p.gn.N{i}(jj);
                
                sl.c_phi{i,j}(:,1) = s.phi{i,j}(:,1) + s.phi{j,i}(:,1);
                
                s.mu_phi{i,j}(:,1) = p.init*ones(p.h,1);
                
                
                sl.c_gf{i,j}(:,1) = p.m.Hc{i,j}{1}*s.u{i}(:,1) + p.m.Hc{i,j}{2}*s.u{j}(:,1);
                
                s.mu_gf{i,j}(:,1) = p.init*ones(p.h,1);
                
                sl.c_psi{i,j}(:,1) = p.m.G{i,j}{1}*s.u{i}(:,1) + p.m.G{i,j}{2}*s.u{j}(:,1) - p.m.g{i,j};
                
                s.lambda_psi{i,j}(:,1) = p.init*ones(2*p.h,1);
                
            end
        end
        % auxiliary var nu and lambda
            
            s.lambda_mg= zeros(2*p.h,1);
            s.lambda_gu= zeros(2*p.h,1);
end