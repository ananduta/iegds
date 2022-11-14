function [s,np] = pfb_hsdm(np)
%% Distributed algorithm (pFB + HSDM) for GNE selection
% P2P market
% W. Ananduta
% 15/07/2021

        %% INITIALIZATION
        % Assigning parameters of the algorithm
        np = alg_param_pfb(np);
        
        
        np.t_max = 1000;
        np.er_max = 1e-1; 

        % Generate matrices for cost function and constraints       
        np = build_mat_exP2P(np);
        
        
        
        %% Initialization of the variables
        
        t=1;
        s.sigma_mg = np.sumPd;
        s.sigma_gu = np.sumGd;
        %s.sigma_mg = zeros(np.h,1);
        for i=1:np.n
            % decision variables
            %s = initialize_u_quadp(np,s,i);
            s.u{i}(:,1) = np.init*ones(size(np.A_ineq{i},2),1);
            s.p_di{i}(:,1) = np.Sdi{i}*s.u{i}(:,1);
            s.p_st{i}(:,1) = np.Sst{i}*s.u{i}(:,1);
            s.p_mg{i}(:,1) = np.Smg{i}*s.u{i}(:,1);
            s.d_gu{i}(:,1) = np.Sgu{i}*s.u{i}(:,1);
            for jj = 1:length(np.N{i})
                j = np.N{i}(jj);
                s.p_tr{i,j} = np.Str{i,j}*s.u{i}(:,1);
            end
            
            s.sigma_mg(:,1) = s.sigma_mg(:,1) + s.p_mg{i}(:,1);
            s.sigma_gu(:,1) = s.sigma_gu(:,1) + s.d_gu{i}(:,1);
            % auxiliary var nu and lambda
            s.nu{i} = zeros(2*np.h,1);
            s.lambda_mg{i}= zeros(2*np.h,1);
            s.nu_gu{i} = zeros(2*np.h,1);
            s.lambda_gu{i}= zeros(2*np.h,1);
        end
        for i=1:np.n
            for jj=1:length(np.N{i})
                j = np.N{i}(jj);
                sl.c_tr{i,j}(:,1) =  s.p_tr{i,j}(:,1) + s.p_tr{j,i}(:,1);              
                s.mu_tr{i,j}(:,1) = np.init*ones(np.h,1);
            end
        end
        
        
        


        %% Iteration
        while 1
            tic
            t
            

            %s.ph_mg(:,t+1) = 0;
            s.sigma_mg_t(:,t) = np.sumPd;%s.sum_p_pd;
            s.sigma_gu_t(:,t) = np.sumGd;
            %s.sigma_mg(:,t+1) = zeros(np.h,1);
            
           for i=1:np.n
                
                % primal update of prosumer (quadratic prog.)
                %s = loc_opt_qprog(np,s,sl,t,i);
                
                nabx_L = gradx_L(np,s.u{i}(:,t),s.sigma_mg,s.lambda_mg,s.sigma_gu,s.lambda_gu,s.mu_tr,t,i);
                y = s.u{i}(:,t) - np.rho(i)*nabx_L;

                u_i = projXi(y,np,i);
                s.u{i}(:,t+1) = u_i;
                s.p_di{i}(:,t+1) = np.Sdi{i}*s.u{i}(:,t+1);
                s.p_st{i}(:,t+1) = np.Sst{i}*s.u{i}(:,t+1);
                s.p_mg{i}(:,t+1) = np.Smg{i}*s.u{i}(:,t+1);
                s.d_gu{i}(:,t+1) = np.Sgu{i}*s.u{i}(:,t+1);
                for jj = 1:length(np.N{i})
                    j = np.N{i}(jj);
                    s.p_tr{i,j}(:,t+1) = np.Str{i,j}*s.u{i}(:,t+1);
                end
                %s.sigma_mg_t(:,t+1) = s.sigma_mg_t(:,t) + s.p_mg_t{i}(:,t);
                
                s.nu{i}(:,t+1) = s.nu{i}(:,t);
                s.nu_gu{i}(:,t+1) = s.nu_gu{i}(:,t);
                % dual variables (reciprocity constraint) Prosumers
                for jj=1:length(np.N{i})
                        j = np.N{i}(jj);

                        
                        % Auxiliary variable update:
                        s.nu{i}(:,t+1) = s.nu{i}(:,t+1) - np.gamma(i)*(s.lambda_mg{i}(:,t)-s.lambda_mg{j}(:,t));
                        
                        s.nu_gu{i}(:,t+1) = s.nu_gu{i}(:,t+1) - np.gamma_gu(i)*(s.lambda_gu{i}(:,t)-s.lambda_gu{j}(:,t));

                   
                end
                
                
                
            end

            res2(t) = 0; % Extra: compute norm of residual
            dres2(t) = 0;

            s.sigma_mg(:,t+1) = np.sumPd;
            s.sigma_gu(:,t+1) = np.sumGd;
            for i=1:np.n
                
                s.lambda_mg{i}(:,t+1) = s.lambda_mg{i}(:,t);
                s.lambda_gu{i}(:,t+1) = s.lambda_gu{i}(:,t);
                % dual variables (reciprocity constraint) Prosumers
                for jj=1:length(np.N{i})
                        j = np.N{i}(jj);

                        % Dual variables (reciprocity constraints) update:
                        
                        % reciprocity constraint
                        sl.c_tr{i,j}(:,t+1) = s.p_tr{i,j}(:,t+1) + s.p_tr{j,i}(:,t+1);
                        s.mu_tr{i,j}(:,t+1) = s.mu_tr{i,j}(:,t) + np.tau_tr(i,j)*(2*sl.c_tr{i,j}(:,t+1)-sl.c_tr{i,j}(:,t));
                                              
                        % grid constraint update:
                        s.lambda_mg{i}(:,t+1) = s.lambda_mg{i}(:,t+1) +  np.tau_mg(i)*(2*s.nu{i}(:,t+1)-2*s.nu {j}(:,t+1)...
                                                   -(s.nu{i}(:,t)-s.nu{j}(:,t)-s.lambda_mg{i}(:,t) + s.lambda_mg{j}(:,t)) );
                        
                        % gas load constraint update:
                        s.lambda_gu{i}(:,t+1) = s.lambda_gu{i}(:,t+1) +  np.tau_gu(i)*(2*s.nu_gu{i}(:,t+1)-2*s.nu_gu{j}(:,t+1)...
                                                   -(s.nu_gu{i}(:,t)-s.nu_gu{j}(:,t)-s.lambda_gu{i}(:,t) + s.lambda_gu{j}(:,t)) );
                        
                        % Extra: compute norm of residual
                        res = sl.c_tr{i,j}(:,t+1);
                        res2(t) = norm([res2(t);res],inf);

                   
                end
                
                % dual variable (grid constraint) update (cont'd):
                
                s.lambda_mg{i}(:,t+1) = max(0, s.lambda_mg{i}(:,t+1) + np.tau_mg(i)*(2*[s.p_mg{i}(:,t+1);-s.p_mg{i}(:,t+1)] ...
                                            -[s.p_mg{i}(:,t);-s.p_mg{i}(:,t)] - [np.pmg_max/np.n*ones(np.h,1);-np.pmg_min/np.n*ones(np.h,1)]));
                
                s.lambda_gu{i}(:,t+1) = max(0, s.lambda_gu{i}(:,t+1) + np.tau_gu(i)*(2*[s.d_gu{i}(:,t+1);-s.d_gu{i}(:,t+1)] ...
                                            -[s.d_gu{i}(:,t);-s.d_gu{i}(:,t)] - [np.dg_max/np.n*ones(np.h,1);-np.dg_min/np.n*ones(np.h,1)]));
                
%                 % HSDM step
%                 beta_hsdm = 1/t;
%                 nab_phi_i = grad_sel_hsdm(s.u{i}(:,t+1),np,i);
%                 s.u{i}(:,t+1) = s.u{i}(:,t+1) - beta_hsdm*nab_phi_i;
%                 s.p_di{i}(:,t+1) = np.Sdi{i}*s.u{i}(:,t+1);
%                 s.p_st{i}(:,t+1) = np.Sst{i}*s.u{i}(:,t+1);
%                 s.p_mg{i}(:,t+1) = np.Smg{i}*s.u{i}(:,t+1);
%                 s.d_gu{i}(:,t+1) = np.Sgu{i}*s.u{i}(:,t+1);
%                 for jj = 1:length(np.N{i})
%                     j = np.N{i}(jj);
%                     s.p_tr{i,j}(:,t+1) = np.Str{i,j}*s.u{i}(:,t+1);
%                 end
                
                s.sigma_mg(:,t+1) = s.sigma_mg(:,t+1) + s.p_mg{i}(:,t+1);
                s.sigma_gu(:,t+1) = s.sigma_gu(:,t+1) + s.d_gu{i}(:,t+1);
                % compute error

                res1(i,t+1)= norm(s.u{i}(:,t+1)-s.u{i}(:,t),inf);
            end

            
            
            
            % Extra: stopping criterion
            s.res(t) = (res2(t));
            
            %s.dres(t) = sqrt(dres2(t));
            
            
            
            s.res1(t) = norm(res1(:,t+1),inf);
            [s.res(t);s.res1(t)]
            
           
            if  s.res(t) < np.er_max && s.res1(t) < np.er_max
                break
            end
            

            t= t+1;
            toc
        end
end