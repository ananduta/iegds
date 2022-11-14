function [s,sl,np] = dist_alg_eP2P(np)
%% Fully Distributed algorithm
% extended P2P market
% G. Belgioioso & W. Ananduta
% 12/11/2020

        %% INITIALIZATION
        % Assigning parameters of the algorithm
        %np = alg_param(np);
        np = alg_param_D(np);
        np.t_max = 20000;
        np.er_max = 0.1; 

        % % Generate matrices S_i^mg and S_ij^tr
        %[np.Sdg,np.Sst,np.Smg,np.Str] = gen_Smat(np);

        % Generate matrices for cost function and constraints
        np = build_mat_exP2P(np);                                             
        % np = build_mat_exDSO(np); %EDIT
        
        %% Initialization of the variables
        %s.ph_mg(:,1) = zeros(np.h,1);
        
        %s.sigma_mg = zeros(np.h,1);
        s=struct;
        for i=1:np.n
            % decision variables
            s = initialize_u_quadp(np,s,i);

            sl.b{i}(:,1) = np.Pd(i,1:np.h)' - s.p_di{i}(:,1) - s.p_st{i}(:,1);

            %s.sigma_mg(:,1) = s.sigma_mg(:,1) + s.p_mg{i}(:,1);
        end
        

        for i=1:np.n
            for jj=1:length(np.N{i})
                j = np.N{i}(jj);
                sl.c_tr{i,j}(:,1) = s.p_tr{i,j}(:,1) + s.p_tr{j,i}(:,1);
                
                sl.mu_tr{i,j}(:,1) = zeros(np.h,1);
            end
        end

        s.sum_p_pd = zeros(np.h,1);
        s.sum_sigma(:,1) = zeros(np.h,1);
        t = 0;
        for y=1:np.b
            
            % Primal variable
            s = projUy(s,0,np,y,t); 

            % aggregation update 
            s.sigma_mg{y}(:,t+1) = np.Pd(np.n+np.B_n(i),1:np.h)';
            sl.by{y}(:,t+1) = zeros(np.h,1);
            for ii = 1:length(np.N_b{y})
                i = np.N_b{y}(ii);
                s.sigma_mg{y}(:,t+1) = s.sigma_mg{y}(:,t+1) + s.p_mg{i}(:,t+1);
                sl.by{y}(:,t+1) = sl.by{y}(:,t+1) + sl.b{i}(:,t+1);
            end

            %auxiliary variables

            sl.w_mg{y}(:,t+1) = zeros(2*np.h,1);
            sl.w_tg{y}(:,t+1) = zeros(np.h,1);

            % dual update (global grid constraints)
            sl.c_mg{y}(:,t+1) = kron([1;-1], s.sigma_mg{y}(:,t+1)) -1/np.b*[np.pmg_max*ones(np.h,1);-np.pmg_min*ones(np.h,1)] + sl.w_mg{y}(:,t+1);
            sl.lambda_mg{y}(:,t+1) = zeros(2*np.h,1);
            sl.c_tg{y}(:,t+1) = s.sigma_mg{y}(:,t+1) - s.p_tg{y}(:,t+1) + sl.w_tg{y}(:,t+1);
            sl.mu_tg{y}(:,t+1) = zeros(np.h,1);

            % dual update (power balance on bus y)
            sl.c_pb{y}(:,t+1) = np.Pd(np.n+y,1:np.h)' - s.p_tg{y}(:,t+1)+sl.by{y}(:,t+1);
            for zz = 1:length(np.B{y})
                z = np.B{y}(zz);
                sl.c_pb{y}(:,t+1) = sl.c_pb{y}(:,t+1) + s.p_l{y,z}(:,t+1);
            end
            sl.mu_pb{y}(:,t+1) = zeros(np.h,1);
            
            s.sum_p_pd = s.sum_p_pd + np.Pd(np.n+y,1:np.h)';
            s.sum_sigma(:,t+1) = s.sum_sigma(:,t+1) + s.sigma_mg{y}(:,t+1);
        end
            
                % dual update (power flow equation)
        for y = 1:np.b
            for zz = 1:length(np.B{y})
                z = np.B{y}(zz);

                bp = s.p_l{y,z}(:,t+1) - s.p_l{z,y}(:,t+1);
                sl.c_p{y,z}(:,t+1) = 2*np.Bnet(y,z)*(s.theta{y}(:,t+1)-s.theta{z}(:,t+1)) - 2*np.Gnet(y,z)*(s.v{y}(:,t+1)-s.v{z}(:,t+1))-bp;
                sl.mu_p{y,z}(:,t+1) = zeros(np.h,1);

                s.d_p{y,z}(:,t+1) = s.p_l{y,z}(:,t+1) + s.p_l{z,y}(:,t+1);
                sl.mu_pc{y,z}(:,t+1) = zeros(np.h,1);

                bq = s.q_l{y,z}(:,t+1) - s.q_l{z,y}(:,t+1);
                sl.c_q{y,z}(:,t+1) = 2*np.Gnet(y,z)*(s.theta{y}(:,t+1)-s.theta{z}(:,t+1)) + 2*np.Bnet(y,z)*(s.v{y}(:,t+1)-s.v{z}(:,t+1))-bq;
                sl.mu_q{y,z}(:,t+1) = zeros(np.h,1);

                s.d_q{y,z}(:,t+1) = s.q_l{y,z}(:,t+1) + s.q_l{z,y}(:,t+1);
                sl.mu_qc{y,z}(:,t+1) = zeros(np.h,1);
            end
        end
        


        t = 1;
        %% Iteration
        while 1
            tic
            t
            % 1) Strategy update
           
            
           for i=1:np.n
                
                % primal update of prosumer (quadratic prog.)
                s = loc_opt_qprog_D(np,s,sl,t,i);
                
                % local load imbalance of prosumer i
                sl.b{i}(:,t+1) = np.Pd(i,1:np.h)' - s.p_di{i}(:,t+1) - s.p_st{i}(:,t+1);
                
                % compute error

                res1(i,t+1)= norm(s.u{i}(:,t+1)-s.u{i}(:,t),inf);
            end

            res2(t) = 0; % Extra: compute norm of residual
            dres2(t) = 0;

            % dual variables (reciprocity constraint) Prosumers
            for i=1:np.n

                for j=1:np.n
                    if np.Adj(i,j) == 1

                        % Dual variables (reciprocity constraints) update:
                        
                        % Aux vector
                        sl.c_tr{i,j}(:,t+1) = s.p_tr{i,j}(:,t+1) + s.p_tr{j,i}(:,t+1);
                        
                        % Reflected dual ascent
                        sl.mu_tr{i,j}(:,t+1) = sl.mu_tr{i,j}(:,t) + np.beta_tr(i,j)*(2*sl.c_tr{i,j}(:,t+1) - sl.c_tr{i,j}(:,t));

                        % Extra: compute norm of residual
                        res = sl.c_tr{i,j}(:,t+1);
                        res2(t) = res2(t)+res'*res;

                         
                    end
                end
            end


            % DSO=========================================================================================================================================
             %res2_DSO(t)
            s.sum_sigma(:,t+1) = zeros(np.h,1);
            for y=1:np.b
                
                % Primal update
                a_thet = zeros(np.h,1);
                a_v = zeros(np.h,1);
                a_p =[]; 
                a_q = []; 
                for zz = 1:length(np.B{y})
                    z = np.B{y}(zz);                    
%                    a_thet = a_thet + 4*(np.Bnet(y,z)*sl.mu_p{y,z}(:,t) + np.Gnet(y,z)*sl.mu_q{y,z}(:,t));
%                    a_v = a_v + 4*(np.Bnet(y,z)*sl.mu_q{y,z}(:,t) - np.Gnet(y,z)*sl.mu_p{y,z}(:,t));
                    a_thet = a_thet + 1/2*(sl.mu_p{y,z}(:,t)/np.Bnet(y,z) + sl.mu_q{y,z}(:,t)/np.Gnet(y,z));
                    a_v = a_v + 1/2*(sl.mu_q{y,z}(:,t)/np.Bnet(y,z) - sl.mu_p{y,z}(:,t)/np.Gnet(y,z));
                    a_p = [a_p; sl.mu_pc{y,z}(:,t) - sl.mu_pb{y}(:,t) - sl.mu_p{y,z}(:,t)];
                    a_q = [a_q; sl.mu_qc{y,z}(:,t) - sl.mu_q{y,z}(:,t)];
                end
                a_DSO = [a_thet; a_v; -sl.mu_tg{y}(:,t)-sl.mu_pb{y}(:,t);a_p;a_q];
                
                s = projUy(s,a_DSO,np,y,t); 
                
                % residual
                res1(np.n+y,t+1)= norm(s.u_DSO{y}(:,t+1)-s.u_DSO{y}(:,t),inf);
                
                % aggregation update 
                s.sigma_mg{y}(:,t+1) = np.Pd(np.n+y,1:np.h)';
                sl.by{y}(:,t+1) = zeros(np.h,1);
                for ii = 1:length(np.N_b{y})
                    i = np.N_b{y}(ii);
                    s.sigma_mg{y}(:,t+1) = s.sigma_mg{y}(:,t+1) + s.p_mg{i}(:,t+1);
                    sl.by{y}(:,t+1) = sl.by{y}(:,t+1) + sl.b{i}(:,t+1);
                end
                s.sum_sigma(:,t+1) = s.sum_sigma(:,t+1) + s.sigma_mg{y}(:,t+1);
                
                %auxiliary update
                lam_mg = zeros(2*np.h,1);
                mu_tg = zeros(np.h,1);
                for zz = 1:length(np.B{y})
                    z = np.B{y}(zz);
                    lam_mg = lam_mg + sl.lambda_mg{z}(:,t);
                    mu_tg = mu_tg + sl.mu_tg{z}(:,t);
                end
                sl.w_mg{y}(:,t+1) = sl.w_mg{y}(:,t) + np.delta_mg*(length(np.B{y})*sl.lambda_mg{y}(:,t) - lam_mg );
                sl.w_tg{y}(:,t+1) = sl.w_tg{y}(:,t) + np.delta_tg*(length(np.B{y})*sl.mu_tg{y}(:,t) - mu_tg );
                
                % dual update (global grid constraints)
                sl.c_mg{y}(:,t+1) = kron([1;-1], s.sigma_mg{y}(:,t+1)) -1/np.b*[np.pmg_max*ones(np.h,1);-np.pmg_min*ones(np.h,1)] - sl.w_mg{y}(:,t+1);
                sl.lambda_mg{y}(:,t+1) = max(0, sl.lambda_mg{y}(:,t) + np.gamma_mg(y)*(2*sl.c_mg{y}(:,t+1)-sl.c_mg{y}(:,t)));
                
                sl.c_tg{y}(:,t+1) = s.sigma_mg{y}(:,t+1) - s.p_tg{y}(:,t+1) - sl.w_tg{y}(:,t+1);
                sl.mu_tg{y}(:,t+1) = sl.mu_tg{y}(:,t)+ np.beta_tg(y)*(2*sl.c_tg{y}(:,t+1)-sl.c_tg{y}(:,t));
                
                % dual update (power balance on bus y)
                sl.c_pb{y}(:,t+1) = np.Pd(np.n+y,1:np.h)' + sl.by{y}(:,t+1) - s.p_tg{y}(:,t+1);
                for zz = 1:length(np.B{y})
                    z = np.B{y}(zz);
                    sl.c_pb{y}(:,t+1) = sl.c_pb{y}(:,t+1) - s.p_l{y,z}(:,t+1);
                end
                sl.mu_pb{y}(:,t+1) = sl.mu_pb{y}(:,t) + np.beta_pb(y)*(2*sl.c_pb{y}(:,t+1)- sl.c_pb{y}(:,t));
                
                
            end
            
                % dual update (power flow equation)
            for y=1:np.b
                for zz = 1:length(np.B{y})
                    z = np.B{y}(zz);
                    
                    bp = s.p_l{y,z}(:,t+1) - s.p_l{z,y}(:,t+1);
                    sl.c_p{y,z}(:,t+1) = 2*np.Bnet(y,z)*(s.theta{y}(:,t+1)-s.theta{z}(:,t+1)) - 2*np.Gnet(y,z)*(s.v{y}(:,t+1)-s.v{z}(:,t+1))-bp;
                    sl.mu_p{y,z}(:,t+1) = sl.mu_p{y,z}(:,t) + np.beta_p(y,z)*(2*sl.c_p{y,z}(:,t+1)-sl.c_p{y,z}(:,t));
                    
                    s.d_p{y,z}(:,t+1) = s.p_l{y,z}(:,t+1) + s.p_l{z,y}(:,t+1);
                    sl.mu_pc{y,z}(:,t+1) = sl.mu_pc{y,z}(:,t)+ np.beta_pc*(2*s.d_p{y,z}(:,t+1)-s.d_p{y,z}(:,t));
                    
                    bq = s.q_l{y,z}(:,t+1) - s.q_l{z,y}(:,t+1);
                    sl.c_q{y,z}(:,t+1) = 2*np.Gnet(y,z)*(s.theta{y}(:,t+1)-s.theta{z}(:,t+1)) + 2*np.Bnet(y,z)*(s.v{y}(:,t+1)-s.v{z}(:,t+1))-bq;
                    sl.mu_q{y,z}(:,t+1) = sl.mu_q{y,z}(:,t) + np.beta_q(y,z)*(2*sl.c_q{y,z}(:,t+1)-sl.c_q{y,z}(:,t));
                    
                    s.d_q{y,z}(:,t+1) = s.q_l{y,z}(:,t+1) + s.q_l{z,y}(:,t+1);
                    sl.mu_qc{y,z}(:,t+1) = sl.mu_qc{y,z}(:,t)+ np.beta_qc*(2*s.d_q{y,z}(:,t+1)-s.d_q{y,z}(:,t));
                    
                end
            end
           
                                    

            % Extra: stopping criterion
            s.res(t) = sqrt(res2(t));
            
            %s.dres(t) = sqrt(dres2(t));
            
            
%            s.res_dso(t) = sqrt(res_dso2(t));
            
 %           s.dres_dso(t) = sqrt(dres_dso2(t));
            
            s.res1(t) = norm(res1(:,t+1),inf);%sqrt(res1(t+1));
            [s.res1(t); s.res(t)]
            %[s.res(t);s.dres(t);s.dres_dso(t);s.res_dso(t);s.res1(t)]
            %[s.res(t);s.dres(t)]
  %          if s.res(t) < np.er_max && s.dres(t) < np.er_max && s.res_dso(t) < np.er_max && s.dres_dso(t) < np.er_max && s.res1(t)<np.er_max %&& all(s.ph_mg(:,t) >= np.pmg_min) && all(s.ph_mg(:,t) <= np.pmg_max*ones(np.h,1)) 
           % if  s.res(t) < np.er_max && s.dres(t) < np.er_max
            if s.res1(t) < np.er_max
                disp('solution obtained')
                break
            end
            
            if t > np.t_max
                disp('max iteration reached')
                break
            end

            t= t+1;
            toc
        end
        
end