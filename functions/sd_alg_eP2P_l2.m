function [s,sl,np] = sd_alg_eP2P_l2(np)
%% Distributed algorithm
% extended P2P market
% W. Ananduta
% 12/11/2020

        %% INITIALIZATION
        % Assigning parameters of the algorithm
        %np = alg_param(np);
        np = alg_param_37b(np);
        np.t_max = 1000;
        np.er_max = 1e-6; 
        flag_eps2 = 0;
        % % Generate matrices S_i^mg and S_ij^tr
        %[np.Sdg,np.Sst,np.Smg,np.Str] = gen_Smat(np);

        % Generate matrices for cost function and constraints
        np = build_mat_exP2P(np);
        
        %np = build_mat_exP2P_tr(np);
        np = build_mat_exDSO(np);
        
        %% Initialization of the variables
        %s.ph_mg(:,1) = zeros(np.h,1);
        t=1;
        s.sigma_mg = zeros(np.h,1);
        for i=1:np.n
            % decision variables
            %s = initialize_u_quadp(np,s,i);
            s.u{i}(:,1) = np.init*ones(size(np.A_ineq{i},2),1);
            s.p_di{i}(:,1) = np.Sdi{i}*s.u{i}(:,1);
            s.p_st{i}(:,1) = np.Sst{i}*s.u{i}(:,1);
            s.p_mg{i}(:,1) = np.Smg{i}*s.u{i}(:,1);
            for jj = 1:length(np.N{i})
                j = np.N{i}(jj);
                s.p_tr{i,j} = np.Str{i,j}*s.u{i}(:,1);
            end
            sl.b{i}(:,1) = np.Pd(i,1:np.h)' - s.p_di{i}(:,1) - s.p_st{i}(:,1);

            s.sigma_mg(:,1) = s.sigma_mg(:,1) + s.p_mg{i}(:,1);
        end
        s.lambda_mg= zeros(2*np.h,1);

        for i=1:np.n
            for jj=1:length(np.N{i})
                j = np.N{i}(jj);
                sl.c_tr{i,j}(:,1) = s.p_tr{i,j}(:,1) + s.p_tr{j,i}(:,1);
                
                s.mu_tr{i,j}(:,1) = np.init*ones(np.h,1);
            end
        end

        % Initialization of DSO's decision variables         
        s = projDSO(s,0,np,0);                                     % 
        
        s.sum_p_pd = zeros(np.h,1);
        for y = 1:np.b


            sl.c_pb{y}(:,1) = np.Pd(np.n+y,1:np.h)' - s.p_tg{y}(:,1);
            for ii = 1:length(np.N_b{y})
                i = np.N_b{y}(ii);
                sl.c_pb{y}(:,1) = sl.c_pb{y}(:,1) + sl.b{i}(:,1);
            end

            for zz = 1:length(np.B{y})
                z = np.B{y}(zz);
                sl.c_pb{y}(:,1) = sl.c_pb{y}(:,1) - s.p_l{y,z}(:,1);
            end

            sl.mu_pb{y}(:,1) = np.init*ones(np.h,1);

            s.sum_p_pd = s.sum_p_pd + np.Pd(np.n+y,1:np.h)';

        end
        s.sigma_mg(:,1) = s.sigma_mg(:,1) + s.sum_p_pd;

        sl.c_tg(:,1) = s.sigma_mg(:,1);
        for yy = 1:length(np.B_mg)
            y = np.B_mg(yy);
            sl.c_tg(:,1) = sl.c_tg(:,1) - s.p_tg{y}(:,1);                
        end
        sl.mu_tg(:,1) = np.init*ones(np.h,1);



        %% Iteration
        while 1
            tic
            t
            % 1) Strategy update
            %s = loc_opt_c_qprog(np,s,t);

            %s.ph_mg(:,t+1) = 0;
            s.sigma_mg(:,t+1) = s.sum_p_pd;
            
            
           for i=1:np.n
                
                % primal update of prosumer (quadratic prog.)
                s = loc_opt_qprog(np,s,sl,t,i);
                
                % local load imbalance of prosumer i
                sl.b{i}(:,t+1) = np.Pd(i,1:np.h)' - s.p_di{i}(:,t+1) - s.p_st{i}(:,t+1);

                % forward p_mg to DSO
                s.sigma_mg(:,t+1) = s.sigma_mg(:,t+1) + s.p_mg{i}(:,t+1);
                
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
                        s.mu_tr{i,j}(:,t+1) = s.mu_tr{i,j}(:,t) + np.beta_tr(i,j)*(2*sl.c_tr{i,j}(:,t+1) - sl.c_tr{i,j}(:,t));

                        % Extra: compute norm of residual
                        res = sl.c_tr{i,j}(:,t+1);
                        res2(t) = norm([res2(t);res],inf);

                        % EXTRA: Compute dual residual                        
                        dres2(t) = dres2(t) + (s.mu_tr{i,j}(:,t+1) - s.mu_tr{i,j}(:,t))'*(s.mu_tr{i,j}(:,t+1) - s.mu_tr{i,j}(:,t)); 
                    end
                end
            end


            % DSO=========================================================================================================================================
            % Primal update 
            
            for y=1:np.b
                a_DSO{y} = [zeros(2*np.h,1); -sl.mu_tg(:,t)-sl.mu_pb{y}(:,t);kron(ones(length(np.B{y}),1),-sl.mu_pb{y}(:,t));zeros(np.h*length(np.B{y}),1)]; 
            end
           
            s = projDSO(s,a_DSO,np,t);                                         
                                   
            % Dual variable (grid constraints)
            sl.b_DSO(:,t+1) = kron([1;-1],2*s.sigma_mg(:,t+1)-s.sigma_mg(:,t))-[np.pmg_max*ones(np.h,1);-np.pmg_min*ones(np.h,1)];
            s.lambda_mg(:,t+1) = max(0, s.lambda_mg(:,t) + np.gamma_mg*sl.b_DSO(:,t+1));        

            % Dual variable (local power balance of busses)
            res_dso2(t)=0;
            dres_dso2(t)=0;
            for y=1:np.b
                sl.c_pb{y}(:,t+1) = np.Pd(np.n+y,1:np.h)' - s.p_tg{y}(:,t+1);
                for ii = 1:length(np.N_b{y})
                    i = np.N_b{y}(ii);
                    sl.c_pb{y}(:,t+1) = sl.c_pb{y}(:,t+1) + sl.b{i}(:,t+1);
                end

                for zz = 1:length(np.B{y})
                    z = np.B{y}(zz);
                    sl.c_pb{y}(:,t+1) = sl.c_pb{y}(:,t+1) - s.p_l{y,z}(:,t+1);
                end
                sl.mu_pb{y}(:,t+1) = sl.mu_pb{y}(:,t) + np.beta_pb(y)*(2*sl.c_pb{y}(:,t+1) - sl.c_pb{y}(:,t)) ;
                
                % error
                res_dso2(t) = res_dso2(t) + sl.c_pb{y}(:,t+1)'*sl.c_pb{y}(:,t+1);
                dres_dso2(t) = dres_dso2(t) + (sl.mu_pb{y}(:,t+1)-sl.mu_pb{y}(:,t))'*(sl.mu_pb{y}(:,t+1)-sl.mu_pb{y}(:,t));
                
                res1(np.n+y,t+1)= norm(s.u_DSO{y}(:,t+1)-s.u_DSO{y}(:,t),inf);
            end

            % Dual variable (trading with the main grid)
            s.sigma_tg(:,t+1) = zeros(np.h,1) ;
            for yy = 1:length(np.B_mg)
                y = np.B_mg(yy);
                s.sigma_tg(:,t+1) = s.sigma_tg(:,t+1) + s.p_tg{y}(:,t+1);                
            end
            sl.c_tg(:,t+1) = s.sigma_mg(:,t+1) - s.sigma_tg(:,t+1);
            sl.mu_tg(:,t+1) = sl.mu_tg(:,t) + np.beta_tg*(2*sl.c_tg(:,t+1)-sl.c_tg(:,t));

            % Extra: stopping criterion
            s.res(t) = (res2(t));
            
            s.dres(t) = sqrt(dres2(t));
            
            
            s.res_dso(t) = sqrt(res_dso2(t));
            
            s.dres_dso(t) = sqrt(dres_dso2(t));
            
            s.res1(t) = norm(res1(:,t+1),inf);
            [s.res(t);s.dres(t);s.dres_dso(t);s.res_dso(t);s.res1(t)]
            
            %if s.res(t) < np.er_max && s.dres(t) < np.er_max && s.res_dso(t) < np.er_max && s.dres_dso(t) < np.er_max && s.res1(t)<np.er_max %&& all(s.ph_mg(:,t) >= np.pmg_min) && all(s.ph_mg(:,t) <= np.pmg_max*ones(np.h,1)) 
            %if  s.res(t) < np.er_max && s.dres(t) < np.er_max
            if s.res1(t) < 0.01 && s.res(t) < 0.01 && flag_eps2==0
                [s,o] = com_cost(s,np);
                o.err = s.res1;
                save(['sim37b_eps-1_5_B_',date,'_'],'o' )
                flag_eps2 = 1;
                
            end
            if s.res1(t) < np.er_max && s.res(t) < np.er_max
                break
            end
            

            t= t+1;
            toc
        end
end