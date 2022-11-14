function [s,np] = d_alg_qp(np)
%% Distributed algorithm
% extended P2P market
% W. Ananduta
% 17/09/2020

            %% INITIALIZATION
            % Assigning parameters of the algorithm
            np = alg_param(np);
            np.t_max = 1000;
            np.er_max = 0.01;
            
            % % Generate matrices S_i^mg and S_ij^tr
            %[np.Sdg,np.Sst,np.Smg,np.Str] = gen_Smat(np);

            % Generate matrices for cost function and constraints
            np = build_mat(np);                                             % MUST BE EDITED, ADJUSTING CURRENT STRUCTURE OF u AND INCLUDE DSO

            %% Initialization of the variables
            %s.ph_mg(:,1) = zeros(np.h,1);
            t=1;
            s.sigma = zeros(np.h,1);
            for i=1:np.n
                % decision variables
                s = initialize_u_quadp(np,s,i);
                
                s.b{i}(:,1) = np.pd(:,i) - s.p_di{i}(:,1) - s.p_st{i}(:,1);
                
                % auxiliary and dual variables
                s.w{i} = zeros(2*np.h,1);
                for j=1:np.n
                    s.mu{i,j} = zeros(np.h,1);
                end
                

                s.sigma(:,1) = s.sigma(:,1) + s.pmg{i}(:,1);
            end
            s.lambda= zeros(2*np.h,1);
            
            for i=1:np.n
                for jj=1:length(np.N_{i})
                    s.c{i,j}(:,1) = s.p_tr{i,j}(:,1) + s.p_tr{j,i}(:,1);
                end
            end
            
            % Initialization of DSO's decision variables
            s = initialize_uDSO(np,s);                                      % MUST BE CREATED
            
            s.sum_p_pd = zeros(np.h,1);
            for y = 1:np.b
                 
                
                s.d{y}(:,1) = np.p_pd(:,y) - s.p_tg{y}(:,1);
                for ii = 1:length(np.N_b{y})
                    i = np.N_b{y}(ii);
                    s.d{y}(:,1) = s.d{y}(:,1) + s.b{i}(:,1);
                end
                
                for zz = 1:length(np.B_{y})
                    z = np.B_{y}(zz);
                    s.d{y}(:,1) = s.d{y}(:,1) + s.p_l{y,z}(:,1);
                end
                
                s.nu{y}(:,1) = zeros(np.h,1);
                
                s.sum_p_pd = s.sum_p_pd + np.p_pd(:,y);
                
            end
            s.sigma(:,1) = s.sigma(:,1) + s.sum_p_pd;
            
            s.ph_tg(:,1) = s.sigma(:,1);
            for yy = 1:length(np.B_mg)
                y = np.B_mg(yy);
                s.ph_tg(:,1) = s.ph_tg(:,1) - s.p_tg{y}(:,1);                
            end
            s.rho(:,1) = zeros(np.h,1);

            

            %% Iteration
            while 1
                tic
                t
                % 1) Strategy update
                %s = loc_opt_c_qprog(np,s,t);
                
                %s.ph_mg(:,t+1) = 0;
                s.sigma(:,t+1) = s.sum_p_pd;
                
                for i=1:np.n
                    mu = [];
                    for jj = 1:length(np.N_{i})
                        j = np.N_{i}(jj);
                        mu = [mu;s.mu{i,j}(:,t)];
                    end
                    s.a{i}(:,t) = [-s.nu{np.bus_of(i)}(:,t); -s.nu{np.bus_of(i)}(:,t); (kron([1;-1],eye(np.h))'*s.lambda(:,k)+s.rho(:,k); mu]; 
                    s = loc_opt_qprog(np,s,t,i);
                    
                    s.b{i}(:,t+1) = np.pd(:,i) - s.p_di{i}(:,t+1) - s.p_st{i}(:,t+1);
                    
                    
                end
                
                % Collect p_i^mg(t)
                s.sigma(:,t+1) = s.sigma(:,t+1) + s.pmg{i}(:,t+1);
                    
                    
                end

                %res2(t) = 0; % Extra: compute norm of residual
                %dres2(t) = 0;
                
                % dual variables (reciprocity constraint) Prosumers
                for i=1:np.n

                %    s.w{i}(:,t+1) = s.w{i}(:,t); %4) Auxiliary variable update

                    for j=1:np.n
                        if np.Adj(i,j) == 1

                            % 3) Dual variables (reciprocity constraints) update:
                            s.c{i,j}(:,t+1) = s.p_tr{i,j}(:,t+1) + s.p_tr{j,i}(:,t+1);
                            s.mu{i,j}(:,t+1) = s.mu{i,j}(:,t) + np.beta{i,j}*(2*s.c{i,j}(:,t+1) - s.c{i,j}(:,t));

                            % Extra: compute norm of residual
                            %res = np.Str{i,j}*s.u{i}(:,t+1)+np.Str{j,i}*s.u{j}(:,t+1);
                            %res2(t) = res2(t)+res'*res;
                            
                            % EXTRA: Compute dual residual                        
                            %dres2(t) = dres2(t) + (s.mu{i,j}(:,t+1) - s.mu{i,j}(:,t))'*(s.mu{i,j}(:,t+1) - s.mu{i,j}(:,t)); 
                        end
                    end
                
                    
                % DSO
                % Agregate strategy update
                
                % Dual variable (grid constraints)
                s.lambda(:,t+1) = max(0, s.lambda(:,t) + np.gamma*(kron([1;-1],2*s.sigma(:,t+1)-s.sigma(:,t))-[np.pmg_max*ones(np.h,1);-np.pmg_min*ones(np.h,1)]));
                
                % Update of the physical variables 
                s.a{np.n+1}(:,t+1) = [zeros(2*np.h,1); -s.rho(:,t)];
                for y=1:np.b
                    s.a{np.n+1}(:,t+1) = [s.a{np.n+1}(:,t+1);repmat([s.nu{y}(:,k); zeros(np.h,1)],length(np.B_{y}),1)]; 
                end
                s = proj2U_DSO(s,np);                                         % EITHER BY ALTERNATING PROJECTION METHOD OR DIRECTLY SOLVE A CONVEX PROBLEM
                
                
                % Dual variable (local power balance of busses)
                for y=1:np.b
                    s.d{y}(:,t+1) = np.p_pd(:,y) - s.p_tg{y}(:,t+1);
                    for ii = 1:length(np.N_b{y})
                        i = np.N_b{y}(ii);
                        s.d{y}(:,t+1) = s.d{y}(:,t+1) + s.b{i}(:,t+1);
                    end

                    for zz = 1:length(np.B_{y})
                        z = np.B_{y}(zz);
                        s.d{y}(:,t+1) = s.d{y}(:,t+1) + s.p_l{y,z}(:,t+1);
                    end
                    s.nu{y}(:,t+1) = s.nu{y}(:,t+1) + np.delta(y)*(2*s.d{y}(:,t+1) - s.d{y}(:,t)) 
                
                end
                
                % Dual variable (trading with the main grid)
                s.ph_tg(:,t+1) = s.sigma(:,t+1);
                for yy = 1:length(np.B_mg)
                    y = np.B_mg(yy);
                    s.ph_tg(:,t+1) = s.ph_tg(:,t+1) - s.p_tg{y}(:,t+1);                
                end
                s.rho(:,t+1) = s.rho(:,t) + np.eta*(2*s.ph_tg(:,t+1)-s.ph_tg(:,t));
                
                % Extra: stopping criterion
                %s.res(t) = sqrt(res2(t));
                %s.res(t)
                %s.dres(t) = sqrt(dres2(t));
                %s.dres(t)
                
                %if s.res(t) < np.er_max && s.dres(t) < np.er_max %&& all(s.ph_mg(:,t) >= np.pmg_min) && all(s.ph_mg(:,t) <= np.pmg_max*ones(np.h,1)) 
                %    break
                %end


                t= t+1;
                toc
            end
end