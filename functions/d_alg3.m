function [s,np] = d_alg3(np)

%% Distributed algorithm
            %% INITIALIZATION

            % % Generate matrices S_i^mg and S_ij^tr
            %[np.Sdg,np.Sst,np.Smg,np.Str] = gen_Smat(np);

            % Generate matrices for cost function and constraints
            np = build_mat(np);

            % Initialization of the variables
            s.ph_mg(:,1) = zeros(np.h,1);
            t=1;


            for i=1:np.n

                % decision variables
                s = initialize_u(np,s,i);
                % Collect p_i^mg(1)
                s.ph_mg(:,1) = s.ph_mg(:,1) + s.pmg{i}(:,1);

                % auxiliary and dual variables
                s.w{i} = zeros(2*np.h,1);
                for j=1:np.n
                    s.mu{i,j} = zeros(np.h,1);
                end
                s.lambda{i}= zeros(2*np.h,1);


            end

            % Assigning parameters of the algorithm
            np = alg_param(np);
            np.t_max = 1000;
            np.er_max = 0.01;

            %% Iteration
            while 1
                tic
                t
                % 1) Strategy update
                s = loc_opt_c(np,s,t);
%                 s.ph_mg(:,t+1) = 0;
%                 for i=1:np.n
%                     s = loc_opt(np,s,i,t);
% 
%                     % Collect p_i^mg(t)
%                     s.ph_mg(:,t+1) = s.ph_mg(:,t+1) + s.pmg{i}(:,t+1);
%                 end

                res2(t) = 0; % Extra: compute norm of residual
                dres2(t) = 0;
                
                for i=1:np.n

                    s.w{i}(:,t+1) = s.w{i}(:,t); %4) Auxiliary variable update

                    for j=1:np.n
                        if np.Adj(i,j) == 1

                            % 3) Dual variables (reciprocity constraints) update:
                            s.mu{i,j}(:,t+1) = s.mu{i,j}(:,t) - np.beta(i,j)*(np.Str{i,j}*s.u{i}(:,t) + np.Str{j,i}*s.u{j}(:,t)...
                                                - 2*np.Str{i,j}*s.u{i}(:,t+1) - 2*np.Str{j,i}*s.u{j}(:,t+1));

                            % Extra: compute norm of residual
                            res = np.Str{i,j}*s.u{i}(:,t+1)+np.Str{j,i}*s.u{j}(:,t+1);
                            res2(t) = res2(t)+res'*res;
                            
                            % EXTRA: Compute dual residual                        
                            dres2(t) = dres2(t) + (s.mu{i,j}(:,t+1) - s.mu{i,j}(:,t))'*(s.mu{i,j}(:,t+1) - s.mu{i,j}(:,t)); 
                            % 4) Auxiliary variable update:
                            s.w{i}(:,t+1) = s.w{i}(:,t+1)+ np.gamma(i)*(s.lambda{i}(:,t)-s.lambda{j}(:,t));

                        end
                    end

                    %5) Dual variable (grid constraints) update:        
                    s.lambda{i}(:,t+1) = max(0, s.lambda{i}(:,t)+ np.delta(i)*(2*[np.Smg{i}*s.u{i}(:,t+1);-np.Smg{i}*s.u{i}(:,t+1)]...
                                                -[np.Smg{i}*s.u{i}(:,t);-np.Smg{i}*s.u{i}(:,t)]-[np.pmg_max/np.n*ones(np.h,1);-np.pmg_min/np.n*ones(np.h,1)]...
                                                -2*s.w{i}(:,t+1) + s.w{i}(:,t)));
                    
                    
                end
                % Extra: stopping criterion
                s.res(t) = sqrt(res2(t));
                s.res(t)
                s.dres(t) = sqrt(dres2(t));
                s.dres(t)
                
                if s.res(t) < np.er_max && s.dres(t) < np.er_max %&& all(s.ph_mg(:,t) >= np.pmg_min) && all(s.ph_mg(:,t) <= np.pmg_max*ones(np.h,1)) 
                    break
                end


                t= t+1;
                toc
            end
end