function [s,p,o] = iterate_GNEPc_pen(s,p)
% Solve convexified GNEP of ED in IEGDS
% Via PPP method
% W. Ananduta
% 04/10/2021

    p.t_max = 1e6;
    p.er_max =  1e-3; 
        k=1;
        %% ITERATION
        for k = 1:p.t_max
            %tic
            
            
            
            s.sigma_mg(:,k+1) = p.en.sumPd(1:p.h);
            s.sigma_gu(:,k+1) = p.gn.sumGd(1:p.h);
            
            
            % primal update of prosumers 
            for i=1:p.n
                
                % (quadratic prog.)
                s = loc_opt_qprog_pen(p,s,k,i); 
                
                s.sigma_mg(:,k+1) = s.sigma_mg(:,k+1) + s.p_mg{i}(:,k+1);
                s.sigma_gu(:,k+1) = s.sigma_gu(:,k+1) + s.d_gu{i}(:,k+1);
                
                
                % compute error
                res1(i,k+1)= norm(s.u{i}(:,k+1)-s.u{i}(:,k),inf);
            end
            
            
            res2(k) = 0; % Extra: compute norm of residual
            res3(k) = 0; % Extra: compute norm of residual
            res4(k) = 0;
            
            % dual variable updates
            for i=1:p.n
                
                
                for jj=1:length(p.en.N{i})
                j = p.en.N{i}(jj);
                
                % dual variable for power flow constraints
                sl.c_pl{i,j}(:,k+1) = p.m.PF{i,j}{1}*s.u{i}(:,k+1) + p.m.PF{i,j}{2}*s.u{j}(:,k+1);
                
                s.mu_pl{i,j}(:,k+1) = s.mu_pl{i,j}(:,k) + p.beta_pl(i)*(2*sl.c_pl{i,j}(:,k+1) - sl.c_pl{i,j}(:,k));
                
                resc = sl.c_pl{i,j}(:,k+1);
                
                res4(k) = norm([res4(k);resc],inf);
                
                end
                
                
                for jj = 1:length(p.gn.N{i})
                    j = p.gn.N{i}(jj);
                    
                    % dual variable for reciprocity of \phi
                    sl.c_phi{i,j}(:,k+1) = s.phi{i,j}(:,k+1) + s.phi{j,i}(:,k+1);

                    s.mu_phi{i,j}(:,k+1) = s.mu_phi{i,j}(:,k) + p.beta_phi(i)*(2*sl.c_phi{i,j}(:,k+1) - sl.c_phi{i,j}(:,k));

                    % dual variable for gas-flow equation
                    sl.c_gf{i,j}(:,k+1) = p.m.Hc{i,j}{1}*s.u{i}(:,k+1) + p.m.Hc{i,j}{2}*s.u{j}(:,k+1);

                    s.mu_gf{i,j}(:,k+1) = s.mu_gf{i,j}(:,k) + p.beta_gf(i)*(2*sl.c_gf{i,j}(:,k+1) - sl.c_gf{i,j}(:,k));
                    
                    % dual variable for psi coupling const
                    sl.c_psi{i,j}(:,k+1) = p.m.G{i,j}{1}*s.u{i}(:,k+1) + p.m.G{i,j}{2}*s.u{j}(:,k+1) - p.m.g{i,j};

                    s.lambda_psi{i,j}(:,k+1) = max(0, s.lambda_psi{i,j}(:,k) + p.beta_psi(i)*(2*sl.c_psi{i,j}(:,k+1) - sl.c_psi{i,j}(:,k)) );
                    
                    % Extra: compute norm of residual
                        res = sl.c_gf{i,j}(:,k+1);
                        resb = sl.c_phi{i,j}(:,k+1);
                        %res = sign(sl.c_psi{i,j}(:,k+1));
                        res2(k) = norm([res2(k);res],inf);
                        %res2(k) = max([res2(k);res]);
                        res3(k) = norm([res3(k);resb],inf);
                end
                
                                
            end
            
            % dual variable (grid constraint) update
            s.lambda_mg(:,k+1) = max(0, s.lambda_mg(:,k) + p.delta*(2*[s.sigma_mg(:,k+1);-s.sigma_mg(:,k+1)]-[s.sigma_mg(:,k);-s.sigma_mg(:,k)] ...
                                        - [p.en.pmg_max*ones(p.h,1);-p.en.pmg_min*ones(p.h,1)]));

            % dual variable (gas coupling constraint) update
            s.lambda_gu(:,k+1) = max(0, s.lambda_gu(:,k) + p.delta_gu*(2*[s.sigma_gu(:,k+1);-s.sigma_gu(:,k+1)]-[s.sigma_gu(:,k);-s.sigma_gu(:,k)] ...
                                        - [p.gn.dg_max*ones(p.h,1);-p.gn.dg_min*ones(p.h,1)]));
            
                                    
            % Extra: stopping criterion
            s.res(k) = (res2(k));
            s.resb(k) = (res3(k));
            s.resc(k) = res4(k);
            %s.dres(t) = sqrt(dres2(t));
            
            
            
            s.res1(k) = norm(res1(:,k+1),inf);
            
            s.residual(:,k) = [s.res(k);s.resb(k);s.resc(k);s.res1(k)];
            if mod(k,100) == 0
                k
                 s.residual(:,k)
%                   [s.res(k);s.res1(k)]
                 r.residual = s.residual;
                 save(['iteration_ppp'],'r')
            end
            
          
            
           
           if  norm(s.residual(:,k),inf) < p.er_max 
               % if  s.res1(k) < p.er_max
               break
           end
 %           

            %k= k+1;
            %toc
        end
        o.n_iter = k;
        o.residual = s.residual(:,end);
        for i = 1:p.n
            o.u{i} = s.u{i}(:,end);
            
            o.p_di{i} = p.m.Sdg{i}*o.u{i};
            o.p_ch{i} = p.m.Sch{i}*o.u{i};
            o.p_dh{i} = p.m.Sdh{i}*o.u{i};
            o.p_mg{i} = p.m.Smg{i}*o.u{i};
            o.d_gu{i} = p.m.Sgu{i}*o.u{i};
            o.th{i} = p.m.Sth{i}*o.u{i}/p.scaling;
            o.v{i} = p.m.Sv{i}*o.u{i}/p.scaling;
            for jj = 1:length(p.en.N{i})
                j = p.en.N{i}(jj);
                o.p_l{i,j} = p.m.Spl{i,j}*o.u{i};
            end

        % y
            o.psi{i} = p.m.Spsi{i}*o.u{i};
            
            for jj = 1:length(p.gn.N{i})
                j = p.gn.N{i}(jj);
                o.phi{i,j} = p.m.Sphi{i,j}*o.u{i};
                o.pen{i,j} = p.m.Spen{i,j}*o.u{i};
            end
        end
        
        p.u_init = o.u;
        
end