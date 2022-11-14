function Ji_best = loc_opt_NI_c(o,p,i)
% Local optimization for MI-NI
% W. Ananduta
% 10/11/2021

    c = zeros(p.nu(i),1);
    c(1) = p.en.c_dg(i);
    c(2) = p.en.c_st(i);
    

    for jj = 1:length(p.en.N{i})
        j = p.en.N{i}(jj);
        c(4+jj,1) = p.en.c_tr(i,j);
        
    end
    
    c1 = [];
    
    
    
    for h = 1:p.h
        c_1 = c;
        c_1(3) = p.en.d_mg*(o.sigma_mg(h)-o.p_mg{i}(h));
        c_1(4) = p.gn.d_gu*(o.sigma_gu(h)-o.d_gu{i}(h));
        c1 = [c1;c_1];
    end 
       
    f = c1;
    
    xi = sdpvar(p.nx(i),1);
    yi = sdpvar(p.ny(i),1);
    zi = sdpvar(p.nz(i),1);
    
    ui = [xi;yi;zi];
    
    J = ui'*p.m.Qh{i}*ui + f'*ui;
    
    


    A = p.m.Aineq{i};
    b = p.m.bineq{i};
    Aeq = p.m.Aeq{i};
    beq = p.m.beq{i};
    cons = [];
    cons = [cons, A*ui <= b];
    cons = [cons, Aeq*ui == beq];
    

    
    for jj=1:length(p.en.N{i})
        j = p.en.N{i}(jj);

        cons = [cons,p.m.Str{i,j}*ui + o.p_tr{j,i} == 0];

    end


    for jj = 1:length(p.gn.N{i})
        j = p.gn.N{i}(jj);

        % dual variable for reciprocity of \phi
        cons = [cons,p.m.Sphi{i,j}*ui + o.phi{j,i}==0];

        uj = o.u{j};
        uj(p.nx(j)+1:p.nx(j)+p.h,1) = o.psi{j};
        
        % dual variable for gas-flow equation
        cons = [cons,p.m.Hc{i,j}{1}*ui + p.m.Hc{i,j}{2}*uj ==0];

        
        % dual variable for psi coupling const
        cons = [cons,p.m.G{i,j}{1}*ui + p.m.G{i,j}{2}*uj - p.m.g{i,j}<=0];

    end
    
    options = sdpsettings('verbose',1,'solver','osqp','gurobi.OptimalityTol',1e-9,'osqp.eps_abs',1e-10);
    sol = optimize(cons,J,options);
   
    Ji_best = value(J);
    
    


end