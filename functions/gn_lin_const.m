function p = gn_lin_const(p)
% linear constraints of gas network
% W. Ananduta
% 11/08/2021

gn = p.gn;
for i = 1:p.gn.n
           
    % Bounds on the pressures [(17), GS_for_ED_games_on_IEGS.pdf]
    % upper bound
    E4 = [1 zeros(1,gn.noN(i))];
    E4t = kron(eye(p.h),E4);
    e4t = gn.pc.psi_max(i)*ones(p.h,1);
    
    % lower bound
    E5 = -E4;
    E5t = kron(eye(p.h),E5);
    e5t = -gn.pc.psi_min(i)*ones(p.h,1);
    
    gn.Aineq{i} = [E4t;E5t];
    gn.bineq{i} = [e4t;e5t];
    
    % Bounds on the flow (16)
    for jj = 1:gn.noN(i)
        j = gn.N{i}(jj);
        
        % upper bound
        E6 = zeros(1,1+gn.noN(i));
        E6(1,1+jj) = 1;
        E6t = kron(eye(p.h),E6);
        e6t = gn.pc.phi_max(i)*ones(p.h,1);
        
        % lower bound
        E7 = -E6;
        E7t = kron(eye(p.h),E7);
        e7t = gn.pc.phi_max(i)*ones(p.h,1);
        
        gn.Aineq{i} = [gn.Aineq{i};E6t;E7t];
        gn.bineq{i} = [gn.bineq{i};e6t;e7t];
    end
    
    
    
end

p.gn = gn;
end