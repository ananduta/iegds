function [s,p] = gs_iegds(p)
% Gauss-Southwell for ED in IEGDS
% W. Ananduta
% 11/08/2021


% Define local constraints of prosumers
np = p.tn;
np = build_mat_exP2P(np);
p.tn = np;
clearvars('np');


% Define local constraints of network operator
p = gn_lin_const(p); % linear constraints (16)-(17)
p = pwa_gasflow(p);  % PWA approximation of gas-flow equations (23)-(24)



% Initialization
k=1;

s.sigma_mg(:,1) = zeros(p.h,1);
s.sigma_gu(:,1) = zeros(p.h,1);
for i = 1:p.n
    s.u{i}(:,1)  = init_prosumer(p,i);
    s.sigma_mg(:,1) = s.sigma_mg(:,1) + p.tn.Smg{i}*s.u{i}(:,1);
    s.sigma_gu(:,1) = s.sigma_gu(:,1) + p.tn.Sgu{i}*s.u{i}(:,1);
end
s  = init_no(s,p);

% iteration
kmax = 1e4;
while 1
    for i = 1:p.n
        
        % prosumer i is chosen

        % compute local cost
        Ji_now = local_cost(s,p,i,k);

        % compute a best response and cost value
        
        [u_i,Ji_new] = br_iegds(s,p,i,k);
        

        
        % update decision of the chosen prosumer
        if Ji_now - Ji_new > p.epsi
            s.u{i}(:,k+1) = u_i;
        else
            s.u{i}(:,k+1) = s.u{i}(:,k);
        end
        
        % update decisions of all other players (including NO)
        for j = 1:p.n
            if j ~= i
                s.u{j}(:,k+1) = s.u{j}(:,k);
            end
            s.u_no{j}(:,k+1) = s.u_no{j}(:,k);
        end
        
        % update aggregative term
        s.sigma_mg(:,k+1) = s.sigma_mg(:,k) - p.tn.Smg{i}*s.u{i}(:,k) + p.tn.Smg{i}*s.u{i}(:,k+1);
        s.sigma_gu(:,k+1) = s.sigma_gu(:,k) - p.tn.Sgu{i}*s.u{i}(:,k) + p.tn.Sgu{i}*s.u{i}(:,k+1);
        k = k+1;
    end
    
    % compute local cost
    Jno_now = local_cost_NO(s,p,k);

    % compute a best response and cost value
    [u_no,Jno_new] = br_iegds_NO(s,p,k);
    
    % update decision of NO
    if Jno_now - Jno_new > p.epsi
        for i = 1:p.gn.n
            s.u_no{i}(:,k+1) = u_no{i};
        end
    else
        for i = 1:p.gn.n
            s.u_no{i}(:,k+1) = s.u_no{i}(:,k);
        end
        
    end
    
    % update decision of all prosumers
    for j = 1:p.n
        s.u{j}(:,k+1) = s.u{j}(:,k);
    end
    s.sigma_mg(:,k+1) = s.sigma_mg(:,k);
    s.sigma_gu(:,k+1) = s.sigma_gu(:,k);
    
    k = k+1;
    
end

end