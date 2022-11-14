function np = build_mat_exDSO(np)
% Construct matrices for optimization of DSO
% Extended P2P Market
% 25/09/2020
% W. Ananduta

np.le_u = zeros(np.b,1);
for y = 1:np.b
    np.le_u(y) = 3+2*length(np.B{y});
    
    % lower and upper bounds of voltage angles (eq. 18)
    a1= [1 zeros(1,2+2*length(np.B{y}))];
    A1 = kron([a1;-a1],eye(np.h));
    b1 = kron([np.theta_max(y);-np.theta_min(y)],ones(np.h,1));
    
    % lower and upper bounds of voltage magnitudes (eq. 19)
    a2= [0 1 zeros(1,1+2*length(np.B{y}))];
    A2 = kron([a2;-a2],eye(np.h));
    b2 = kron([np.v_max(y);-np.v_min(y)],ones(np.h,1));
    
    np.Aineq_dso{y} = [A1;A2];
    np.bineq_dso{y} = [b1;b2];
    
end

% coupling constraints (15) and (16)
c=0;
for y=1:np.b
    
    % start and end index of bus y
    y_st = sum(np.le_u(1:y-1))+1;
    y_en = sum(np.le_u(1:y));
        
    for zz = 1:length(np.B{y})
        z=np.B{y}(zz);
                
        % start and end index of bus z
        z_st = sum(np.le_u(1:z-1))+1;
        z_en = sum(np.le_u(1:z));
        
        % (15) eq. coupling const for real power line (y,z)
        c=c+1;
        E = zeros(1,sum(np.le_u));                                    %initialize matrix for (15) coup. const. of real power line
        
        y_au = zeros(1,length(np.B{y}));
        %y_au(zz) = -1/1000;
        y_au(zz) = -1; % for 37-bus system
        
        E(y_st:y_en) = [np.Bnet(y,z) -np.Gnet(y,z) 0 y_au zeros(1,length(np.B{y}))]; %encode components of bus y for (15) coup. const. of real power line          
        
        E(z_st:z_en)= [-np.Bnet(y,z) np.Gnet(y,z) 0 zeros(1,2*length(np.B{z}))]; %encode components of bus z for (15) coup. const. of real power line          
        
        % full matrix over time horizon for (15) coup. const. of real power line
        Ef{c} = kron(E,eye(np.h));                                          
        
        
        % (16) eq. coupling const for reactive power line (y,z)
        c=c+1;
        E = zeros(1,sum(np.le_u));                                    %initialize matrix for (16) coup. const. of reactive power line
                
        E(y_st:y_en) = [np.Gnet(y,z) np.Bnet(y,z) 0 zeros(1,length(np.B{y})) y_au];  %encode components of bus y for (16) coup. const. of reactive power line          
                       
        E(z_st:z_en)= [-np.Gnet(y,z) -np.Bnet(y,z) 0 zeros(1,2*length(np.B{z}))]; %encode components of bus z for (16) coup. const. of reactive power line          
        
        % full matrix over time horizon for (16) coup. const. of reactive power line
        Ef{c} = kron(E,eye(np.h));                                          
        
        
    end
end
np.A_eq_dso = cat(1,Ef{:});
np.A_eq_dso_pinv = pinv(np.A_eq_dso);%'*np.A_eq_dso)*np.A_eq_dso'; 
%np.Aaug_dso_inv = inv([eye(size(np.A_eq_dso,2)) np.A_eq_dso'; np.A_eq_dso zeros(size(np.A_eq_dso,1))]);
%np.Aaug_dso_inv = np.Aaug_dso_inv(1:size(np.A_eq_dso,2),1:size(np.A_eq_dso,2));
np.n_eq_dso = c;
np.A_eq_dso_A = sparse(np.A_eq_dso_pinv*np.A_eq_dso);

end