function s = initialize_uDSO(s,np)
% Update decisions of DSO
% 25-09-2020
% W. Ananduta

yalmip('clear')

cost=0;
cons=[];
for y=1:np.b
    u{y} = sdpvar(np.h*(3+2*length(np.B{y})),1);
    
    % cost function
    
   % cost = cost + 0.5*u{y}'*u{y} + (s.u_DSO{y}(:,t) - np.A_DSO{y}*a)'*u{y};
    
    cons = [cons, np.Aineq_dso{y}*u{y} <= np.bineq_dso{y}];
    
    % initialize p_tg =0
    A_au = [0 0 1 zeros(1,2*length(np.B{y}))];
    A_aux = kron(A_au,eye(np.h));
    cons = [cons, A_aux*u{y} == 0];
    
end
for y=1:np.b
    for zz = 1:length(np.B{y})
        z = np.B{y}(zz);
        
        % (15) eq. coupling const for real power line
        y_au = zeros(1,length(np.B{y}));
        y_au(zz) = -1;
        y_aux = [np.Bnet(y,z) -np.Gnet(y,z) 0 y_au zeros(1,length(np.B{y}))];
        y_auxi = kron(y_aux,eye(np.h));
        
        z_aux = [-np.Bnet(y,z) np.Gnet(y,z) 0 zeros(1,2*length(np.B{z}))];
        z_auxi = kron(z_aux,eye(np.h));
        
        cons =[cons, y_auxi*u{y} + z_auxi*u{z} == 0];
        
        % (16) eq. coupling const for reactive power line
        y_au = zeros(1,length(np.B{y}));
        y_au(zz) = -1;
        y_aux = [np.Gnet(y,z) np.Bnet(y,z) 0 y_au zeros(1,length(np.B{y}))];
        y_auxi = kron(y_aux,eye(np.h));
        
        z_aux = [-np.Gnet(y,z) -np.Bnet(y,z) 0 zeros(1,2*length(np.B{z}))];
        z_auxi = kron(z_aux,eye(np.h));
        
        cons =[cons, y_auxi*u{y} + z_auxi*u{z} == 0];
        
        
        % (17) line capacity constraints
        cons = [cons, u{y}(np.h*(3+zz-1)+1:np.h*(3+zz)).^2 + u{y}(np.h*(3+length(np.B{y})+zz-1)+1: np.h*(3+length(np.B{y})+zz)).^2 <= np.s_bar(y,z)^2*ones(np.h,1)];
    end
end
ops = sdpsettings('verbose',0,'solver','fmincon');
st=optimize(cons,cost,ops)

if st.problem ~= 0
%     disp(st.problem): in case solver cannot find the solution
    disp('not solved properly')
    st.problem
    pause
%     s.u{i}(:,t+1) = s.u{i}(:,t);
%     s.pmg{i}(:,t+1) = s.pmg{i}(:,t);
else
    %Assigning the decisions of each bus
    for y = 1:np.b
        s.u_DSO{y}(:,1) = value(u{y});
        s.p_tg{y}(:,1) = s.u_DSO{y}(np.h*2+1:np.h*3,1);
        for zz = 1:length(np.B{y})
            z = np.B{y}(zz);
            s.p_l{y,z}(:,1) =  s.u_DSO{y}(np.h*(3+zz-1)+1:np.h*(3+zz),1);
        end
    end
end

end