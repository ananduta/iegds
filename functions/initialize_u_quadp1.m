function  s = initialize_u_quadp1(np,s,i)
% Initialization of the algorithm
% W. Ananduta
% 24/09/2020

% Inputs:
% np: Network parameters
% s: simulation results
% i: agent
% t: iteration



%for i = 1:np.n
    
    A = np.A_ineq{i};
    b = np.b_ineq{i};
    Aeq = np.Aeq{i};
    beq = np.beq{i};
    options = optimset('Display','off');

    [u_i,fval,exitflag] = quadprog(zeros(size(A,2)),zeros(size(A,2),1),A,b,Aeq,beq,[],[],[],options);

    if exitflag ~= 1 && exitflag ~= 0
    %     disp(st.problem): in case solver cannot find the solution
        disp('not solved properly')
        st.problem
        pause
    %     s.u{i}(:,t+1) = s.u{i}(:,t);
    %     s.pmg{i}(:,t+1) = s.pmg{i}(:,t);
    else
        %Assigning the decisions of agent i
        
        %dim_prev_ag = 0;
        %for i = 1:np.n
         %   dim_i = np.h*(3+sum(np.Adj(i,:)));
          %  u_i = u_all(dim_prev_ag+1:dim_prev_ag+dim_i);
        s.u{i}(:,1) = u_i;
        s.p_di{i}(:,1) = np.Sdi{i}*u_i;
        s.p_st{i}(:,1) = np.Sst{i}*u_i;
        s.p_mg{i}(:,1) = np.Smg{i}*u_i;
        s.p_mgs{i}(:,1) = np.Smgs{i}*u_i;
        for jj = 1:length(np.N{i})
            j = np.N{i}(jj);
            s.p_tr{i,j} = np.Str{i,j}*u_i;
        end
            
           % dim_prev_ag = dim_prev_ag+dim_i;
        %end
    end

%end

