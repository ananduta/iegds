% Main code 1
% P2P market, single simulation
% W. Ananduta
% 24/07/2019

clear all
close all
clc


% Generate case study
n = 20; %number of agents
sp = 0.5; %sparsity, in (0,1)
ty = 0; %type of case study: (0)heterogenous  or (1)uniform  agents
tc = 1; %transaction costs: (0)none, (1)uniform, or (2)heterogenous
np = gen_case(n,sp,ty,tc);

%% Distributed algorithm
%% INITIALIZATION

% Generate matrices S_i^mg and S_ij^tr
[np.Str,np.Smg] = gen_Smat(np);

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
np.t_max = 2000;
np.er_max = 0.01;

%% Iteration
while 1
    tic
    t
    % 1) Strategy update
%     s.ph_mg(:,t+1) = 0;
%     for i=1:np.n
        s = loc_opt_c(np,s,t);
        
        % Collect p_i^mg(t)
%         s.ph_mg(:,t+1) = s.ph_mg(:,t+1) + s.pmg{i}(:,t+1);
%     end
        
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
                                
                % 4) Auxiliary variable update:
                s.w{i}(:,t+1) = s.w{i}(:,t+1)+ np.gamma(i)*(s.lambda{i}(:,t)-s.lambda{j}(:,t));
                
            end
        end
      
        %5) Dual variable (grid constraints) update:        
        s.lambda{i}(:,t+1) = max(0, s.lambda{i}(:,t)+ np.delta(i)*(2*[np.Smg{i}*s.u{i}(:,t+1);-np.Smg{i}*s.u{i}(:,t+1)]...
                                    -[np.Smg{i}*s.u{i}(:,t);-np.Smg{i}*s.u{i}(:,t)]-[np.pmg_max/np.n*ones(np.h,1);-np.pmg_min/np.n*ones(np.h,1)]...
                                    -2*s.w{i}(:,t+1) + s.w{i}(:,t)));
                                
        % EXTRA: Compute dual residual                        
        dres2(t) = dres2(t) + (s.lambda{i}(:,t+1) - s.lambda{i}(:,t))'*(s.lambda{i}(:,t+1) - s.lambda{i}(:,t)); 
    end
    % Extra: compute residual
    s.res(t) = sqrt(res2(t));
    s.res(t)
    s.dres(t) = sqrt(dres2(t));
    s.dres(t)
    if s.res(t) < np.er_max && s.dres(t) < np.er_max && all(s.ph_mg(:,t) >= np.pmg_min) && all(s.ph_mg(:,t) <= np.pmg_max*ones(np.h,1)) 
        break
    end
    
    
    t= t+1;
    toc
    
end
%%
figure; plot(s.res);% set(gca,'yscale','log');
%%
figure;
subplot(3,1,1)
hold on, grid on, box on,
plot(s.ph_mg(:,end), 'LineWidth',1.1);
plot(20*ones(1,24),'k--','LineWidth',1.2)
title('total imported power from main grid')

subplot(3,1,2)
hold on, grid on, box on,
for i = 1:n
    plot(np.Smg{i}*s.u{i}(:,end),'LineWidth',1.1);
end
title('power trading with main grid')

subplot(3,1,3)
hold on, grid on, box on,
for tt = 1:24
    Pdt(1,tt) = sum(np.Pd(:,tt));
end
plot(Pdt,'LineWidth',1.1);
title('Total load')


%%
figure
subplot(2,1,1)
hold on, grid on, box on,
for i = 1:n
    
    Ni = sum(np.Adj(i,:));
    
    %compute a_ni,3
    a = zeros(3+Ni,1);
    a(2,1) = 1;
    
    %compute Smg
    Sst = kron(eye(np.h),a'); 
    
    
    plot(Sst*s.u{i}(:,end),'LineWidth',1.1);
end
title('power storage')
subplot(2,1,2)
hold on, grid on, box on,
for i = 1:n
    
    Ni = sum(np.Adj(i,:));
    
    %compute a_ni,3
    a = zeros(3+Ni,1);
    a(1,1) = 1;
    
    %compute Smg
    Sst = kron(eye(np.h),a'); 
    
    
    plot(Sst*s.u{i}(:,end),'LineWidth',1.1);
end
title('power dispatchable generator')
%%
figure
hold on, grid on, box on,
for i = 1:n
    Ni = sum(np.Adj(i,:));
    s.x(i,1) = np.x0(i);
    for tt = 1:24
        s.x(i,tt+1) = s.x(i,tt) - s.u{i}((3+Ni)*(tt-1)+2,end);
    end
    plot(s.x(i,:),'LineWidth',1.1);
end    
title('SoC')

%%
figure
hold on; box on; grid on;
title('total power transfer')
        pt_t = zeros(24,1);
        for i = 1:np.n
            for j= 1:np.n
                if j> i && np.Adj(i,j) ==1
                    pt = np.Str{i,j}*s.u{i}(:,end);
                    pt_t = pt_t + abs(pt);
                end
            end
        end
        
        
        plot(pt_t,'LineWidth',1.1);

    %legend('no cost','uniform costs','heteregenous costs')
    xlabel('k [hour]')
    ylabel('Power [kWh]')
    xlim([1 24])