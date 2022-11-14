% Simulations
% P2P market on IEGDS
% test centralization
% 07/02/2021


clear all
close all
clc

% Add path of folder 'functions'
addpath([pwd,'/functions'])
%rng(22922)
%load('sim0911.mat')

% generate case
ty = [0]; %type of case study: (0)heterogenous  or (1)uniform  agents
tc = [1]; %uniform trading cost

% set the number of agents
n_agents = 33;
n_passive = 0;

r_max = 10;
sp_set = 0.7;

%run('gen_iegs_6n.m')
%run('gen_iegs_33b_20n.m')





%%

c_max = 1;
for c = 1:c_max
    disp(['case number ',num2str(c)])
    %run('gen_iegs_33b_20n.m')
    run('gen_iegs_33b_20n.m')
    p_data = p;
    
    
    %% Proposed method
    disp('Proposed method')
    %pwaReg = [20,40];
    pwaReg = [25];
    for rr = 1:length(pwaReg)
        disp(['no. of regions=',num2str(pwaReg(rr))])
        p.r = pwaReg(rr);
        [p1{rr},o{rr},q{rr}] = alg2StageIterative(p);
    end
    
    %% Proposed method with SOCP
    [~,oSOCP,qSOCP] = alg2StageSOCP(p);
    
    %% Fix flow directions based on proposed method
    for i = 1:p.n
        for jj = 1:length(p.gn.N{i})
            j = p.gn.N{i}(jj);
            
            %p.alphaFixed{i,j} = o{end}.delta_psi_s{i,j};
            
            p.alphaFixed{i,j} = abs(o{end}.delta_psi_s{i,j} - ones(p.h,1));
            
        end
            
    end
 
    
    %% SOCP method
    disp('SOCP method')
    pS = p;
    pS.fixedInt_flag = 1; pS.pen = 0;    
    gamma_pen = [0];
    for rr = 1:length(gamma_pen)
        disp(['penalty weight=',num2str(gamma_pen(rr))])
        pS.Gamma_pen = gamma_pen(rr);
        [~,os{rr},qs{rr}] = algSOCP(pS);
        
    end
    rr= length(gamma_pen);
    
    %% SOCP-SCP method
    disp('SOCP method with SCP')
    p.fixedInt_flag = 1;
    [~,os{rr+1},qs{rr+1}] = algSOCP_SCP(p);
       
    %% SOCP method with convex envelope
    disp('SOCP method with convex envelope')
    p.fixedInt_flag = 1; p.pen = 1;    
    p.Gamma_pen = 0;
%     [~,os{rr+2},qs{rr+2}] = algSOCP_ce(p);
    [~,osc,qsc] = algSOCP_ce(p);
    
%    save(['sim_',num2str(c),'_',date,'.mat'],'p_data','o','q','os','qs','oss','qss')



%% Solve MI-SOCP problem
disp('MISOCP method')
    pS.fixedInt_flag = 1; pS.pen = 0;    
    gamma_pen = [0];
    for rr = 1:length(gamma_pen)
        disp(['penalty weight=',num2str(gamma_pen(rr))])
        pS.Gamma_pen = gamma_pen(rr);
        [~,om{rr},qm{rr}] = algMISOCP(pS);
        
    end
    rr= length(gamma_pen);
    
   %% Fix flow directions based on MISOCP method
    for i = 1:p.n
        for jj = 1:length(p.gn.N{i})
            j = p.gn.N{i}(jj);
            
            p.alphaFixed{i,j} = om{rr}.alpha{i,j};
            
            %p.alphaFixed{i,j} = abs(o{end}.delta_psi_s{i,j} - ones(p.h,1));
            
        end
            
    end 
    
    %% SOCP-SCP method
    disp('SOCP method with SCP')
    p.fixedInt_flag = 1;
    [~,om{rr+1},qm{rr+1}] = algSOCP_SCP(p);
       
    %% SOCP method with convex envelope
    disp('SOCP method with convex envelope')
    p.fixedInt_flag = 1; p.pen = 0;   
    p.Gamma_pen = 0;
%     [~,os{rr+2},qs{rr+2}] = algSOCP_ce(p);
    [~,omc,qmc] = algSOCP_ce(p);
    
    
    %% MISOCP with convex envelope
    disp('SOCP method with convex envelope')
    p.fixedInt_flag = 1; p.pen = 0;    
    p.Gamma_pen = 0;
%     [~,os{rr+2},qs{rr+2}] = algSOCP_ce(p);
    [~,omc1,qmc1] = algMISOCP_ce(p);
end

