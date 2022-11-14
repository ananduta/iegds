% Simulations
% P2P market on IEGDS
% test centralized codes
% 03/10/2022


clear all
close all
clc

% Add path of folder 'functions'
addpath([pwd,'/functions'])
rng(61021)
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

c_max = 20;
for c = 1:c_max
    disp(['case number ',num2str(c)])
    %run('gen_iegs_33b_20n.m')
    run('gen_iegs_33b_20n.m')
    p_data{c} = p;
    
    
    %% Proposed method
    disp('Proposed method')
    %pwaReg = [20,40];
    pwaReg = [20,45];
    for rr = 1:length(pwaReg)
        disp(['no. of regions=',num2str(pwaReg(rr))])
        p.r = pwaReg(rr);
        p.gn.r = p.r;
        [~,oPWA{c,rr},qPWA{c,rr}] = alg2StageIterative(p);
    end
    
    % Proposed method with MISOCP
    [~,oMIS{c},qMIS{c}] = alg2StageSOCP(p);
    
    
    
    %% Fix flow directions based on proposed method
    for i = 1:p.n
        for jj = 1:length(p.gn.N{i})
            j = p.gn.N{i}(jj);
            
            %p.alphaFixed{i,j} = o{end}.delta_psi_s{i,j};
            
            p.alphaFixed{i,j} = abs(oPWA{c,end}.delta_psi_s{i,j} - ones(p.h,1));
            
            
            % check differences on binary decisions
            qMIS{c}.DFlowDir{i,j} = p.alphaFixed{i,j} - oMIS{c}.alpha_s{i,j};
        
            qMIS{c}.DFlow{i,j} = oPWA{c,end}.phi{i,j} - oMIS{c}.phi{i,j};
        end
            
    end
 
    
    %% SOCP method
    disp('SOCP method')
    pS = p;
    pS.fixedInt_flag = 1; pS.pen = 0;    
    gamma_pen = [0,50];
    for rr = 1:length(gamma_pen)
        disp(['penalty weight=',num2str(gamma_pen(rr))])
        pS.Gamma_pen = gamma_pen(rr);
        [~,oSOC{c,rr},qSOC{c,rr}] = algSOCP(pS);
        
    end
    rr= length(gamma_pen);
    
    %% SOCP-SCP method
    disp('SOCP method with SCP')
    p.fixedInt_flag = 1;
    [~,oSOC{c,rr+1},qSOC{c,rr+1}] = algSOCP_SCP(p);
       
    %% SOCP method with convex envelope
    disp('SOCP method with convex envelope')
    p.fixedInt_flag = 1; p.pen = 1;    
    p.Gamma_pen = 0;
    for i = 1:p.n
        for jj = 1:p.gn.noN(i)
            j = p.gn.N{i}(jj);
            p.gn.phim{i,j} = -p.gn.phi_max(i,j)*ones(p.h,1);
            p.gn.phiM{i,j} = p.gn.phi_max(i,j)*ones(p.h,1);
        end
    end
    [~,oSOC{c,rr+2},qSOC{c,rr+2}] = algSOCP_ce(p);
    
    
    %% SOCP method with refined convex envelope
    disp('SOCP method with refined convex envelope')
    p.fixedInt_flag = 1; p.pen = 1;    
    p.Gamma_pen = 0;
    
    % define the parameters of the affine function at each region
    p.r = pwaReg(end);
    phiM = p.gn.phi_max;
    phim = -p.gn.phi_max;
    clearvars('p.gn.phi_max');
    for i = 1:p.n
        for jj = 1:p.gn.noN(i)
            j=p.gn.N{i}(jj);
                        
            nf = @(y) y^2/p.gn.cf(i,j);
            
            pwaf = pwa_approx_nf(p.r,phim(i,j),phiM(i,j),nf);
            
            for h = 1:p.h
                
                for m = 1:p.r
                    if oPWA{c,end}.delta_s{i,j}{m}(h,1) == 1
                        p.gn.phim{i,j}(h,1) = pwaf.m(m);
                        p.gn.phiM{i,j}(h,1) = pwaf.M(m);
                    end
                end
            end
        end
    end
    
    
    [~,oSOC{c,rr+3},qSOC{c,rr+3}] = algSOCP_ce(p);
 
    %%
    save(['sim2_',date,'.mat'],'p_data','oPWA','qPWA','oMIS','qMIS','oSOC','qSOC')

end

