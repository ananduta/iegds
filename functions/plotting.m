%% plot effect of storage units

figure
subplot(2,1,1)
dat = 'sim37b_tr2_sce_3.mat';
hold on, grid on, box on
load(dat);

Pdg_t = zeros(np.h,1);
for i=1:np.n
    Ni = sum(np.Adj(i,:));
    
    %compute a_ni,3
    a = zeros(3+Ni,1);
    a(1,1) = 1;
    Sdg = kron(eye(np.h),a');
    Pdg_t = Pdg_t + Sdg*s.u{i}(:,end);
end

Pdt = sum(np.Pd);

bar([s.sigma_mg(:,end) Pdg_t],'stacked')
plot(Pdt,'-o','color','k','LineWidth',1.5)
title('\textbf{Scenario a: No storage units}','Interpreter','latex')
ylabel('Power [KW]','Interpreter','latex')
%ylim([0, 800])

subplot(2,1,2)
hold on, grid on, box on,
dat = 'sim37b_tr2_sce_4.mat';
load(dat);

Pdg_t = zeros(np.h,1);
for i=1:np.n
    Ni = sum(np.Adj(i,:));
    
    %compute a_ni,3
    a = zeros(3+Ni,1);
    a(1,1) = 1;
    Sdg = kron(eye(np.h),a');
    Pdg_t = Pdg_t + Sdg*s.u{i}(:,end);
end

Pdt = sum(np.Pd);
bar([s.sigma_mg(:,end) Pdg_t],'stacked')
plot(Pdt,'-o','color','k','LineWidth',1.5)
title('\textbf{Scenario b: All agents with storage units}','Interpreter','latex')
ylabel('Power [KW]','Interpreter','latex')
xlabel('time step [hour]','Interpreter','latex')
legend({'$\hat{p}_k^{\mathrm{mg}}$','$\sum_{i\in \mathcal{N}} p_{i,k}^{\mathrm{dg}}$','$\sum_{i\in \mathcal{N}} p_{i,k}^{\mathrm{d}}$'},'Interpreter','latex'  )
%ylim([0, 800])

%{
subplot(3,1,2)
hold on, grid on, box on,
dat = 'sim3_n20_st2.mat';
load(dat);

Pdg_t = zeros(np.h,1);
for i=1:np.n
    Ni = sum(np.Adj(i,:));
    
    %compute a_ni,3
    a = zeros(3+Ni,1);
    a(1,1) = 1;
    Sdg = kron(eye(np.h),a');
    Pdg_t = Pdg_t + Sdg*s.u{i}(:,end);
end

Pdt = sum(np.Pd);
bar([s.ph_mg(:,end) Pdg_t],'stacked')
plot(Pdt,'-o','color','r','LineWidth',1.1)
title('\textbf{Scenario b: 8 agents with storage units}','Interpreter','latex')
ylabel('Power [kW]','Interpreter','latex')
ylim([0, 800])

%%



figure
subplot(2,1,1)
hold on, box on, grid on
load('sim2_main2.mat')
t_e = sim.t{1};
for l = 1:10
    plot([0.2:0.1:1],t_e(:,l),'-.','color',[.0, .4, .7])
end
plot([0.2:0.1:1],mean(t_e'),'-o','LineWidth',1.5,'color','r')
ylabel('Number of iterations','Interpreter','latex')
%xlabel('Graph connectivity','Interpreter','latex')
title(' \textbf{a. 10-agent network}','Interpreter','latex')

subplot(2,1,2)
hold on, box on, grid on
%load('sim2_main2.mat')
t_e = sim.t{2};
for l = 1:10
    plot([0.1:0.1:1],t_e(:,l),'-.','color',[.0, .4, .7])
end
plot([0.1:0.1:1],mean(t_e'),'-o','LineWidth',1.5,'color','r')
ylabel('Number of iterations','Interpreter','latex')
xlabel('Graph connectivity','Interpreter','latex')
title(' \textbf{b. 20-agent network}','Interpreter','latex')
%}