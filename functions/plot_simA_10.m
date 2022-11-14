%% plot effect of storage units

%%
figure
subplot(2,1,1)
%dat = 'sim3_n10_st0.mat';
hold on, grid on, box on
%load(dat);
load('sim1.mat')
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
plot(Pdt,'-o','color','r','LineWidth',1.1)
title('\textbf{Scenario a: No storage units}','Interpreter','latex')
ylabel('Power [kW]','Interpreter','latex')

subplot(2,1,2)
hold on, grid on, box on,
load('sim2.mat')
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
plot(Pdt,'-o','color','r','LineWidth',1.1)
title('\textbf{Scenario c: All agents with storage units}','Interpreter','latex')
ylabel('Power [kW]','Interpreter','latex')
xlabel('time step [hour]','Interpreter','latex')
legend({'$\hat{p}_k^{\mathrm{mg}}$','$\sum_{i\in \mathcal{N}} p_{i,k}^{\mathrm{dg}}$','$\sum_{i\in \mathcal{N}} p_{i,k}^{\mathrm{d}}$'},'Interpreter','latex'  )

%%
subplot(3,1,2)
hold on, grid on, box on,
dat = 'sim3_n10_st2.mat';
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
title('\textbf{Scenario b: 4 agents with storage units}','Interpreter','latex')
ylabel('Power [kW]','Interpreter','latex')
