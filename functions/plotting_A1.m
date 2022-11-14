%% plot effect of storage units
clear all
close all

dat{1} = 'sim_jext_B_04-Feb-2021 _1.mat';
dat{2} = 'sim_jext_B_04-Feb-2021 _2.mat';
dat{3} = 'sim_jext_B_04-Feb-2021 _3.mat';
dat{4} = 'case_jext_03-Feb-2021.mat';

load(dat{4});
figure
subplot(4,1,1)

hold on, grid on, box on
load(dat{1});

Pdg_t = zeros(np.h,1);
Pdi_t = zeros(np.h,1);
Pmg_t = np.sumPd;
for i=1:np.n
   Pdi_t = Pdi_t + o.p_di{i}(:,end);
   Pmg_t = Pmg_t + o.p_mg{i}(:,end);
end

Pdt = sum(np.Pd);

bar([Pmg_t Pdi_t],'stacked')
plot(Pdt,'-o','color','k','LineWidth',1.5)
title('\textbf{Scenario a: No storage units}','Interpreter','latex')
ylabel('Power [KW]','Interpreter','latex')
%ylim([0, 800])

% subplot(3,1,2)
% hold on, grid on, box on,
% load(dat{3});
% 
% Pdg_t = zeros(np.h,1);
% Pdi_t = zeros(np.h,1);
% Pmg_t = np.sumPd;
% for i=1:np.n
%    Pdi_t = Pdi_t + o.p_di{i}(:,end);
%    Pmg_t = Pmg_t + o.p_mg{i}(:,end);
% end
% 
% Pdt = sum(np.Pd);
% 
% bar([Pmg_t Pdi_t],'stacked')
% plot(Pdt,'-o','color','k','LineWidth',1.5)
% title('\textbf{Scenario b: Half of agents with storage units}','Interpreter','latex')
% ylabel('Power [KW]','Interpreter','latex')
% xlabel('time step [hour]','Interpreter','latex')
%legend({'$\hat{p}_k^{\mathrm{mg}}$','$\sum_{i\in \mathcal{N}} p_{i,k}^{\mathrm{dg}}$','$\sum_{i\in \mathcal{N}} p_{i,k}^{\mathrm{d}}$'},'Interpreter','latex'  )
%ylim([0, 800])

subplot(4,1,2)
hold on, grid on, box on,
load(dat{3});

Pdg_t = zeros(np.h,1);
Pdi_t = zeros(np.h,1);
Pmg_t = np.sumPd;
for i=1:np.n
   Pdi_t = Pdi_t + o.p_di{i}(:,end);
   Pmg_t = Pmg_t + o.p_mg{i}(:,end);
end

Pdt = sum(np.Pd);

bar([Pmg_t Pdi_t],'stacked')
plot(Pdt,'-o','color','k','LineWidth',1.5)
title('\textbf{Scenario b: All agents with storage units}','Interpreter','latex')
ylabel('Power [KW]','Interpreter','latex')
xlabel('time step [hour]','Interpreter','latex')
legend({'total power from main grid','total power locally generated','total demand'},'Interpreter','latex'  )

subplot(4,1,3)
hold on, grid on, box on,

for c=[1]
    load(dat{c});

    Ptr_t = zeros(24,1);
    
    for i=1:np.n
        for j = 1:np.n
            if ~isempty(o.p_tr{i,j}) 
                for k=1:np.h
                    
                        Ptr_t(k,1) = Ptr_t(k,1) + max(0,o.p_tr{i,j}(k,end));
                        
                        
                end
            end
        end
    end

%    Pdt = sum(np.Pd);

    plot(Ptr_t,'LineWidth',1.5)
end    
%end
%plot(Pdt,'--o','color','k','LineWidth',1.5)
title('\textbf{P2P trading}','Interpreter','latex')
ylabel('Power [KW]','Interpreter','latex')
xlabel('time step [hour]','Interpreter','latex')
legend({'no agents with storage','all agents with storage'},'Interpreter','latex'  )

subplot(4,1,4)
hold on, grid on, box on
 load(dat{1})
 J1 = o.J/100;
 load(dat{3})
 J2 = o.J/100;
 DJ1 = (J2-J1);%./J2;
%  load(dat{2})
%  J3 = s.J;
%  load(dat{5})
%  J4 = s.J;
%  DJ2 = (J4-J3)./J4;
%  load(dat{3})
%  J5 = s.J;
%  load(dat{6})
%  J6 = s.J;
%  DJ3 = (J6-J5)./J6;
%bar([DJ1 DJ2 DJ3])
bar([DJ1])
%legend('no storage', 'half with storage', 'all with storage')
xlabel('Agents')
ylabel('Cost difference ')
title('Cost difference of Scenario a. - Scenario b.')
xlim([0 52])