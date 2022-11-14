
%%
figure
subplot(2,1,1)
hold on; grid on; box on;

load('sim_D_ag_27-Feb-2021.mat')
t = r.iter_tab(1:5,:)';
t_a = mean(t);
t_min = min(t);
t_max = max(t);
 plot([40:10:80],t_a,'-o','LineWidth',2.0,'color','b');
[ph,msg]=jbfill([40:10:80],t_min,t_max,[.0, .4, .7],[.0, .4, .7],0,0.5);

% for i=1:size(t,2)
%     plot([40:10:80],t(i,:),'-.','LineWidth',0.7,'color',[.0, .4, .7])
% end
%plot([40:10:90],t_a,'-o','LineWidth',1.5,'color','r')

ylabel('Number of iterations','Interpreter','latex')
legend({'Average'},'Interpreter','latex')
xlabel('Number of agents','Interpreter','latex')

%%
subplot(2,1,2)
hold on; grid on; box on;
load('sim_D_con_12-Mar-2021.mat')
t = r.iter_tab';
t_a = mean(t);
t_min = min(t);
t_max = max(t);
 plot([0.1:0.1:1],t_a,'-o','LineWidth',2.0,'color','b');
[ph,msg]=jbfill([0.1:0.1:1],t_min,t_max,[.0, .4, .7],[.0, .4, .7],0,0.5);
% t_a1 = mean(t1);
% plot([0.1:0.1:1],t_a1,'-o','LineWidth',2.0,'color','r')
% for i=1:size(t1,2)
%     plot([0.1:0.1:1],t1(i,:),'-.','LineWidth',0.7,'color',[.0, .4, .7])
% end

ylabel('Number of iterations','Interpreter','latex')
legend({'Average'},'Interpreter','latex')
xlabel('Connectivity','Interpreter','latex')
