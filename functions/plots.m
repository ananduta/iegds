%% Imported power from main grid
m = ['s','*','o'];
figure
title('Homogeneous agents')
subplot(1,2,1)
hold on; box on; grid on;

plot(20*ones(1,24),'k--','LineWidth',1.2)
for i = 1:2
    %for ii = 1:3
        dat = ['sim_cs',num2str(i),'_ty',num2str(1),'_tc',num2str(1)];
        
        
        load(dat);
        plot(s.ph_mg(:,end),'Marker',m(i),'LineWidth',1);
    %end
end
legend('lower bound','n=20','n=10','n=50')
xlabel('k [hour]')
ylabel('Power [kWh]')
xlim([1 24])
subplot(1,2,2)
hold on; box on; grid on;
title('Heterogeneous agents')
plot(20*ones(1,24),'k--','LineWidth',1.2)
for i = 1:2
    %for ii = 1:3
        dat = ['sim_cs',num2str(i),'_ty',num2str(0),'_tc',num2str(1)];
        
        
        load(dat);
        plot(s.ph_mg(:,end),'Marker',m(i),'LineWidth',1);
    %end
end
%legend('lower bound','n=20','n=10','n=50')
xlabel('k [hour]')
ylabel('Power [kWh]')
xlim([1 24])
%% total power exchanged
m = ['s','*','o'];
figure
subplot(1,2,1)
hold on; box on; grid on;
title('Homogeneous agents')
% for c = 1:3
    for cc = 1:3
        dat = ['sim_cs',num2str(1),'_ty',num2str(1),'_tc',num2str(cc-1)];
        load(dat);
        pt_t = zeros(24,1);
        for i = 1:np.n
            for j= 1:np.n
                if j> i && np.Adj(i,j) ==1
                    pt = np.Str{i,j}*s.u{i}(:,end);
                    pt_t = pt_t + abs(pt);
                end
            end
        end
        
        
        plot(pt_t,'Marker',m(cc),'LineWidth',1);
    end
    legend('no cost','uniform costs','heteregenous costs')
    xlabel('k [hour]')
    ylabel('Power [kWh]')
    xlim([1 24])

subplot(1,2,2)
hold on; box on; grid on;
title('Heterogeneous agents')
% for c = 1:3
    for cc = 1:3
        dat = ['sim_cs',num2str(1),'_ty',num2str(0),'_tc',num2str(cc-1)];
        load(dat);
        pt_t = zeros(24,1);
        for i = 1:np.n
            for j= 1:np.n
                if j> i && np.Adj(i,j) ==1
                    pt = np.Str{i,j}*s.u{i}(:,end);
                    pt_t = pt_t + abs(pt);
                end
            end
        end
        
        
        plot(pt_t,'Marker',m(cc),'LineWidth',1);
    end
    legend('no cost','uniform costs','heteregenous costs')
    xlabel('k [hour]')
    ylabel('Power [kWh]')
    xlim([1 24])
% end

%%
%m = ['s','*','o'];
cs =2;
figure
subplot(1,2,1)
hold on; box on; grid on;
title('Homogeneous agents')
dat = ['sim_cs',num2str(cs),'_ty',num2str(1),'_tc',num2str(0)];
load(dat);
J_b = s.J;
% for c = 1:3
    for cc = 1:2
        dat = ['sim_cs',num2str(cs),'_ty',num2str(1),'_tc',num2str(cc)];
        load(dat);
               
        dJ(cc,:) =s.J-J_b;

    end
    bar(dJ');
    legend('uniform costs','heteregenous costs')
    xlabel('Agent')
    ylabel('Cost')
    xlim([0 np.n+1])

   subplot(1,2,2)
hold on; box on; grid on;
title('Heterogeneous agents')
dat = ['sim_cs',num2str(cs),'_ty',num2str(0),'_tc',num2str(0)];
load(dat);
J_b = s.J;
% for c = 1:3
    for cc = 1:2
        dat = ['sim_cs',num2str(cs),'_ty',num2str(0),'_tc',num2str(cc)];
        load(dat);
               
        
        dJ(cc,:) =s.J-J_b;

    end
    bar(dJ');
    
    legend('uniform costs','heteregenous costs')
    xlabel('Agent')
    ylabel('Cost')
    xlim([0 np.n+1])