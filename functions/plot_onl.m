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