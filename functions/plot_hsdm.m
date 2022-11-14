function plot_hsdm(o,o1,np)
figure
subplot(2,1,1)
hold on, grid on, box on

Pmg_t = np.sumPd;
Pmg_t1 = np.sumPd;
for i=1:np.n
   Pmg_t = Pmg_t + o.p_mg{i}(:,end);
   Pmg_t1 = Pmg_t1 + o1.p_mg{i}(:,end);
end
plot(Pmg_t,'LineWidth',1.5)
plot(Pmg_t1,'LineWidth',1.5)
plot(np.n*np.pmg_ref*ones(1,np.h),'k','LineWidth',2);

subplot(2,1,2)
hold on, grid on, box on


Ptr_t = zeros(24,1);
Ptr_t1 = zeros(24,1);
for i=1:np.n
    for j = 1:np.n
        if ~isempty(o.p_tr{i,j}) 
            for k=1:np.h

                    Ptr_t(k,1) = Ptr_t(k,1) + max(0,o.p_tr{i,j}(k,end));
                    Ptr_t1(k,1) = Ptr_t1(k,1) + max(0,o1.p_tr{i,j}(k,end));


            end
        end
    end
end


    
plot(Ptr_t,'LineWidth',1.5)
plot(Ptr_t1,'LineWidth',1.5)

end


