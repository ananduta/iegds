% IEGS Network generation
% random gas network

% W. Ananduta
% 18/06/2021

% Generate an integrated electrical and gas (IEGS) network
% Set the parameters of the network
% Set the gas and electrical loads

% time horizon
p.h = 1;
p.n = n_agents;


%% Electrical network (en)
% set time horizon
np.h= p.h;


% Adjacency matrix of trading network
np.n = p.n;
np.Adj = sparse(randconG(np.n,sp_set));

np.pas_ag = n_passive;

[np.N,np.noN] = id_neigh(np.Adj);
% assign components to each agent
n = np.n;
b = np.pas_ag;
np.b = b;
% type of load
% randomly assign the type of load profile   
np.t_lpr = randi([1 6], n+b,1);
% 0 = no load
%np.t_lpr(n+1) = 0;
%np.t_lpr(n+7) = 0;
%np.t_lpr(n+8) = 0;
    
% storage units (no storage units in [Sousa, et. al.,2019] and [Le Cadre, et al, 2019]
%np.st_un = zeros(n,1);
n_st=floor(np.n/2);
 ag_st = randperm(np.n,n_st); 
 np.st_un = zeros(n,1);
 for i = 1:n_st
     np.st_un(ag_st(i)) = 1;
end

% dispatchable units
%np.d_un = randi([0 1], n,1);
%np.d_un = [ones(1, 10) zeros(1,n-10)];

 frac_di = 0.25;
 n_di=floor(np.n*frac_di);
 ag_dg = randperm(np.n,n_di); 
 np.d_un = zeros(n,1);
 for i = 1:n_di
     np.d_un(ag_dg(i)) = 1;
end
%np.d_un = zeros(n,1);

% subset of gas-powered dispatchable units
frac_gu = 0.5;
n_gu = floor(n_di*frac_gu);
ag_gu = ag_dg(1:n_gu);
np.dgu_un = zeros(n,1);
for i = 1:n_gu
    np.dgu_un(ag_gu(i)) = 1;
end

% PV generation units
np.r_un = zeros(n+b,1);
for i=1:n
    if np.st_un == 1
        np.r_un(i) = randi([0 4]);
    end
end



% assign parameters in the local constraints
np = gen_param(np,ty); 



% generate power load and non-dispatchable profiles
[np.Pd,np.Pl,np.Pr] = gen_Pd_ext(np.n,np.t_lpr,np.r_un,np.pas_ag);

np.sumPd = sum(np.Pd(np.n+1:end,:))';



% initial condition of variables
np.init = 0;

% assign per-unit costs
np = gen_cost(np,tc); 

p.en = np;

clearvars('np');




%% Gas network (gn)
% number of nodes
gn.n = n_agents; 
sp = sp_set;

% Adjacency matrix of gas network

Adj = randconG(gn.n,sp);
G = graph(Adj);
G = minspantree(G);
Adj = adjacency(G);

gn.Adj = Adj;

[gn.N,gn.noN] = id_neigh(gn.Adj);
% Associate one node with an agent


gn.Ngu = ag_gu; % Set of nodes in Gas network with gas-powered units
gn.gfg = p.en.dgu_un; % Indicator whether the nodes have gas-powered units


% Set nodes connected to transmission network
%frac_di = 0.25;
% n_di=floor(np.n*frac_di);
n_gt = 1;
 ag_gt = randperm(gn.n,n_gt); 
 gn.gt = zeros(n,1);
 for i = 1:n_gt
     gn.gt(ag_gt(i)) = 1;
 end

% Set constraint parameters (pc)
gn = param_cons_ran(gn);



% generate gas demand
Gd_def = 30;
gn.Gdem = Gd_def*(1-0.2*rand(gn.n+n_passive,p.h));
gn.sumGd = sum(gn.Gdem(gn.n+1:end,:))';

% per-unit cost of gas
gn.d_gu = 0.1412*100/mean(1+sum(gn.Gdem(gn.n+1:end,:)));

% Set number of regions for PWA gas-flow model
p.r = 6;
gn.r = p.r;
p.gn = gn;
clearvars('gn');

% 
% %% Electrical network (37-bus IEEE network)
% en.n = 36;
% 
% en.Adj = zeros(en.n);
% en.B = zeros(en.n);
% en.G = zeros(en.n);
% 
% % operating voltage
% en.v_op = 4.8^2*1000; %kV
% 
% load('case_37bus.mat');
% tab = case_37bus;
% for i= 1:length(tab(:,1))
%     en.Adj(tab(i,1),tab(i,2)) = 1;
%     
%     if  tab(i,4) == 721
%         z_b = 0.2926 + 0.1973i;
%     elseif tab(i,4) == 722
%         z_b = 0.4751 + 0.2973i;
%     elseif tab(i,4) == 723
%         z_b = 1.2936 + 0.6713i;
%     elseif tab(i,4) == 724
%         z_b = 2.0952 + 0.7758i;
%     end
%     z= z_b*tab(i,3)/5280;
%     en.B(tab(i,1),tab(i,2)) = abs(imag(1/z));
%     en.G(tab(i,1),tab(i,2)) = abs(real(1/z));
% end
% en.Adj = en.Adj + en.Adj';
% [en.N,en.noN] = id_neigh(en.Adj);
% 
% en.B = en.B + en.B';
% en.B = en.v_op*en.B;
% en.G = en.G + en.G';
% en.G = en.v_op*en.G;
% 
% % en.rateA = 50*en.Adj;
% % en.rateA(4,5) = 60;
% % en.rateA(5,4) = 60;
% % en.rateA(10,11) = 60;
% % en.rateA(11,10) = 60;
% 
% % voltage limits
% en.v_max = 1.1;
% en.v_min = 0.9;
% en.theta_max = 15*pi/180;
% en.theta_min = -15*pi/180;
% 
% 
% % Set parameters of generators
% 
% % set details [pmax; q1; q2; q3]
% gen_det = [1000   800   500;
%            0.0004 0.001 0.005;
%            33.51  42.63 37.7;
%            176.95 129.97 137.41];           
% 
% % number of generators
% n_dg = n_gfg + floor(en.n*0.3);
% ag_dg = randperm(en.n,n_dg); 
% % type of generator (3 types, see gen_det)
% en.dg = zeros(en.n,1);
% for i = 1:n_dg
%     en.dg(ag_dg(i)) = randi([1 3]);
% end
% 
% % bounds
% en.pmin = zeros(en.n,1);
% en.pmax = zeros(en.n,1);
% 
% % cost function/ constraint parameters
% en.q1 = zeros(en.n,1);
% en.q2 = zeros(en.n,1);
% en.q3 = zeros(en.n,1);
% 
% 
% for i = 1:en.n
%     if en.dg(i) > 0
%         
%         en.pmax(i) = gen_det(1,en.dg(i));
%         en.q1(i)   = gen_det(2,en.dg(i));
%         en.q2(i)   = gen_det(3,en.dg(i));
%         en.q3(i)   = gen_det(4,en.dg(i));
%     end
% end
% 
% 
% % Set of nodes in Electrical network with gas-powered units 
% en.Ngu = ag_dg(1:n_gfg);
% 
% 
% % Set power load
% n_el=floor(en.n*0.7);
% ag = randperm(en.n,n_el); 
% en.el = zeros(en.n,1);
% for i = 1:n_el
%     en.el(ag(i)) = 1;
% end
% en.d = zeros(p.h,en.n);
% for i =1:en.n
%     if en.el(i) ==1
%         en.d(1:p.h,i) = 75*(1+ 0.4*rand(p.h,1));
%     end
% end
% 
% p.en = en;
% clearvars('en')
% 
% % Pair of nodes in gas and electrical network
% p.pair = [p.gn.Ngu; p.en.Ngu]';
          

%% Identify dimensions of decision variables (per time step, h)
for i = 1:p.n
p.nx(i) = 4+ p.en.noN(i);
p.ny(i) = 2+(2+p.gn.r)*p.gn.noN(i);
p.nz(i) = (1+3*p.gn.r)*p.gn.noN(i);
p.nu(i) = p.nx(i) + p.ny(i) + p.nz(i);
end
