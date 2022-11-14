% IEGS Network generation
% random gas network

% W. Ananduta
% 18/06/2021

% Generate an integrated electrical and gas (IEGS) network
% Set the parameters of the network
% Set the gas and electrical loads

% time horizon
p.h = 1;

%
%n_agents = 10;
p.n = n_agents;

% scaling
p.scale = 1;

%% Electrical network (en)
% set time horizon
np.h= p.h;

np.n = p.n;
maxBranch = 3;

np.Adj = randTreeGraph(p.n,maxBranch);
[np.N,np.noN] = id_neigh(np.Adj);

% Set B and G
% based on impedance in IEEE 13-bus network
np.Bnet = zeros(np.n);
np.Gnet = zeros(np.n);

data_tab3 = [300, 500, 800, 1000, 2000];

for i = 1:np.n
    for jj = 1:np.noN(i)
        j = np.N{i}(jj);
        if j > i
            tab4 = randi([601,607]);
            indexTab3 = randi([1,5]);
            tab3 = data_tab3(indexTab3);

            if  tab4 == 603
                z_b = 1.3294 + 1.3471i;
            elseif tab4 == 602
                z_b = 0.7526 + 1.1814i;
            elseif tab4 == 601
                z_b = 0.3465 + 1.0179i;
            elseif tab4 == 607
                z_b = 1.3425 + 0.5124i;
            elseif tab4 == 605
                z_b = 1.3292 + 1.3475i;
            elseif tab4 == 604
                z_b = 1.3238 + 1.3569i;
            elseif tab4 == 606
                z_b = 0.7982 + 0.4463i;
            end

            z= z_b*tab3/5280; % miles to feet
            np.Bnet(i,j) = abs(imag(1/z));
            np.Gnet(i,j) = abs(real(1/z));
        end
        
    end
end

np.v_op = 4.16^2*1000; %kV

np.scaling = 1e4;
np.Bnet = sparse(np.Bnet + np.Bnet');
np.Bnet = np.v_op*np.Bnet/np.scaling;
np.Gnet = sparse(np.Gnet + np.Gnet');
np.Gnet = np.v_op*np.Gnet/np.scaling;


np.pas_ag = 0;



% Set nodes connected to transmission network
%frac_di = 0.25;
% n_di=floor(np.n*frac_di);
n_et = 1;
 ag_et = 1; 
 np.et = zeros(np.n,1);
 for i = 1:n_et
     np.et(ag_et(i)) = 1;
 end
 
 
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
np.st_un = zeros(n,1);
% n_st=floor(np.n/2);
%  ag_st = randperm(np.n,n_st); 
%  np.st_un = zeros(n,1);
%  for i = 1:n_st
%      np.st_un(ag_st(i)) = 1;
% end

% dispatchable units
%np.d_un = randi([0 1], n,1);
%np.d_un = [ones(1, 10) zeros(1,n-10)];

 
 ag_dg = [1,3,4]; 
 np.d_un = zeros(n,1);
 for i = 1:length(ag_dg)
     np.d_un(ag_dg(i)) = 1;
end
%np.d_un = zeros(n,1);

% subset of gas-powered dispatchable units

ag_gu = [1,3];
np.dgu_un = zeros(n,1);
for i = 1:length(ag_gu)
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
np.scale=p.scale;
np = gen_param(np,ty); 



% generate power load and non-dispatchable profiles
%[np.Pd,np.Pl,np.Pr] = gen_Pd_ext(np.n,np.t_lpr,np.r_un,np.pas_ag);

Pd_data = [45, 60, 90, 100, 120, 150, 200, 210, 420]; % based on 33-bus system

np.Pd = zeros(p.h,p.n);
for i=1:np.n
    idx_Pd_data = randi([1,length(Pd_data)]);
    
    np.Pd(1:p.h,i) = Pd_data(idx_Pd_data)*ones(p.h,1);
end
np.Pd(1:p.h,1) = zeros(p.h,1);   % zero load on the pcc 

np.Pd = np.Pd';
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
%sp = sp_set;

gn.scale = p.scale;

% Adjacency matrix of gas network
maxBranch = 3;

gn.Adj = randTreeGraph(gn.n,maxBranch);
[gn.N,gn.noN] = id_neigh(gn.Adj);


% Associate one node with an agent


gn.Ngu = ag_gu; % Set of nodes in Gas network with gas-powered units
gn.gfg = p.en.dgu_un; % Indicator whether the nodes have gas-powered units


% Set nodes connected to transmission network
%frac_di = 0.25;
% n_di=floor(np.n*frac_di);
n_gt = 1;
 ag_gt = 1; 
 gn.gt = zeros(n,1);
 for i = 1:n_gt
     gn.gt(ag_gt(i)) = 1;
 end


% generate gas demand
gn.Gd_max = 6/gn.scale;
gn.Gdem = zeros(gn.n,p.h);

for i = 1:gn.n
    if rand <= 0.6
        gn.Gdem(i,:) = gn.Gd_max - 5*rand/gn.scale;
    end
end
gn.Gdem(1:p.h,1) = zeros(p.h,1);

gn.sumGd = sum(gn.Gdem(gn.n+1:end,:))';

% Set constraint parameters (pc)
gn = param_cons_gasNetwork_ran(gn);

% per-unit cost of gas
gn.d_gu = 0.1412*gn.scale*0.5/mean(1+sum(gn.Gdem(gn.n+1:end,:)));

% Set number of regions for PWA gas-flow model
p.r = 14;
gn.r = p.r;
p.gn = gn;
clearvars('gn');



%% Identify dimensions of decision variables (per time step, h)
for i = 1:p.n
p.nx(i) = 8+ p.en.noN(i);
p.nt(i) = p.gn.noN(i);
p.ny(i) = 2+(2+p.gn.r)*p.gn.noN(i)+p.nt(i);
p.nz(i) = (1+3*p.gn.r)*p.gn.noN(i);
p.nu(i) = p.nx(i) + p.ny(i) + p.nz(i);
end
p.id_dg = 1;
p.id_ch = 2;
p.id_dh = 3;
p.id_mg = 4;
p.id_gu = 5;
p.id_th = 6;
p.id_v = 7;
p.id_eg = 8;
