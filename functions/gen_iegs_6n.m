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
n_agents = 6;
p.n = n_agents;

% scaling
p.scale = 100;

%% Electrical network (en)
% set time horizon
np.h= p.h;


% % Adjacency matrix of trading network
% % A=1, B=2, C=4,D=6, E=5, F=3,
% 
 np.n = p.n;
% 
% Adj = zeros(np.n);
% 
% Adj(1,2) = 1;
% Adj(1,6) = 1;
% Adj(2,6) = 1;
% Adj(2,4) = 1;
% Adj(3,4) = 1;
% Adj(3,5) = 1;
% Adj(5,6) = 1;
% Adj = Adj + Adj';
% 
% np.Adj = Adj;


% IEEE 13-bus feeder

load('case_13bus.mat');
np.v_op = 4.16^2*1000; %kV
tab = lines_13_bus;
np.Adj = zeros(11);
np.Bnet = np.Adj;
np.Gnet = np.Adj;
for i= 1:length(tab(:,1))
    np.Adj(tab(i,1),tab(i,2)) = 1;
    
    if  tab(i,4) == 603
        z_b = 1.3294 + 1.3471i;
    elseif tab(i,4) == 602
        z_b = 0.7526 + 1.1814i;
    elseif tab(i,4) == 601
        z_b = 0.3465 + 1.0179i;
    elseif tab(i,4) == 607
        z_b = 1.3425 + 0.5124i;
    elseif tab(i,4) == 605
        z_b = 1.3292 + 1.3475i;
    elseif tab(i,4) == 604
        z_b = 1.3238 + 1.3569i;
    elseif tab(i,4) == 606
        z_b = 0.7982 + 0.4463i;
    end
    z= z_b*tab(i,3)/5280; % miles to feet
    np.Bnet(tab(i,1),tab(i,2)) = abs(imag(1/z));
    np.Gnet(tab(i,1),tab(i,2)) = abs(real(1/z));
end
np.Adj = sparse(np.Adj + np.Adj');
np.scaling = 1e4;
p.scaling = np.scaling;
np.Adj = np.Adj(1:p.n,1:p.n);
np.Bnet = sparse(np.Bnet + np.Bnet');
np.Bnet = np.v_op*np.Bnet/np.scaling;
np.Gnet = sparse(np.Gnet + np.Gnet');
np.Gnet = np.v_op*np.Gnet/np.scaling;

np.pas_ag = 0;

[np.N,np.noN] = id_neigh(np.Adj);

% Set nodes connected to transmission network
%frac_di = 0.25;
% n_di=floor(np.n*frac_di);
n_et = 1;
 ag_et = 3; 
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
[np.Pd,np.Pl,np.Pr] = gen_Pd_ext(np.n,np.t_lpr,np.r_un,np.pas_ag);

np.Pd = zeros(p.h,p.n);
np.Pd(1:p.h,3) = 75;
np.Pd(1:p.h,4) = 100;
np.Pd(1:p.h,5) = 100;
np.Pd(1:p.h,1) = 75;
np.Pd(1:p.h,2) = 100;
np.Pd(1:p.h,6) = 100;

np.Pd = np.Pd';
np.Pd = (0.5+0.7*rand)*np.Pd;
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

gn.scale = p.scale;

% Adjacency matrix of gas network

Adj = zeros(gn.n);
Adj(1,2) = 1;
Adj(2,4) = 1;
Adj(2,5) = 1;
Adj(3,5) = 1;
Adj(5,6) = 1;

Adj = Adj + Adj';

gn.Adj = Adj;

[gn.N,gn.noN] = id_neigh(gn.Adj);
% Associate one node with an agent


gn.Ngu = ag_gu; % Set of nodes in Gas network with gas-powered units
gn.gfg = p.en.dgu_un; % Indicator whether the nodes have gas-powered units


% Set nodes connected to transmission network
%frac_di = 0.25;
% n_di=floor(np.n*frac_di);
n_gt = 1;
 ag_gt = 6; 
 gn.gt = zeros(n,1);
 for i = 1:n_gt
     gn.gt(ag_gt(i)) = 1;
 end

% Set constraint parameters (pc)
gn = param_cons_6n(gn);



% generate gas demand
Gd_def = 250/gn.scale;
gn.Gdem = zeros(gn.n+n_passive,p.h);
gn.Gdem(1, :) = Gd_def;%*(1-0.2*rand(1,p.h));
gn.Gdem(3, :) = Gd_def;%*(1-0.2*rand(1,p.h));
%gn.Gdem(2, :) = Gd_def;%*(1-0.2*rand(1,p.h));
gn.Gdem(4, :) = Gd_def;%*(1-0.2*rand(1,p.h));
%gn.Gdem(5, :) = Gd_def*(1-0.2*rand(1,p.h));
%gn.Gdem(6, :) = Gd_def*(1-0.2*rand(1,p.h));
gn.sumGd = sum(gn.Gdem(gn.n+1:end,:))';

% per-unit cost of gas
gn.d_gu = 0.1412*gn.scale*0.5/mean(1+sum(gn.Gdem(gn.n+1:end,:)));
gn.d_gu_l = 0.001;
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
