function np = gen_case(n,sp,ty,tc)
% Generate a case study
% 19/07/2019
% W. Ananduta

% Inputs:
% n = number of agents
% sp = sparsity of the network [in (0,1)]
% ty = type of case study: (0)heterogenous  or (1)uniform  agents
% tc = transaction costs: (0)none, (1)uniform, or (2)heterogenous

% Output:
% np = network's parameters (struct, see below for details) 

%% assign default value of the arguments:
% sparsity=1, ty=1, tc=1
if nargin<2
    sp = 1;
    ty = 1;
    tc = 1;
end

if nargin<3
    ty = 1;
    tc = 1;
end

if nargin<4
    tc = 1;
end

%% break conditions
if ty > 1
    disp('choose type of case study: 0=heterogenous or 1=uniform');
    return
end

if tc > 2
    disp('choose type of transaction costs: 0=none, 1=uniform or 2=heterogenous');
    return
end



%% generate the network randomly
np.n = n;
np.Adj = randconG(n,sp);

% homogenous case study: small residential agents, with PV and storage unit
np.t_lpr = ones(n,1);
np.r_un = ones(n,1);
np.d_un = zeros(n,1);
np.st_un = ones(n,1);

% generate heterogenous case study
if ty==0
    % randomly assign the type of load profile
    % 1 = small residential
    % 2 = large residential
    % 3 = industrial
    np.t_lpr = randi([1 3], n,1);

    % randomly assign renewable generation units
    % 1 = only PV
    % 2 = only wind
    % 3 = both PV and wind
    %np.r_un = randi([1 3], n,1);
    np.r_un = ones(n,1);
    
    % randomly assign dispatchable units (only industrial or large residential 
    % agents might have dispatchable units)
    % 0 = no dispatchable unit
    % 1 = has dispatchable unit
    np.d_un = randi([0 1],n,1);
    for i = 1:n
        if np.t_lpr(i) == 0
            np.d_un(i) = 0;
        end
    end
    
    % randomly assign storage units
    % 0 = no storage unit
    % 1 = has storage unit
    np.st_un = randi([0 1],n,1);
end

% assign parameters in the local constraints
np = gen_param(np,ty);

% assign per-unit costs
np = gen_cost(np,tc);

% generate load and non-dispatchable profiles
[np.Pd,np.Pl,np.Pr] = gen_Pd(np.n,np.t_lpr,np.r_un);
np.h = size(np.Pd,2); % time horizon
end


