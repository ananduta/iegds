function Pr = gen_rprod(r_un,Ts,l)
% Generating a non-dispatchable production profile
% W. Ananduta
% 22/07/2019

% Inputs:
% r_un = type of non-dispatchable unit,
%           1 = solar
%           2 = wind
%           3 = both
% Ts = sampling time (in hour) <= 1
% l = generation length (in days)

% Output:
% load profile [1x (l*24/Ts) matrix]

% break conditions
if Ts > 1
    disp('Set Ts<=1')
    return
end

if Ts <= 0
    disp('Set Ts > 0 and Ts<=1')
    return
end

% Generating max PV profile
pvp = [0*ones(1,4/Ts) 0.05*ones(1,1/Ts) 0.1*ones(1,1/Ts) 0.15*ones(1,1/Ts) 0.25*ones(1,1/Ts) 0.4*ones(1,1/Ts) 0.6*ones(1,1/Ts)...
       0.8*ones(1,1/Ts) 0.9*ones(1,1/Ts) 0.85*ones(1,1/Ts) 0.8*ones(1,1/Ts) 0.7*ones(1,1/Ts) 0.5*ones(1,1/Ts) 0.3*ones(1,1/Ts)...
       0.15*ones(1,1/Ts) 0.05*ones(1,1/Ts) 0*ones(1,5/Ts)];

% Generating max Wind profile
wp = zeros(1,24/Ts);

% Range of maximum production
rmp(1,:) = [4 6]; %pv
rmp(2,:) = [10 15]; %wind



% Generate nominal production profile
if r_un == 1
    % Determine maximum production
    mp = rand*(rmp(r_un,2)-rmp(r_un,1))+rmp(r_un,1);
    % Nominal profile of 1 day
    ndp_1d = mp*pvp;
elseif r_un == 2
    % Determine maximum production
    mp = rand*(rmp(r_un,2)-rmp(r_un,1))+rmp(r_un,1);
    % Nominal profile of 1 day
    ndp_1d = mp*wp;
else
    mp_pv = rand*(rmp(1,2)-rmp(1,1))+rmp(1,1);
    mp_wp = rand*(rmp(2,2)-rmp(2,1))+rmp(2,1);
    ndp_1d = mp_pv*pvp + mp_wp*wp;
end
for ll = 1:l
    ndp(ll,:) = ndp_1d; 
end
Ndp = reshape(ndp',1,[]);

% Add some randomness (+/-5% of maximum production)
rl = 0.1*mp*(1-2*rand([1,length(Ndp)]))/2;
Pr = max(0,Ndp + rl);

end