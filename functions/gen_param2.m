function np = gen_param(np,ty,Ts)
% Set the parameters in the constraints
% W. Ananduta
% 22/07/2019
% Edit for extension: 30/09/2020

% Inputs:
% np = struct
% ty = type of case

% Output:
% np = struct

if nargin < 3
    Ts=1;
end

% homogenous case study, ty==1
% local constraints
np.pdg_min = zeros(1,np.n);
np.pdg_max = 20*ones(1,np.n);
np.p_dh = 3*ones(1,np.n);
np.p_ch = 3*ones(1,np.n);
np.x_min = 20*ones(1,np.n);
np.x_max = 100*ones(1,np.n);

% storage unit
T = 60*60*Ts; % s
J2Wh = 1/3600;
np.a_st = ones(1,np.n);
np.b_st = T*J2Wh*ones(1,np.n);
np.x0 = np.x_min+50;

% power transfer constraints
np.pt_max = 0*np.Adj;

% power to main grid constraints
np.pmg_min = -30;
np.pmg_max = np.n*100;

% heterogenous case study, ty==0
if ty==0
    for i=1:np.n
        if np.d_un(i) ==1
            np.pdg_max(i) = np.t_lpr(i)^2*10;
        else
            np.pdg_max(i) = 0;
        end
        if np.st_un(i) == 0
            np.p_dh(i) = 0;
            np.p_ch(i) = 0;
            np.x0(i) = 0;
            np.x_min(i) = 0;
            np.x_max(i) = 0;
        else
            np.p_dh(i) = 2*np.t_lpr(i);
            np.p_ch(i) = 2*np.t_lpr(i);
            np.x_max(i) =100*(1+np.t_lpr(i)/max(np.t_lpr));
        end
            
    end
end

% parameters of the physical network
% MUST BE ADJUSTED PROPERLY TO ACCOUNT FOR LINEARIZATION
np.theta_max = 45/180*pi*ones(np.b,1);
np.theta_min = -45/180*pi*ones(np.b,1);
np.theta_max(1) = 0;
np.theta_min(1) = 0;
np.v_max = 1.05*ones(np.b,1);
np.v_min = 0.95*ones(np.b,1);
np.v_max(1) = 1;
np.v_min(1) = 1;
np.s_bar = 80^2*ones(np.b);

end