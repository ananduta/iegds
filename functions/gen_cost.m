function np = gen_cost(np,tc)
% Set the per-unit costs
% W. Ananduta
% 05/08/2021

% Inputs:
% np = struct
% ty = type of case

% Output:
% np = struct

% homogenous case

% per-unit cost of dispatchable generation
np.q_dg = zeros(1,np.n);
np.c_dg = zeros(1,np.n);

% per-unit cost of using storage
np.q_st = 0*ones(1,np.n);
np.c_st = 0*ones(1,np.n);

% per-unit cost of transferring power
np.c_tr = 0.02*100*np.Adj;

% per-unit cost of exchanging with main grid
np.d_mg = 0.1412/10/mean(1+sum(np.Pd(np.n+1:end,:)));
np.d_mg_l = 0.0005;
%np.d_mgs = 0.7/np.n;




for i=1:np.n
    if np.d_un(i) == 1
        % per-unit cost of dispatchable generation
        np.q_dg(i) = 0.01;%20*(1/np.t_lpr(i));
        np.c_dg(i) = 0.039*100;%20*(1/np.t_lpr(i));
        if np.dgu_un(i) ==1
            np.q_dg(i) = 0;
            np.c_dg(i) = 0;
        end
         
    end
    if np.st_un(i) == 0
        np.q_st(i) = 0;
        np.c_st(i) = 0;
    end
end

% no transaction costs
if tc==0
    np.c_tr = 0*np.c_tr;
end

% heterogenous transaction costs
if tc==2
    r_ctr = rand(np.n);
    r_ctr = r_ctr + r_ctr';
    np.c_tr = r_ctr;
end

% penalty of p_tr
np.q_tr = 0.01*100;
end