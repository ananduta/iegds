function p = pwa_gasflow(p)
% pwa approximation of Weymouth gas-flow equation
% W. Ananduta
% 11/08/2021


gn = p.gn;

eps = 1e-7;

for i = 1:gn.n
    
    
    
    
    for jj = 1:gn.noN(i)
        j = gn.N{i}(jj);
        
        % define the parameters of the affine function at each region
        nf = @(y) y^2/gn.pc.cf(i,j);
        pwa = pwa_approx_nf(gn.r,-gn.pc.phi_max(i,j),gn.pc.phi_max(i,j),nf);
        gn.pwa = pwa;
        
        
        %(14)
        GFu = zeros(1,1+gn.noN(i));
        %GFu(2) = gn.pc.cf(i,j);
        GFu(1) = 1;%gn.pc.cf(i,j);
        
        GFz = zeros(1,(1+gn.r)*gn.noN(i));
        %GFz((1+gn.r)*(jj-1) + 1)= -2*gn.pc.cf(i,j);
        GFz((1+gn.r)*(jj-1) + 1)= -2;%*gn.pc.cf(i,j);
        
        GFd = zeros(1,(1+3*gn.r)*gn.noN(i));
        
        GFuj = zeros(1,1+gn.noN(j));
        %GFuj(2) = gn.pc.cf(i,j);
        GFuj(1) = 1;%gn.pc.cf(i,j);
        
        GFzj = zeros(1,(1+gn.r)*gn.noN(j));
        ii = find(gn.N{j}==i);
       % GFzj((1+gn.r)*(ii-1) + 1)= -2*gn.pc.cf(i,j);
        GFzj((1+gn.r)*(ii-1) + 1)= -2;%*gn.pc.cf(i,j);
               
        
        
        %(16)
        H16u = zeros(1,1+gn.noN(i));
        H16u(1) = -1;
        
        H16z = zeros(1,(1+gn.r)*gn.noN(i));
        
        
        H16d = zeros(1,(1+3*gn.r)*gn.noN(i));
        H16d((1+3*gn.r)*(jj-1) + 1) = -(gn.pc.psi_min(i)-gn.pc.psi_max(j));
        
        h16 = -(gn.pc.psi_min(i)-gn.pc.psi_max(j));
        
        H16j = zeros(1,1+gn.noN(j));
        H16j(1) = 1;
        
        %(17)
        H17u = zeros(1,1+gn.noN(i));
        H17u(1) = 1;
        
        H17z = zeros(1,(1+gn.r)*gn.noN(i));
        
        
        H17d = zeros(1,(1+3*gn.r)*gn.noN(i));
        H17d((1+3*gn.r)*(jj-1) + 1) = -(gn.pc.psi_max(i)-gn.pc.psi_min(j))-eps;
        
        h17 = -eps;
        
        H17j = zeros(1,1+gn.noN(j));
        H17j(1) = -1;
        
        Huij = [H16u;H17u];
        Hzij = [H16z;H17z];
        Hdij = [H16d;H17d];
        hj = [h16;h17];
        
        Hj = [H16j;H17j];
        
        gn.Hui{i,j} = kron(eye(p.h),Huij);
        gn.Hzi{i,j} = kron(eye(p.h),Hzij);
        gn.Hdi{i,j} = kron(eye(p.h),Hdij);
        gn.hij{i,j} = kron(ones(p.h,1),hj);
        
        gn.Huj{i,j} = kron(eye(p.h),Hj); %for u_j
        
        
        % Reciprocity constraints        
        GRu = zeros(1,1+gn.noN(i));
        GRu(1+jj) = 1;
        
        GRz = zeros(1,(1+gn.r)*gn.noN(i));
        
        GRd = zeros(1,(1+3*gn.r)*gn.noN(i));
        
        
        GRuj = zeros(1,1+gn.noN(j));
        ii = find(gn.N{j}==i);
        GRuj(1+ii) = 1;
        
        gn.GRcu{i,j} = kron(eye(p.h),GRu);
        gn.GRcz{i,j} = kron(eye(p.h),GRz);
        gn.GRcd{i,j} = kron(eye(p.h),GRd);
        gn.GRcuj{i,j} = kron(eye(p.h),GRuj);


        % Relates flow and delta_psi
        GRu1 = zeros(1,1+gn.noN(i));
        GRu1(1,1+jj) = 1;
        
        GRz1 = zeros(1,(1+gn.r)*gn.noN(i));
        
        GRd1 = zeros(1,(1+3*gn.r)*gn.noN(i));
        GRd1(1,(1+3*gn.r)*(jj-1)+1) = gn.pc.phi_max(i,j);
        
        gR1 = gn.pc.phi_max(i,j);
        
        
        GRu2 = zeros(1,1+gn.noN(i));
        GRu2(1+jj) = -1;
        
        GRz2 = zeros(1,(1+gn.r)*gn.noN(i));
        
        GRd2 = zeros(1,(1+3*gn.r)*gn.noN(i));
        GRd2((1+3*gn.r)*(jj-1)+1) = -gn.pc.phi_max(i,j)-eps;
        
        gR2 = -eps;
        
        GRu = [GRu1;GRu2];
        GRz = [GRz1;GRz2];
        GRd = [GRd1;GRd2];
        gR = [gR1;gR2];
        
        gn.GRu{i,j} = kron(eye(p.h),GRu);
        gn.GRz{i,j} = kron(eye(p.h),GRz);
        gn.GRd{i,j} = kron(eye(p.h),GRd);
        gn.gR{i,j} = kron(ones(p.h,1),gR);
        
        % define the simplex constraints of \delta_(i,j)^m
        
        
        Gsd = zeros(1,(1+3*gn.r)*gn.noN(i));
        gsd = 1;
        
        % define local constraints
        
        
        for m = 1:gn.r
            
            % indexes
            %delta
            id_delta_psi = (1+3*gn.r)*(jj-1) + 1;
            id_alpha_m = (1+3*gn.r)*(jj-1) + 1+3*(m-1) + 1;
            id_beta_m = (1+3*gn.r)*(jj-1) + 1+3*(m-1) + 2;
            id_delta_m = (1+3*gn.r)*(jj-1) + 1+3*(m-1) + 3;
            
            %(14) Gas-flow (cont'd)
            GFz((1+gn.r)*(jj-1) + 1 + m) = pwa.a(m);
            
            GFd(1,id_delta_m) = pwa.b(m);
            
            %(15) simplex constraint (cont'd)
            
            Gsd(1,id_delta_m) = 1;
            
            
            %(18)
            G1u = zeros(1,1+gn.noN(i));
            G1u(1+jj) = 1;
        
            G1z = zeros(1,(1+gn.r)*gn.noN(i));
            
            G1d = zeros(1,(1+3*gn.r)*gn.noN(i));
            G1d(1,id_alpha_m) = gn.pc.phi_max(i,j) - pwa.M(m);
            
            g1 = gn.pc.phi_max(i,j);
            
            %(19)
            G2u = zeros(1,1+gn.noN(i));
            G2u(1+jj) = -1;
        
            G2z = zeros(1,(1+gn.r)*gn.noN(i));
            
            G2d = zeros(1,(1+3*gn.r)*gn.noN(i));
            G2d(1,id_alpha_m) = -gn.pc.phi_max(i,j) - pwa.M(m)-eps;
            
            g2 = - pwa.M(m)-eps;
            
            %(20)
            G3u = zeros(1,1+gn.noN(i));
            G3u(1+jj) = -1;
        
            G3z = zeros(1,(1+gn.r)*gn.noN(i));
            
            G3d = zeros(1,(1+3*gn.r)*gn.noN(i));
            G3d(1,id_beta_m) = gn.pc.phi_max(i,j) + pwa.m(m);
            
            g3 = gn.pc.phi_max(i,j);
            
            %(21)
            G4u = zeros(1,1+gn.noN(i));
            G4u(1+jj) = 1;
        
            G4z = zeros(1,(1+gn.r)*gn.noN(i));
            
            G4d = zeros(1,(1+3*gn.r)*gn.noN(i));
            G4d(1,id_beta_m) = -gn.pc.phi_max(i,j) + pwa.m(m)-eps;
            
            g4 =  pwa.m(m)-eps;
            
            %(22)a
            G5u = zeros(1,1+gn.noN(i));
        
            G5z = zeros(1,(1+gn.r)*gn.noN(i));
            
            G5d = zeros(1,(1+3*gn.r)*gn.noN(i));
            G5d(1,id_alpha_m) = -1;
            G5d(1,id_delta_m) = 1;
            
            g5 =  0;
            
            %(22)b
            G6u = zeros(1,1+gn.noN(i));
        
            G6z = zeros(1,(1+gn.r)*gn.noN(i));
            
            G6d = zeros(1,(1+3*gn.r)*gn.noN(i));
            G6d(1,id_beta_m) = -1;
            G6d(1,id_delta_m) = 1;
            
            g6 =  0;
            
            %(22)c
            G7u = zeros(1,1+gn.noN(i));
        
            G7z = zeros(1,(1+gn.r)*gn.noN(i));
            
            G7d = zeros(1,(1+3*gn.r)*gn.noN(i));
            G7d(1,id_alpha_m) = 1;
            G7d(1,id_beta_m) = 1;
            G7d(1,id_delta_m) = -1;
            
            g7 =  1;
            
            %(23)a
            G8u = zeros(1,1+gn.noN(i));
        
            G8z = zeros(1,(1+gn.r)*gn.noN(i));
            G8z((1+gn.r)*(jj-1) + 1 + m) = -1;
            
            G8d = zeros(1,(1+3*gn.r)*gn.noN(i));
            G8d(1,id_delta_m) = -gn.pc.phi_max(i,j);
            
            g8 =  0;
            
            %(23)b
            G9u = zeros(1,1+gn.noN(i));
            G9u(1+jj) = -1;
        
            G9z = zeros(1,(1+gn.r)*gn.noN(i));
            G9z((1+gn.r)*(jj-1) + 1 + m) = 1;
            
            G9d = zeros(1,(1+3*gn.r)*gn.noN(i));
            G9d(1,id_delta_m) = gn.pc.phi_max(i,j);
            
            g9 =  gn.pc.phi_max(i,j);
            
            %(24)a
            G10u = zeros(1,1+gn.noN(i));
        
            G10z = zeros(1,(1+gn.r)*gn.noN(i));
            G10z((1+gn.r)*(jj-1) + 1 + m ) = 1;
            
            G10d = zeros(1,(1+3*gn.r)*gn.noN(i));
            G10d(1,id_delta_m) = -gn.pc.phi_max(i,j);
            
            g10 =  0;
            
            %(24)b
            G11u = zeros(1,1+gn.noN(i));
            G11u(1+jj) = 1;
        
            G11z = zeros(1,(1+gn.r)*gn.noN(i));
            G11z((1+gn.r)*(jj-1) + 1 + m) = -1;
            
            G11d = zeros(1,(1+3*gn.r)*gn.noN(i));
            G11d(1,id_delta_m) = gn.pc.phi_max(i,j);
            
            g11 =  gn.pc.phi_max(i,j);
            
            Gu{m}= [G1u;G2u;G3u;G4u;G5u;G6u;G7u;G8u;G9u;G10u;G11u];
            Gz{m}= [G1z;G2z;G3z;G4z;G5z;G6z;G7z;G8z;G9z;G10z;G11z];
            Gd{m}= [G1d;G2d;G3d;G4d;G5d;G6d;G7d;G8d;G9d;G10d;G11d];
            
            gm{m} = [g1;g2;g3;g4;g5;g6;g7;g8;g9;g10;g11];
            
            
            
        end
        Guj = cat(1,Gu{:});
        Gzj = cat(1,Gz{:});
        Gdj = cat(1,Gd{:});
        
        gj = cat(1,gm{:});
        
        %clearvars('Gu','Gz','Gd','gm');
        
        %(25)a
        G12u = zeros(1,1+gn.noN(i));

        G12z = zeros(1,(1+gn.r)*gn.noN(i));
        G12z((1+gn.r)*(jj-1) + 1) = -1;

        G12d = zeros(1,(1+3*gn.r)*gn.noN(i));
        G12d(1,id_delta_psi) = gn.pc.psi_min(i);

        g12 =  0;

        %(25)b
        G13u = zeros(1,1+gn.noN(i));
        G13u(1) = -1;

        G13z = zeros(1,(1+gn.r)*gn.noN(i));
        G13z((1+gn.r)*(jj-1) + 1) = 1;

        G13d = zeros(1,(1+3*gn.r)*gn.noN(i));
        G13d(1,id_delta_psi) = -gn.pc.psi_min(i);

        g13 =  -gn.pc.psi_min(i);

        %(26)a
        G14u = zeros(1,1+gn.noN(i));

        G14z = zeros(1,(1+gn.r)*gn.noN(i));
        G14z((1+gn.r)*(jj-1) + 1) = 1;

        G14d = zeros(1,(1+3*gn.r)*gn.noN(i));
        G14d(1,id_delta_psi) = -gn.pc.psi_max(i);

        g14 =  0;

        %(26)b
        G15u = zeros(1,1+gn.noN(i));
        G15u(1) = 1;

        G15z = zeros(1,(1+gn.r)*gn.noN(i));
        G15z((1+gn.r)*(jj-1) + 1) = -1;

        G15d = zeros(1,(1+3*gn.r)*gn.noN(i));
        G15d(1,id_delta_psi) = gn.pc.psi_max(i);

        g15 =  gn.pc.psi_max(i);
        

        
        
        Gua{jj} = [Guj;G12u;G13u;G14u;G15u];%G21u;G22u];
        Gza{jj} = [Gzj;G12z;G13z;G14z;G15z];%G21z;G22z];
        Gda{jj} = [Gdj;G12d;G13d;G14d;G15d];%G21d;G22d];
        
        ga{jj} = [gj;g12;g13;g14;g15];%g21;g22];
        
        %(14) Gas-flow eq. (cont'd)
        gn.GFui_eq{i,j} = kron(eye(p.h),GFu);
        gn.GFzi_eq{i,j} = kron(eye(p.h),GFz);
        gn.GFdi_eq{i,j} = kron(eye(p.h),GFd);
        gn.GFuj_eq{i,j} = kron(eye(p.h),GFuj);
        gn.GFzj_eq{i,j} = kron(eye(p.h),GFzj);
        
        %(15) Simplex const
        
        gn.Gsdi{i,j}= kron(eye(p.h),Gsd);
        gn.gsdi{i,j} = kron(ones(p.h,1),gsd);
        
    end
    Gu_a = cat(1,Gua{:});
    Gz_a = cat(1,Gza{:});
    Gd_a = cat(1,Gda{:});
    g_a = cat(1,ga{:});
    clearvars('Gua','Gza','Gda','ga');
    

% 
%     %(5)
%     G25u = zeros(1,1+gn.noN(i));
%     G25u(1) = 1;
% 
%     G25z = zeros(1,(1+gn.r)*gn.noN(i));
% 
%     G25d = zeros(1,(1+3*gn.r)*gn.noN(i));
% 
%     g25 =  gn.pc.psi_max(i);
%     
%     G26u = zeros(1,1+gn.noN(i));
%     G26u(1) = -1;
% 
%     G26z = zeros(1,(1+gn.r)*gn.noN(i));
% 
%     G26d = zeros(1,(1+3*gn.r)*gn.noN(i));
% 
%     g26 =  -gn.pc.psi_min(i);
    
    
    Gu_a = [Gu_a];%G23u;G24u;G25u;G26u];
    Gz_a = [Gz_a];%G23z;G24z;G25z;G26z];
    Gd_a = [Gd_a];%G23d;G24d;G25d;G26d];
    g_a = [g_a];%g23;g24;g25;g26];
    
    
    Gut{i} = kron(eye(p.h),Gu_a);
    Gzt{i} = kron(eye(p.h),Gz_a);
    Gdt{i} = kron(eye(p.h),Gd_a);
    
    gt{i} = kron(ones(p.h,1),g_a);
    
%     %(2) Gas balance equation
%     GBu = [1 0 ones(1,gn.noN(i))];
%     GBz = zeros(1,(1+gn.r)*gn.noN(i));
%     GBd = zeros(1,(1+3*gn.r)*gn.noN(i));
%     
%     gn.GBu_eq{i} = kron(eye(p.h),GBu);
%     gn.GBz_eq{i} = kron(eye(p.h),GBz);
%     gn.GBd_eq{i} = kron(eye(p.h),GBd);
%     gn.gb_eq{i} = gn.d(1:p.h,i);
end

gn.Gu = Gut;
gn.Gz = Gzt;
gn.Gd = Gdt;
gn.gt = gt;


p.gn = gn;
end