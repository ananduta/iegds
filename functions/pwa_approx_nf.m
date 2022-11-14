function pwa = pwa_approx_nf(r,m,M,nf)
% pwa approximation of a nonlinear function

% W. Ananduta
% 14/06/2021

% Input:
% r= number of evenly divided regions (non-negative regions)
% M= maximum
% m= minimum
% nf = nonlinear function

% Output:
% pwa.r : number of regions
% pwa.a : gradient for each region
% pwa.b : constant for each region


% define regions
pwa.l(1,1) = m;
for i=1:r
    pwa.l(i+1,1) = m + (M-m)*i/r;
end

pwa.m = pwa.l(1:end-1,1);
pwa.M = pwa.l(2:end,1);

x(1,1) = pwa.l(1,1);
y(1,1) = nf(x(1,1));
for i = 1:r
    x(i+1,1) = pwa.M(i);
%    y(i+1,1) = x(i+1,1)^2;
    y(i+1,1) = nf(x(i+1,1));
    pwa.a(i,1) = (y(i+1)-y(i))/(x(i+1)-x(i));
    pwa.b(i,1) = y(i)-pwa.a(i)*x(i);
   
end
pwa.r =r;




