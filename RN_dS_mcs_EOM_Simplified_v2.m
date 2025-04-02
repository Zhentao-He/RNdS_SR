function k = RN_dS_mcs_EOM_Simplified_v2(q,mu,D,R,init,D_tmp,L2,La)
% 本函数计算RN_dS时空下带电标量场的运动方程
% msc为massive charged scalar的缩写
% v2代表该版本使用B'演化B

N = length(R);

psi = init(1:N);
Pi = init(N+1:2*N);
chi = init(2*N+1:3*N);
H = init(3*N+1:4*N);

psiDz = D*psi;
tmp = H.*Pi + psiDz;
tmp2 = real( Pi.*conj(psiDz) );

B_z1 = init(end);
Bp = - 2*q*L2./R.^2.*imag(Pi.*conj(psi));
B = D_tmp\Bp;
B = B - B(end) + B_z1;
kB_z1 = -2*q*L2*imag(conj(psi(end)).*tmp(end))./R(end).^2;


% Gauge
Ap = B.*exp(-chi);
A = D_tmp\Ap;
%A = A - A(1);
 A = A - A(end);

kpsi = Pi + H.*psiDz + 1i*q*A.*psi;

kPi = -mu^2*L2*exp(-chi)./R.^2.*psi + 1i*q*A.*Pi + D*(tmp)  - 2*tmp./R;

kchi = 2*H.*(D*chi) + 4*R .* tmp2 ;

 switch nargin 
     case 6+1
        kH = (1-H.^2).*(D*chi) - 2*R.* ( abs(Pi).^2 + abs(psiDz).^2 + 2*H.*tmp2);
     case 7+1
        kH = -D*(1-H.^2) + 3*(1-H.^2)./R - R.*exp(-chi) + exp(-chi)./R .* ( La +2*mu^2*abs(psi).^2 ) + R.^3.*exp(-chi).*B.^2;
     otherwise
         error('输入参数个数错误')
end

k = [kpsi;kPi;kchi;kH;kB_z1];
end