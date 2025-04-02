% ------------------------------------------------------------
% Nonlinear Evolution of RN-dS with Hyperboloidal Formalism
%                       Zhen-Tao He
%                    arXiv: 2411.03193
% ------------------------------------------------------------
clear; close all; 
format long
%% 1. parameters setting
%  1.1 parameters of black hole 
M = 0.2; % mass
La = 1e-3/M^2; % cosmological constant
Q = 0.3*M; % charge

%  1.2 parameters of scalar field
q = 0.87/Q; % charge
% q = 1;
mu = 2e-4/M; % mass
k = 7.964e-3; % amplitude

%  1.3 parameters of scalar field

B_mode = 1; % 1表示用B'演化B，否则用\dot{B}
initial_type = "Gauss"; % "kZ","poly" ,"Gauss"
RK4_flag = 0; % Set 1 to evolve via myRK4 instead of ode45 

% ------------------------------------------------------
% 检查M,Q的值是否合适
if 9 - 8*Q^2/M^2 < 0
   error('The parameters M and Q are inappropriate!')
end
% 计算L^2
L2 = compute_L2(M,Q,La);
% 检查k值是否足够小
if k/sqrt(Q/q/L2)>0.01
    warning('The amplitude k is a bit large!')
    % pause
end
% 场B和Q差个L2,还有个负号
B = -Q/L2;
% ======== horizon positions of of RN-dS ======
% 作为求解初始视界位置的迭代初值
% as the iterative initial data 
% for solving the initial horizon position
z0_tmp = L2*(3*M - sqrt(9*M^2 - 8*Q^2))/4/Q^2; % H的零点
syms x
Horizon_eq = x^2/L2 - 2*M*x^3/L2^2 + Q^2*x^4/L2^3 - La*L2/3 ==0;
Horizon_tmp = double(solve(Horizon_eq));
Co_H_tmp = Horizon_tmp(2); %宇宙学视界
Ap_H_tmp = Horizon_tmp(3); %表观视界
Ca_H_tmp = Horizon_tmp(4); %柯西视界
% ======================================================
width = (Ap_H_tmp - Co_H_tmp)/5; % 波包宽度
% 波包中心位置
tmp = 0.5;
centre = tmp*Ap_H_tmp + (1-tmp)*Co_H_tmp; 
% used of compact support 'ploy'
zl = centre - width/2;
zr = centre + width/2;
% =========== 求解初始Horizon位置 & 确定演化左右边界 ==========
% 根据ψ初值类型指定f(z)
switch initial_type
    case "kZ"
        fun = @(z) exp(-k^2*z.^2)/L2 .* (L2^2*La./z.^4 + (L2^2*2*mu^2*k^2-1)./z.^2 + B^2);
    case "poly"
        psi_init_fun(x) = piecewise( (x>=zl) & (x<=zr) , k*x.^3.*(1-x/zl).^3.*(x/zr-1).^3 , 0);
        psi_init_fun_tmp(x) = k*x.^3.*(1-x/zl).^3.*(x/zr-1).^3;
        psidz_init_fun(x) = diff(psi_init_fun_tmp,x);
        chi_init_fun_tmp(x) = int(2*x.*psidz_init_fun(x).^2,x);
        chi_init_fun(x) = piecewise((x>=zl) & (x<=zr) ,int(2*x.*psidz_init_fun(x).^2,x),x>zr,chi_init_fun_tmp(zr),x<zl,chi_init_fun_tmp(zl));
        fun = @(z) double(exp(-chi_init_fun(z))) /L2 .* (L2^2*(La+2*mu^2*abs(double(psi_init_fun(z)).^2))./z.^4  -1./z.^2 + B^2);
    case "Gauss"
        psi_init_fun(x) = k*exp(-((x-centre)/width).^2);
        psidz_init_fun(x) = diff(psi_init_fun,x);
        chi_init_fun(x) = int(2*x.*psidz_init_fun(x).^2,x);
        fun = @(z) double(exp(-chi_init_fun(z))) /L2 .* (L2^2*(La+2*mu^2*abs(double(psi_init_fun(z)).^2))./z.^4  -1./z.^2 + B^2);
    otherwise
        error("Please choose the right type of initial data")
end

% 确定演化左右边界 with given initial scalar perturbation
% See Appendix.A2
tic
z0 = fzero(@(z)1/z^3 + z*fun(z)/3,z0_tmp);
C_fun = @(z1) - integral(fun,z1,z0,'RelTol',1e-12,'AbsTol',1e-16) - z0*fun(z0)/3;
Co_H0 = fzero(C_fun,Co_H_tmp );
z1 = 0.9*Co_H0; %左边界
C = C_fun(z1);
z2_fun = @(z2) C + integral(fun,z1,z2,'RelTol',1e-12,'AbsTol',1e-16);
Ap_H0 = fzero(z2_fun,1.01*Ap_H_tmp);
z2 = 1.01*Ap_H0; %右边界
toc
% Cauchy_H0 = fzero(z2_fun,1.01*Cauchy_H_tmp);
% z2 = 1.01*max(Cauchy_H0,Cauchy_H_tmp);
% CH0/AH0
% Mh_init = L2/2/Ap_H_tmp; % initial M_irr

% Only when initial scalar perbation is far enough for the BH,
% M almost equal to M_init [1610.08352].
% M_init = Mh_init + Q^2/4/Mh_init -4*La/3*Mh_init.^3
%% 2.Grids
% Z方向格点参数
RN = 80; 
N = RN+1; 
[D,R] = cheb(RN); RL=z2-z1; % 演化区间长度
R = (RL/2)*(R+1)+z1; % R的区间 
D = D/(RL/2); %切比雪夫谱求导矩阵
% D_tmp = D(1:N-1,1:N-1);
D_tmp=D; D_tmp(end,:)=zeros(1,RN+1); D_tmp(end,end)=1;

% =========== 初值=========== 
H_init = RNdS_H_init(R,z0,C,fun);
switch initial_type
    case "kZ"
        chi_init = k^2*R.^2;
        psi_init = k*R;
        Pi_init = -H_init*k;
    case {"poly","Gauss"}
        psi_init = double(psi_init_fun(R));
        Pi_init = -H_init.*D*psi_init;
        chi_init = double(chi_init_fun(R));
    otherwise
        error("Please choose the right type of initial data")
end
if B_mode ~=1
    B_init = B*ones(RN+1,1);
end
clear fun
% =========== checks before evolutions =========== 
% Mass function
% Mass_init = .5*L2./R .*(1 + Q^2*R.^2/L2^2 - L2^2*La/3./R.^2 - L2*(1-H_init.^2)./R.^2 );
% plot(R,Mass_init)

% 计算rescaled Misner sharp质量
[Mass_MS,M_variation] = compute_MS(L2,R,chi_init,H_init,M,Q,La,z1,z2,Co_H0,Co_H_tmp);

% 检查初值约束
psiDz_init = D*psi_init;
kchi = 2*H_init.* (D*chi_init) + 4*R.*real(Pi_init.*conj(D*psi_init));
kH = (1 - H_init.^2).*(D*chi_init) - 2*R.*(abs(Pi_init).^2 + abs(psiDz_init).^2 + 2*H_init.*real(Pi_init.*conj(psiDz_init)));
constraint2 = (1-H_init.^2).*(D*chi_init)./R  - 2*H_init.*(D*H_init)./R - 3*(1-H_init.^2)./R.^2 + exp(-chi_init)/L2 ...
            - (...
              2*(abs(Pi_init).^2 + abs(psiDz_init).^2 + 2*H_init.*real(Pi_init.*conj(psiDz_init)) ) + L2*exp(-chi_init)./R.^2.*(La + 2*mu^2*abs(psi_init).^2) + R.^2.*exp(-chi_init).*B^2/L2 ...
              );

% convergence_initial_check(psi_init,chi_init,H_init)


%% 3.Evolution
% disp('按任意键开始演化');pause
disp('开始演化！')
% 方程演化
tic
if RK4_flag == 1
    dT = 0.5*(R(1)-R(2))/2;
    T_cyc = ceil(1/dT);
    TL = 4000;
    init = [psi_init;Pi_init;B_init;chi_init;H_init];
    [T,sol] = myRK4(@(T,init)RN_dS_mcs_EOM_Simplified(q,mu,D,R,init,D_tmp,L2),TL,T_cyc,dT,init);
    psi_sol = sol(1:N,:);
    Pi_sol = sol(N+1:2*N,:);
    B_sol = sol(2*N+1:3*N,:);
    chi_sol = sol(3*N+1:4*N,:);
    H_sol = sol(4*N+1:5*N,:);
else
    dT = 50; 
    TL = 5000; % ode45得到的解可能依赖于TL
    T = [0:99 100:dT:TL];
    options = odeset('InitialStep',1e-4,'RelTol',1e-8,'AbsTol',1e-14);
    if B_mode==1
        init = [psi_init;Pi_init;chi_init;H_init;B];
        [T,sol] = ode45(@(T,init)RN_dS_mcs_EOM_Simplified_v2(q,mu,D,R,init,D_tmp,L2),T,init,options);
        psi_sol = sol(:,1:N).';
        Pi_sol = sol(:,N+1:2*N).';
        chi_sol = sol(:,2*N+1:3*N).';
        H_sol = sol(:,3*N+1:4*N).';
        B_z1_sol = sol(:,end).';
        Bp_sol = -2*q*L2./R.^2.*imag(Pi_sol.*conj(psi_sol));
        B_sol = D_tmp\Bp_sol;
        B_sol = B_sol - B_sol(end,:) + B_z1_sol;
        clear B_z1_sol
    else
        init = [psi_init;Pi_init;B_init;chi_init;H_init];
        [T,sol] = ode45(@(T,init)RN_dS_mcs_EOM_Simplified(q,mu,D,R,init,D_tmp,L2),T,init,options);
        psi_sol = sol(:,1:N).'; % '对于复数数组是复共轭转置！！！
        Pi_sol = sol(:,N+1:2*N).';
        B_sol = sol(:,2*N+1:3*N)';
        chi_sol = sol(:,3*N+1:4*N)';
        H_sol = sol(:,4*N+1:5*N)';
    end
end
toc
clear sol
disp('演化完成。')

% solution of the gauge field A
% Ap_sol = B_sol.*exp(-chi_sol);
% A_sol = D_tmp\Ap_sol;
% A_sol = A_sol - A_sol(end,:);
% plot(R,A_sol(:,end))

%% 4.calution after evolutions

% some scripts:
% RNdS_error_check 
% spectral_convergence_check
% RNdS_Quantities
% RNdS_psi_h

% clear unnecessary data before save
clear C_fun z2_fun fun H_fun Q_fun Qtot_fun rho_fun
% save data 
% for the critical phenomena project 
save("muM_2e-4_LaM2_1e-3_qQ.87_k_"+num2str(k*1e3,'%g')+"e-3.mat")
%% 5.local funcions
function L2 = compute_L2(M,Q,La) 
Delta = sqrt(9*M^2 - 8*Q^2);
L2 = - 96*Q^6 / ( 3*(27*M^4-36*M^2*Q^2 + 8*Q^4 - 9*M^3*Delta + 8*M*Q^2*Delta) + 32*La*Q^6 );
% Q很小时舍入误差非常大,适当情况下用泰勒展开近似计算
if Q < 0.1*M
    L2_tmp = 27*M^2/(1-9*M^2*La) - 9*Q^2/(1-9*M^2*La)^2 - (1-36*M^2*La)*Q^4/M^2/(1-9*M^2*La)^3;
    % 判断直接计算的结果相对近似计算结果是否偏差太大
    if abs( (L2_tmp - L2)/L2_tmp ) > 0.5
        warning("L2直接计算的值"+num2str(L2)+"被取代为"+num2str(L2_tmp));
        L2 = L2_tmp;
    end
end
end
function [Mass_MS,M_variation] = compute_MS(L2,R,chi_init,H_init,M,Q,La,z1,z2,Co_H0,Co_H_tmp) 
    Mass_MS = L2/2./R .* (1 - L2*exp(chi_init).*(1-H_init.^2)./R.^2- La*L2^2/3./R.^2) ; 
    MS_coef = real_to_cheb(Mass_MS); 
    BH_Mass = (M-Q^2/2/L2*Co_H_tmp);
    % 质量增加的值
    M_variation = cheb_interpolate(MS_coef,z1,z2,Co_H0) - BH_Mass;
    % 切换成质量增加的百分比
    M_variation = M_variation/BH_Mass;
    if M_variation > 0.1
        error("The initial perturbation is too large, whose mass ratio is "+num2str(M_variation*100)+"%.")
    end
end
function H = RNdS_H_init(Z,z0,C,fun)
H = 0*Z;
Z1 = min(Z);
for i = 1:length(Z)
    if Z(i) < z0
        H(i) = sqrt(1 - ( C + integral(fun,Z1,Z(i)) ) *Z(i)^3 );
    else
        H(i) = - sqrt(1 - ( C + integral(fun,Z1,Z(i)) ) *Z(i)^3 );
    end
end
end
function [] = convergence_initial_check(psi_init,chi_init,H_init)
    close all
    figure
    plot(abs(real_to_cheb(psi_init)));set(gca, 'YScale', 'log')
    figure
    plot(abs(real_to_cheb(H_init)));set(gca, 'YScale', 'log')
    figure
    plot(abs(real_to_cheb(chi_init)));set(gca, 'YScale', 'log')
end
function [T,sol] = myRK4(fun,TL,T_cyc,dT,init)

% fun(T,x) 给出导数（列向量），x也为列向量形式
% dT时间步长，每T_cyc个步长(周期)记录数据，运行TL个周期
T = 0:TL;
T = T*T_cyc*dT;

NN = length(init);
sol = zeros(NN,TL+1);
sol(:,1) = init;
for ti = 1:TL
    for tj = 1:T_cyc
        K1 = fun( ((ti-1)*T_cyc+(tj-1))*dT , init );
        K2 = fun( ((ti-1)*T_cyc+(tj-1)+0.5)*dT , init+K1*dT/2 );
        K3 = fun( ((ti-1)*T_cyc+(tj-1)+0.5)*dT , init+K2*dT/2 );
        K4 = fun( ((ti-1)*T_cyc+(tj-1)+1)*dT , init+K3*dT );
        init = init + dT/6*(K1+2*K2+2*K3+K4);
    end
    sol(:,ti+1) = init;
end
end