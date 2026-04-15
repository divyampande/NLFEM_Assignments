
function [E,phi_0,alfa,phi_f,nu, R, m_w, rf, rho_0,rho_f,W_f1,W_f2, W_f3, W_f4,...
    E_a1, E_a2,  E_a3,  E_a4, A_01, A_02, A_03, A_04,n_1, n_2, n_3, n_4 ,...
    E_all, A_all, n_all, w0_all, wf_all] = material_prop_fun()



E=1.379*10^9;
phi_0=0.02;
alfa=1;
phi_f=0.2;
nu=0.2;
R=8.314; 
m_w=0.03;
rf=0.33; 
rho_0=1500;
rho_f=1300;
W_f1=0;
W_f2=0;
W_f3=0.290;
W_f4=0.190;
E_a1=88764.4;
E_a2=117236.0;
E_a3= 211443.5;
E_a4=272155.0;
A_01=1.207305*10^10; 
A_02= 4.0575*10^9;
A_03=3.857777*10^14;
A_04=5.583611*10^15;
n_1=3.5;
n_2=6.5;
n_3=6.5; 
n_4=3.3;
E_all=[88764.4, 117236.0, 211443.5, 272155.0 ];                      %activation energy for each reaction-
A_all=[1.207305*10^10, 4.0575*10^9, 3.857777*10^14, 5.583611*10^15]; % pre-exponential factor for each reaction
n_all=[ 3.5, 6.5, 6.5 3.3];                                          % order of chemical reaction for each reaction
w0_all=[0.015 0.095 0.590 0.300 ];
wf_all=[0  0  0.290, 0.190 ];
