clc

clear 

close all



%import files

A = importdata("/davinci-1/work/openfoam/rhoEnergyFoam/TestCaseV2006/TBL_NASA_M02_REFOpt_ModeA_CFL1_Spalart_WF/postProcessing/sampleDict/0.029/sample_T_p_rho_thermo:mu.xy");

B = importdata("/davinci-1/work/openfoam/rhoEnergyFoam/TestCaseV2006/TBL_NASA_M02_REFOpt_ModeA_CFL1_Spalart_WF/postProcessing/sampleDict/0.029/sample_U.xy");



%A = importdata("x02neg_T_p_rho_thermo:mu_RCF.xy");

%B = importdata("x02neg_U_RCF.xy");



%A = importdata("x02neg_T_p_rho_thermo:mu_NASA9.xy");

%B = importdata("x02neg_U_NASA9.xy");



%A = importdata("grandezze_ml_3_6.xy");

%B = importdata("U_ml_3_6.xy");



%A = importdata("grandezze_ml_1.xy");

%B = importdata("U_ml_1.xy");



%rename variables

y (:,:) = B(:,1);

%y = y_b-y_b(1);



U_x (:,:) = B(:,2);     %U component on x  

%U_x = U_xb-U_xb(1);

U_y (:,:) = B(:,3);     %U component on y

U_z (:,:) = B(:,4);     %U component on z

    

T   (:,:) = A(:,2);     %Temperature    

p   (:,:) = A(:,3);     %Pressure   

rho (:,:) = A(:,4);     %density    

mu  (:,:) = A(:,5);     %dynamic viscosity

nu = mu./rho;           %kinematic viscosity



%defining gradient

%delta_u = -22*U_x(1)+36*U_x(2)-18*U_x(3)+4*U_x(4);

%delta_y = -22*y(1)+36*y(2)-18*y(3)+4*y(4);



delta_u = U_x(2)-U_x(1);

delta_y = y(2)-y(1);



dudy = delta_u/delta_y;



%wall shear stress 

tau_w = mu.*dudy;



%defining shear velocity u_tau

u_tau = sqrt(abs(tau_w)./rho);



%defining functions in wall units

U = sqrt(U_x.^2 + U_y.^2 + U_z.^2);

yPlus = (y.*u_tau)./nu;



U_Plus = (1./0.41)*log(yPlus)+5;                %log law

uPlus = U./u_tau;

%uPlus = uPlus/2.5;



%friction coefficient

L_x = linspace(0,1.06,39);

C_f = tau_w ./ (0.5*rho.*(U).^2);



%defining Reynolds at wall

Re_tau = u_tau/nu;



%evaluate delta99

L = length(U_x);

u_99 = 0.99*U_x(L);

u_98 = 0.98*U_x(L);

F = find(U_x < u_99);

L_f = length(F);

delta99 = y(L_f);



%Plots

figure(1)

%a0 = plot(U_x, y, LineWidth=3);

semilogx(yPlus,uPlus, LineWidth=3);

hold on

semilogx(yPlus,U_Plus, LineWidth=3);

hold on

semilogx(yPlus, yPlus, LineWidth=3);

legend('TBL','Log law', 'U+ = y+', Location='northwest')

xlabel("y+")

ylabel("U+")

xlim([0.1 10^3])

ylim([0 30])

title('Law of the Wall')

grid on



% figure(2)
% 
% plot(L_x, C_f, LineWidth=3)
% 
% xlabel("plate length")
% 
% ylabel("Cf")
% 
% xlim([0 2])
% 
% ylim([0 1])
% 
% title('Friction coefficient')
% 
% grid on
