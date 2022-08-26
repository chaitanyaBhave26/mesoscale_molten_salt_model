%Uses an ideal solution fit to CALPHAD data and electrode potentials to
%calculate the equilibrium concentrations of Ni and Cr in FCC alloy +
%molten FLiBe system

clc; close all; clear;
syms x_Ni_metal x_Ni_melt x_Cr_metal x_Cr_melt;
syms w_Ni w_Cr;
kB = 8.617333262145e-5; %eV/K
T = 700+273; %K
F = 96485.33212; %C/mol Faraday's constant
Gxs =  -1.56448695e+04 +1.56011217*T;%3.6122583e+03 + 4.45532057e+00*T;%


G0_f_Ni_metal = (-5179.159 + 117.854*T - 22.096*T*log(T) - (4.8407e-3)*T^2)/F; %eV/atom
G0_f_Cr_metal = (-1572.94 + 157.643*T - 26.908*T*log(T) + 1.89435E-3*T^2 - 1.47721E-6*T^3 + 139250/T + Gxs)/F; %eV/atom
G0_f_Va_metal =  1.56 - 3.3*kB*T;

k_metal = 1;%1.770%1.89;

F_metal = G0_f_Ni_metal*x_Ni_metal + G0_f_Cr_metal*x_Cr_metal + G0_f_Va_metal*(1-x_Ni_metal-x_Cr_metal) + k_metal*kB*T*(x_Ni_metal*log(x_Ni_metal) + x_Cr_metal*log(x_Cr_metal) + (1-x_Ni_metal-x_Cr_metal)*log(1-x_Ni_metal-x_Cr_metal) );

E0_Cr = -0.39;
E0_Ni = 0.473;
E0_F  = 2.871;
E_F = -(-155*4184)/2/F; %3.64; %V -> F/F- potential of salt in experiment
n = 2;
k_melt = 1;

G0_f_Cr_Cr = (-8856.94 + 157.48*T - 26.908*T*log(T) + 1.89435E-3*T^2 - 1.47721E-6*T^3 + 139250/T)/F;
G0_f_Ni_melt = G0_f_Ni_metal + n*(E0_Ni + (-E0_F+E_F) );
G0_f_Cr_melt = G0_f_Cr_Cr  + n*(E0_Cr + (-E0_F+E_F) );

F_melt = G0_f_Ni_melt*x_Ni_melt + G0_f_Cr_melt*x_Cr_melt + k_melt*kB*T*(x_Ni_melt*log(x_Ni_melt) + x_Cr_melt*log(x_Cr_melt));

omega_metal = G0_f_Va_metal - kB*T*k_metal*log( 1 + exp( (w_Ni - (G0_f_Ni_metal  - G0_f_Va_metal))/kB/T/k_metal ) + exp( (w_Cr - (G0_f_Cr_metal - G0_f_Va_metal))/kB/T/k_metal ) );
c_Ni_metal  = exp( (w_Ni - (G0_f_Ni_metal - G0_f_Va_metal))/kB/T/k_metal )/( 1 + exp( (w_Ni - (G0_f_Ni_metal - G0_f_Va_metal))/kB/T/k_metal ) + exp( (w_Cr - (G0_f_Cr_metal - G0_f_Va_metal))/kB/T/k_metal ));
c_Cr_metal  = exp( (w_Cr - (G0_f_Cr_metal - G0_f_Va_metal))/kB/T/k_metal )/( 1 + exp( (w_Ni - (G0_f_Ni_metal - G0_f_Va_metal))/kB/T/k_metal ) + exp( (w_Cr - (G0_f_Cr_metal - G0_f_Va_metal))/kB/T/k_metal ));

c_Ni_melt = exp(-1+ (w_Ni - G0_f_Ni_melt )/kB/T/k_melt );
c_Cr_melt = exp(-1+ (w_Cr - G0_f_Cr_melt )/kB/T/k_melt );
omega_melt = simplify(-kB*T*k_melt*(c_Ni_melt+c_Cr_melt) );

del_omega = omega_metal - omega_melt;

% mu_Ni_vals = [];%zeros(N);
% mu_Cr_vals = [];%zeros(N);
% Ni_metal_vals = [];%zeros(N);
% Cr_metal_vals = [];%zeros(N);
% Ni_melt_vals = [];%zeros(N);
% Cr_melt_vals = [];%zeros(N);
% 
% %
% N = 100; min = -1.0; max = 1.0;
% for i = 1:N
%     mu = min + (max-min)*i/N;
%     F_temp =subs(del_omega,w_Ni,mu);
% %     assume(w_Cr > -10.0 & w_Cr < 10.0);
%     S = vpasolve(F_temp==0);
%     if (isempty(S)==0 && abs(vpa(S)) < 10.0 )
%         mu_Ni_vals = [mu_Ni_vals,mu]; mu_Cr_vals = [mu_Cr_vals,vpa(S)];
%         Ni_metal_vals = [Ni_metal_vals,vpa(subs(subs(c_Ni_metal,w_Ni,mu),w_Cr,vpa(S)) )];
%         Cr_metal_vals = [Cr_metal_vals,vpa(subs(subs(c_Cr_metal,w_Ni,mu),w_Cr,vpa(S)) )];
%         Ni_melt_vals = [Ni_melt_vals,vpa(subs(subs(c_Ni_melt,w_Ni,mu),w_Cr,vpa(S)) )];
%         Cr_melt_vals = [Cr_melt_vals,vpa(subs(subs(c_Cr_melt,w_Ni,mu),w_Cr,vpa(S)) )];
% %         plot(mu,vpa(S),'*');
%     end
% end
% 
% plot(mu_Ni_vals,mu_Cr_vals,'r-','LineWidth',2);
% title('Equilibrium \mu_{Cr} vs \mu_{Ni}');
% xlabel('\mu_{Ni} (eV/atom)');
% ylabel('\mu_{Cr} (eV/atom)');
% ax= gca;
% ax.FontSize = 16;
% 
% figure(2);
% hold on;
% plot(mu_Ni_vals,Ni_metal_vals,'r-','LineWidth',2);
% plot(mu_Ni_vals,Cr_metal_vals,'b-','LineWidth',2);
% plot(mu_Ni_vals,1-Ni_metal_vals-Cr_metal_vals,'k-','LineWidth',2);
% 
% title('Equilibrium mole fractions of Ni and Cr in metal phase vs \mu_{Ni}');
% xlabel('\mu_{Ni}');
% ylabel('Component mole fraction');
% legend({'x_{Ni}','x_{Cr}','x_{Va}'});
% ax= gca;
% ax.FontSize = 16;
% set(gca, 'YScale', 'log')
% 
% hold off;
% 
% figure(3);
% hold on;
% plot(mu_Ni_vals,Ni_melt_vals,'r-','LineWidth',2);
% plot(mu_Ni_vals,Cr_melt_vals,'b-','LineWidth',2);
% title('Equilibrium mole fractions of Ni and Cr in molten FLiBe vs \mu_{Ni}');
% xlabel('\mu_{Ni}');
% ylabel('Component mole fraction');
% legend({'x_{Ni^{2+}}','x_{Cr^{2+}}'});
% 
% ax= gca;
% ax.FontSize = 16;
% set(gca, 'YScale', 'log')
% hold off



x_Va = exp(-(G0_f_Va_metal)/kB/T); %Equilibrium vacancy concentration in Ni

x_Cr = 0.0561;
x_Ni = 1-x_Cr-x_Va;

ni_flibe = 3.6948e-11; %Same chemical potential as Ni in Ni-5Cr
cr_flibe = 25e-6;

mu_ni_metal = vpa(subs(subs(diff(F_metal,x_Ni_metal),x_Ni_metal,x_Ni ),x_Cr_metal,x_Cr))
mu_cr_metal = vpa(subs(subs(diff(F_metal,x_Cr_metal),x_Ni_metal,x_Ni ),x_Cr_metal,x_Cr))
mu_ni_melt = vpa(subs(subs(diff(F_melt,x_Ni_melt),x_Ni_melt, ni_flibe ),x_Cr_melt,cr_flibe))
mu_cr_melt = vpa(subs(subs(diff(F_melt,x_Cr_melt),x_Ni_melt, ni_flibe ),x_Cr_melt,cr_flibe))
%
del_omega_metal = double(subs(subs(del_omega,w_Ni,mu_ni_metal),w_Cr,mu_cr_metal))
del_omega_melt = double(subs(subs(del_omega,w_Ni,mu_ni_melt),w_Cr,mu_cr_melt))

% del_omega_int = del_omega_metal-del_omega_melt

% h = fsurf(F_metal,[0 1 0 1],'EdgeColor','r');
% camlight(110,70)
% % brighten(0.6)
% % h.EdgeColor = 'none';
% h.AmbientStrength = 0.6;
% % colormap winter;
%
% hold on;
% h = fsurf(F_melt,[0 1 0 1],'EdgeColor','b');
% camlight(110,70)
% % brighten(0.6)
% % h.EdgeColor = 'none';
% h.AmbientStrength = 0.6;
% colormap hot;
