clc;clear;close all;

x_C = [0.004,0.015,0.07]*1e-2; X = 0:1e-5:1e-3;

f = figure(1); ax= gca;
ax.FontSize = 16;%set(gcf,'units','inches','position',[10,10,12,9])


%Volume diffusivity of Cr based on C conc
D0_V = [5.1e-4,1.0e-3,1.8e-2]*1e12;

subplot(2,2,1);
plot(x_C,(D0_V),'k*','MarkerSize',8);hold on;
P1 = polyfit(x_C,log(D0_V),1);
plot(X,exp(P1(1)*X+P1(2)),'r--','LineWidth',2 );
set(gca,'Yscale','log');
title('D^V_0');ylabel('D^V_0 (\mum^2/s)');xlabel('Carbon mass fraction');legend('Experimental data','Fit');
ax = gca; 
ax.FontSize = 16;

% figure(2);
subplot(2,2,2);
E0_V = [286e3,300e3,340e3];
plot(x_C,E0_V,'k*','MarkerSize',8);hold on;
P2 = polyfit(x_C,E0_V,1);
plot(X,P2(1)*X+P2(2),'r--','LineWidth',2);
title('E^0_V');ylabel('E^0_V (kJ/mol)');xlabel('Carbon mass fraction');legend('Experimental data','Fit');
ax = gca; 
ax.FontSize = 16;

%Grain boundary diffusivity of Cr based on C conc
D0_GB = [4.8e-12,1.9e-9,8.2e-8]* 1e18;
subplot(2,2,3);
plot(x_C,D0_GB,'k*','MarkerSize',8);hold on;
P3 = polyfit(x_C,log(D0_GB),1);
plot(X,exp(X*P3(1)+P3(2)),'r--','LineWidth',2);
set(gca,'Yscale','log');
title('D^{GB}_0');ylabel('D^{GB}_0 (\mum^2/s)');xlabel('Carbon mass fraction');legend('Experimental data','Fit');
ax = gca; 
ax.FontSize = 16;
E0_GB = [203e3,227e3,335e3];
subplot(2,2,4);
plot(x_C,E0_GB,'k*','MarkerSize',8);hold on;
P4 = polyfit(x_C,E0_GB,1);
plot(X,X*P4(1)+P4(2),'r--','LineWidth',2);
title('E^0_{GB}');ylabel('E^0_{GB} (kJ/mol)');xlabel('Carbon mass fraction');legend('Experimental data','Fit');
ax = gca; 
ax.FontSize = 16;

syms xc T; %carbon concentration
R = 8.314; %T = 973; 

D0 = P1(1);D1=P1(2);E0=P2(1);E1=P2(2);
D_V = simplify(exp(D0*xc+D1)*exp(-(E0*xc+E1)/R/T));%exp(P1(1)*xc+P1(2))*exp(-(P2(1)*xc+P2(2))/R/T);

D0 = P3(1);D1=P3(2);E0=P4(1);E1=P4(2);
D_GB = simplify(exp(D0*xc+D1)*exp(-(E0*xc+E1)/R/T));%exp(P3(1)*xc+P3(2))*exp(-(P4(1)*xc+P4(2))/R/T);

%Approach 1
d_GB = 5e-4;G=4.64
f = d_GB/G;
D = double(subs(subs(D_V*(1-f) + D_GB*f/d_GB,T,973),xc,0))

%Approach 2
f = 1/G;
D = double(subs(subs(D_V*(1-f) + D_GB*f,T,973),xc,0))

