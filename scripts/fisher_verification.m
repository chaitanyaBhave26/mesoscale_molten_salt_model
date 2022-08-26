clear;clc;close all;
ColorSet =lines(10);
figure(9)
figure9 = figure(9);
axes9 = axes('Parent',figure9,...
'Box'         , 'on'     , ...
'TickDir'     , 'in'     , ...
'TickLength'  , [.02 .02] , ...
'XMinorTick'  , 'off'      , ...
'YMinorTick'  , 'off'      , ...
'YGrid'       , 'off'      , ...
'XColor'      , 'k', ...
'YColor'      , 'k', ...
'FontWeight' , 'bold'     ,...
'FontSize'   , 14,...
'LineWidth'   , 1.5   );
set(gca, 'ColorOrder', ColorSet);
hold all
d_GB = sym(5e-4); %0.5 um
syms f;

%%Ni-5Cr
D = 2.072e-7;D_GB = (7.26e-4/d_GB);   %Literature diffusivity
 
c0 = 1;%M0 = 8.7e3;

syms x y t c_plot(x,y);

eta = y/sqrt(D*t);
eps = (x-d_GB/2)/sqrt(D*t);
beta = (D_GB/D)*d_GB/2/sqrt(D*t);

%%%Fisher solution
c = c0*exp(-pi^0.25*eta/beta^0.5)*(erfc(0.5*eps));

t_max = 1e6;
x_max = 6;%(solve(subs(eps,t,t_max) - 6));
y_max = 150;

c_t_gb = int(int(c,x,d_GB/2,x_max),y,0,y_max);
c_t = 2*( c_t_gb + x_max*2*c0*sqrt(D*t/pi));

fplot(t,c_t,[0,t_max],'b','LineWidth',2.0);
 
xlabel("Time (hours)");
ylabel("Mass loss (mg/cm^2)");

c_y = vpa(subs(c,x,d_GB/2),3);
% hold on;
% for i = 5e-4:1e-3:1e-1
%     temp = solve(subs(subs(c_y,t,3.6e6),d_GB,i)-0.1);
%     plot(i,temp,'*');
% end


% cot_phi = (eta*beta/4)^0.3333 - (1/3)*(1/2/eta/beta)^0.3333;

% phi = double(acot(subs(cot_phi,t,3.6e6) )*180/pi)
% hold on;
% for f = 0:4
%     c_t_temp = subs(c_t,D_GB,D*10^f);
%     plot(10^f,subs(c_t_temp,t,3.6e6),'r*','MarkerSize',12);
% end
% 
% xlabel("f (D_{GB} = 10^f D) ");
% ylabel("dC");
% set(gca,'FontSize',24);
% set(gca,'yscale','log');
% set(gca,'xscale','log');




% c_fit = 2*(d_GB/1)*(D_GB/D)^0.75*(D*t)^0.75; %0.8573 depends on D_GB

% fplot((t),(c_t),[0,3.6e6],'r');hold on;
% fplot((t),(c_fit),[0,3.6e6],'k--');

% latex(vpa(simplify(c_t),2))

% fsurf(c_plot,[0 4 0 50]);
% xlabel('x');
% ylabel('y');

% log_dM = []
% log_t = []
% hold on;
% for t_val = 1e4:1e5:1e6
%     c_plot(x,y) = subs(c,t,t_val);
%     x_int = 2*vpaintegral(c_plot,0,x_max);
%     dM = vpaintegral(x_int,0,y_max);
%     log_dM = [log_dM,log(dM)];
%     log_t = [log_t,log(t_val)];
% end
% 
% plot(log_t,log_dM);

%%%Whipple asymptotic 

