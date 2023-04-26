function [t_vec, y_mat] = rsv_EI_solver(IE,II)
% The function solves the 6-var r-s-v model for E-I network
% Example: [t_vec,y_mat] = rsv_EI_solver(-0.2,0);
% Input: IE -- external drive to the excitatory population
%        II -- external drive to the inhibitory population
% Output: t_vec -- column vector, discretization in time
%         y_mat -- length(t_vec)*6, time courses of rE,sE,vE,rI,sI,vI 
%                  in columns 

%% set parameters
p = load("rsv_EI_INGpara.mat"); %load parameters into a structure p
y0 = [0.1 0.1 0.1 0.8 0.8 0.5]; % initial conditions r0, s0, v0
tspan = [0 200]; % simulation time

%% solve the r,s,v system
ODE = @(t, y) eq_rsv(t, y, IE, II, p);
[t_vec, y_mat] = ode45(ODE, tspan, y0); % solve the ODE

%% plot the time courses
red1 = '#bc3333';
red2 = '#ca6f6f';
red3 = '#d8abac';
blue1 = '#142896'; 
blue2 = '#5a68b1';
blue3 = '#a0a7cd';

figure
subplot(2,1,1)
title(['I_E = ',num2str(IE),', I_I = ',num2str(II)])
hold on
plot(t_vec,y_mat(:,1),'Color',red1,'LineStyle','-','LineWidth',2)
plot(t_vec,y_mat(:,2),'Color',red2,'LineStyle','--','LineWidth',2)
plot(t_vec,y_mat(:,3),'Color',red3,'LineStyle','-.','LineWidth',2)
legend('r_E','s_E','v_E')
box on
set(gca,'LineWidth',2)
set(gca,'Fontsize',20)
ylim([-0.1,1.1])
xticklabels([])
subplot(2,1,2)
hold on
plot(t_vec,y_mat(:,4),'Color',blue1,'LineStyle','-','LineWidth',2)
plot(t_vec,y_mat(:,5),'Color',blue2,'LineStyle','--','LineWidth',2)
plot(t_vec,y_mat(:,6),'Color',blue3,'LineStyle','-.','LineWidth',2)
legend('r_I','s_I','v_I')
xlabel('t (ms)')
box on
set(gca,'LineWidth',2)
set(gca,'Fontsize',20)
ylim([-0.1,1.1])
set(gcf,'unit','centimeters','position',[0,10,15,18])

%% functions
    function dydt = eq_rsv(t, y, IE, II, p)
        dydt = zeros(6,1);
        dydt(1) = (-y(1)+fE(y(3),p))/p.taurE;
        dydt(2) = (-y(2)+p.gammaE*q(y(1),p)*(1-y(2))+p.sE0)/p.tausE;
        dydt(3) = (-y(3)+p.gEE*y(2)*(p.vEbar-y(3))+p.gEI*y(5)*(p.vIbar-y(3))+IE)/p.tauvE;
        dydt(4) = (-y(4)+fI(y(6),p))/p.taurI;
        dydt(5) = (-y(5)+p.gammaI*q(y(4),p)*(1-y(5))+p.sI0)/p.tausI;
        dydt(6) = (-y(6)+p.gIE*y(2)*(p.vEbar-y(6))+p.gII*y(5)*(p.vIbar-y(6))+II)/p.tauvI;
    end

    function y = fE(x, p)
       y = 1 / (1 + exp((p.vthE - x)/p.kfE));
    end
    
    function y = fI(x, p)
       y = 1 / (1 + exp((p.vthI - x)/p.kfI));
    end

    function y = q(x, p)
       y = 1 / (1 + exp((p.rth - x)/p.kq));
    end

end