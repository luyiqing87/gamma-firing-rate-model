function [t_vec, y_mat] = rsv_solver(II)
% The function solves the r-s-v model by calling the ODE solver ode45
% Example: [t_vec,y_mat] = rsv_solver(1);
% Input: II -- external drive to the inhibitory population
% Output: t_vec -- column vector, discretization in time
%         y_mat -- length(t_vec)*3, time courses of r,s,v in column 1,2,3

%% set parameters
p = load("rsv_parameters.mat"); %load parameters into a structure p
y0 = [0.5 0.1 0.3]; % initial conditions r0, s0, v0
tspan = [0 1100]; % simulation time

%% solve the r,s,v system
ODE = @(t, y) eq_rsv(t, y, II, p);
[t_vec, y_mat] = ode45(ODE, tspan, y0); % solve the ODE

%% plot the time courses
blue1 = '#142896'; 
blue2 = '#5a68b1';
blue3 = '#a0a7cd';

figure
hold on
plot(t_vec,y_mat(:,1),'Color',blue1,'LineStyle','-','LineWidth',2)
plot(t_vec,y_mat(:,2),'Color',blue2,'LineStyle','--','LineWidth',2)
plot(t_vec,y_mat(:,3),'Color',blue3,'LineStyle','-.','LineWidth',2)
legend('r','s','v')
xlabel('t (ms)')
box on
set(gca,'LineWidth',2)
set(gca,'Fontsize',20)
ylim([-0.1,1.1])
xlim([1000,1100])
xticks(1000:25:1100)
xticklabels({'0','25','50','75','100'})
title(['I = ',num2str(II)])
set(gcf,'unit','centimeters','position',[0,10,15,10])

%% functions
    function dydt = eq_rsv(t, y, Ivar,p)
        dydt = zeros(3,1);
        dydt(1) = (-y(1)+f(y(3),p))/p.taur;
        dydt(2) = (-y(2)+p.gamma*q(y(1),p)*(1-y(2))+p.s0)/p.taus;
        dydt(3) = (-y(3)+p.g*y(2)*(p.vibar-y(3))+Ivar)/p.tauv;
    end

    function y = f(x,p)
       y = 1 / (1 + exp((p.thf - x)/p.kf));
    end
    
    function y = q(x,p)
       y = 1 / (1 + exp((p.thq - x)/p.kq));
    end

end