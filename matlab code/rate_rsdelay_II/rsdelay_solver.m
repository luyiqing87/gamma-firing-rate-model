function [t_vec, y_mat] = rsdelay_solver(II)
% The function solves the r-s delay model 
% Example: [t_vec, y_mat] = rsdelay_solver(2);
% Input: II -- external drive to the inhibitory population
% Output: t_vec -- row vector, discretization in time
%         y_mat -- 2*length(t_vec), time courses of r,s in row 1,2

% load parameters
p = load("parameters_rsdelay.mat");

% set up dde solver
lags = p.delta;
T_all = 100;
tspan = [0,T_all];
func = @(t,x,z) ddefun(t,x,z,p);

% solve the dde
sol = dde23(func,lags,@history,tspan);
t_vec = linspace(0,T_all,10000)';
y_mat = deval(sol,t_vec);

% plot the time course    
blue1 = '#142896'; 
blue2 = '#5a68b1';
figure
hold on
plot(t_vec,y_mat(1,:),'Color',blue1,'LineStyle','-','LineWidth',2)
plot(t_vec,y_mat(2,:),'Color',blue2,'LineStyle','--','LineWidth',2)
legend('r', 's')
xlabel('t (ms)')
title(['I = ',num2str(II)])
box on
set(gca,'LineWidth',2)
set(gca,'Fontsize',20)
ylim([-0.1,1.1])
xlim([0,100])
set(gcf,'unit','centimeters','position',[0,10,15,10])
  
% main function for r-s delay model
function dydt = ddefun(t,x,z,p)
    dydt = zeros(2,1);
    dydt(1) = (-x(1)+F(-p.w*x(2)+II,p))/p.taur;
    dydt(2) = (-x(2)+p.gamma*z(1)*(1-x(2))+p.s0)/p.taus;
end

% other functions
function y=history(t)
    y(1) = 0;
    y(2) = 0;
end

function y = F(x,p)
    y = 1./(1+exp((p.thF-x)/p.kF));
end

end