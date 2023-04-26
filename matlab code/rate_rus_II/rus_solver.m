function [t_vec, y_mat, y_last] = rus_solver(y_ic)
% Solving the r-u-s rate model, y = [r;u;s] in column
% Example: 
% (1) Default IC: [t_vec, y_mat, y_last] = rus_solver;
% (2) IC-last: [t_vec, y_mat, y_last] = rus_solver(y_last); 
%     need to run (1) first before running (2)
%
% Input: y_ic -- the initial condition [r0;s0;v0] in column,
%                has a default value when not specified, to use the value 
%                at the end of last trial, set it to be y_last
%
% Output: t_vec -- column vector storing the time
%         y_mat -- matrix storing the r(t), u(t), s(t) in the 1st, 2nd, 
%                  3rd row
%         y_last -- 3*1 vector storing [r;s;t] at the end of simulation
%         par -- 1*12 vector storing all parameters

% load the parameters
p = load('parameters_rus.mat');

% set up the ODE solver
if nargin == 0
    y_ic = [0.5;0.1;0.3]; % default IC when not specified
end
T = 100; % total time
dt = 0.05;
tspan = [0 T];
t_vec = 0:dt:T;

% solve by MATLAB ODE solver ode45
myODE = @(t,y) eq_rsv(t,y,p);
y_sol = ode45(myODE, tspan, y_ic);
y_mat = deval(y_sol,t_vec);
y_last = y_mat(:,end);

% plotting the time course

green = [0.4660, 0.6740, 0.1880];
blue1 = '#142896'; 
blue2 = '#5a68b1';
figure
plot(t_vec,y_mat(1,:),'Color',blue1,'LineWidth',2,'LineStyle','-')
hold on
plot(t_vec,y_mat(2,:),'Color',green,'LineWidth',2,'LineStyle','-.')
plot(t_vec,y_mat(3,:),'Color',blue2,'LineWidth',2,'LineStyle','--')
legend('r','u','s')
xlabel('t (ms)')
ylim([-0.1,1.1])
box on
set(gca,'LineWidth',2)
set(gca,'Fontsize',20)
set(gcf,'unit','centimeters','position',[0,10,15,10])

% main r-u-s ODEs 
    function dydt = eq_rsv(t, y, p)
        dydt = zeros(3,1);
        dydt(1) = (-y(1)+F(p.I-p.w*y(3),p))/p.taur;
        dydt(2) = (-y(2)+y(1))/p.tauu;
        dydt(3) = (-y(3)+p.gamma*q(y(2),p)*(1-y(3))+p.s0)/p.taus;
    end
    
% other functions
    function y = F(x,p)
       y = 1 / (1+exp((p.thF-x)/p.kF));
    end
    
    function y = q(x,p)
       y = 1 / (1+exp((p.thq-x)/p.kq));
    end


end