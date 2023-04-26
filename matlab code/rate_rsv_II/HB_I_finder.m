% solve for the Hopf bifurcations with respect to the input drive I_I 

% load parameters
p = load("rsv_parameters.mat"); %load parameters into a structure p

% setup nonlinear equation solver
fun = @(xx) HB_eq(xx, p);
lb = [0;0;0;0;0];
ub = [1;1;1;10;10];
options = optimoptions('lsqnonlin','Display','off','FunctionTolerance',1e-10);

% left HB
x0 = [0.1;0.1;p.thf;1;1]; 
[xHBl, resnorml] = lsqnonlin(fun, x0, lb, ub, options);
[rHBl, sHBl, vHBl, betaHBl, IHBl] = ...
    deal(xHBl(1),xHBl(2),xHBl(3),xHBl(4),xHBl(5));

% right HB
x0 = [0.9;0.9;p.thf;1;2]; 
[xHBu, resnormu] = lsqnonlin(fun, x0, lb, ub, options);
[rHBu, sHBu, vHBu, betaHBu, IHBu] = ...
    deal(xHBu(1),xHBu(2),xHBu(3),xHBu(4),xHBu(5));

% display I_HB and corresponding frequency
disp(['I_HB1 = ', num2str(IHBl),', f = ', num2str(1000*betaHBl/(2*pi))])
disp(['I_HB2 = ', num2str(IHBu),', f = ', num2str(1000*betaHBu/(2*pi))])

% nonlinear equations for Hopf Bifurcations
function y = HB_eq(xx, p)
    x = xx(1:3);
    beta = xx(4);
    I_HB = xx(5);

    alphar = 1 / p.taur;
    alphas = (1 + p.gamma * q(x(1),p)) / p.taus;
    alphav = (1 + p.g * x(2)) / p.tauv;

    A = p.gamma*p.g*fprime(x(3),p)*qprime(x(1),p)*(1-x(2))*(x(3)-p.vibar)...
        /(p.taur*p.taus*p.tauv);  
    y(1) = A - (alphar+alphas+alphav)*beta^2 + alphar*alphas*alphav;
    y(2) = -beta^2 + (alphar*alphas + alphas*alphav + alphav*alphar);
    y(3) = (-x(1) + f(x(3),p))/p.taur;
    y(4) = (-x(2) + p.gamma*q(x(1),p)*(1-x(2)) + p.s0)/p.taus;
    y(5) = (-x(3) + p.g*x(2)*(p.vibar-x(3)) + I_HB)/p.tauv;

end

% other functions
function y = f(x, p)
   y = 1 / (1+exp((p.thf-x)/p.kf));
end

function y = q(x, p)
   y = 1 / (1+exp((p.thq-x)/p.kq));
end

function y = fprime(x, p)
    A = exp( (p.thf-x) / p.kf );
    y = A / (p.kf * (1+A)^2);
end

function y = qprime(x, p)
    A = exp( (p.thq-x) / p.kq );
    y = A / (p.kq * (1+A)^2);
end
