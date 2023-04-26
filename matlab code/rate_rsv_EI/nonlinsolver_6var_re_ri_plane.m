% Plot the rE and rI-'nullclines' in the projected rE-rI plane
% load parameters
p = load("rsv_EI_INGpara.mat");

% calculate nullclines and fixed point
[rE_ncl1,rI_ncl1] = rEnlc_solver(p.IE,p);
[rE_ncl2,rI_ncl2] = rInlc_solver(p.II,p);

% plot the projected rE-rI plane
figure
hold on
rEnull = plot(rE_ncl1,rI_ncl1,'Color','#bc3333','LineWidth',2);
rInull = plot(rE_ncl2,rI_ncl2,'Color','#142896','LineWidth',2);
legend([rEnull rInull],'dr_E/dt=0','dr_I/dt=0')
ylim([-0.1,1.1])
xlim([-0.1,1.1])
title(['I_{exc}=',num2str(IE),', I_{inh}=',num2str(II)])
set(gca,'Fontsize',18)
set(gca,'LineWidth',2)
xlabel('r_E')
ylabel('r_I')
axis square
box on
set(gcf,'unit','centimeters','position',[0,10,12,12])

%% functions
function v = vE(r, p)
    v = p.kfE * log(r./(1-r))+p.vthE;
end

function v = vI(r, p)
    v = p.kfI * log(r./(1-r))+p.vthI;
end

function q_out = q(r, p)
    q_out = (1 + exp(-(r-p.rth)./p.kq)).^(-1);
end

function s = sE(r, p)
    s = 1 - (1-p.sE0)./(p.gammaE*q(r, p)+1);
end

function s = sI(r, p)
    s = 1 - (1-p.sI0)./(p.gammaI*q(r, p)+1);
end

function out = dvEdt(rE,rI,IE,p)
    out = -vE(rE,p) + p.gEE*sE(rE, p).*(p.vEbar-vE(rE,p)) ...
        + p.gEI*sI(rI,p).*(p.vIbar-vE(rE,p)) + IE;
end

function out = dvIdt(rE,rI,II,p)
    out = -vI(rI,p) + p.gIE*sE(rE,p).*(p.vEbar-vI(rI,p)) ...
        + p.gII*sI(rI,p).*(p.vIbar-vI(rI,p)) + II;
end

function out = dvdt(r,IE,II,p)
    rE = r(1);
    rI = r(2);
    out(1) = -vE(rE,p)+p.gEE*sE(rE,p).*(p.vEbar-vE(rE,p))+p.gEI*sI(rI,p).*(p.vIbar-vE(rE,p))+IE;
    out(2) = -vI(rI,p)+p.gIE*sE(rE,p).*(p.vEbar-vI(rI,p))+p.gII*sI(rI,p).*(p.vIbar-vI(rI,p))+II;
end

function [rE_ncl,rI_ncl] = rEnlc_solver(IE,p)
    n = 500;
    delta = 5e-3; 
    rE_ncl = linspace(delta,1-delta,n);
    rI_ncl = rE_ncl*0;
    x0 = 0.2;
    for i=1:n
        fun = @(rI) dvEdt(rE_ncl(i),rI,IE,p);
        rI_ncl(i) = fsolve(fun,x0,optimoptions('fsolve','Display','off'));
    end
end

function [rE_ncl,rI_ncl] = rInlc_solver(II,p)
    n = 500;
    delta = 1e-15;
    rI_ncl = linspace(delta,1-delta,n);
    rE_ncl = rI_ncl*0;
    x0 = 0.5;
    for i=1:n
        fun = @(rE) dvIdt(rE,rI_ncl(i),II,p);
        rE_ncl(i) = fsolve(fun,x0,optimoptions('fsolve','Display','off'));
    end
end



