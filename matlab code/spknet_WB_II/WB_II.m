function [spktime, spkcell, tvec, vvec, svec] = WB_II(Idrive, Isigma, save_ON)
% The code simulates an inhibitory network of N Wang-Buzsaki neurons
% Example (parameters used to generate Fig 2A):
% [spktime, spkcell, tvec, vvec, svec] = WB_II(2, 0.02, 0);
%
% Input: Idrive -- mean applied current to each cell
%        Isigma -- standard deviation of applied current to each cell
%        save_ON -- set to 1 to save the data and figure
% Output: spkcell -- row vector, the indices of cells that fire
%         spktime -- row vector, the corresponding spike time
%         tvec -- row vector, discretization of time in simulation 
%         vvec -- N * length(tvec) matrix, each row stores the time course 
%         of the membrane potential of a cell
%         svec -- N * length(tvec) matrix, each row stores the time course 
%         of the synaptic gating variable of a cell

%% set parameters
% network setup
N = 200; % total number of I cells
Msyn = 60; % the number of presynaptic cells for each I cell
ConMat = zeros(N); % connectivity matrix
for j = 1:N
    index = randperm(N);
    indexcut = index(1:Msyn+1);
    if any(indexcut==j)
        id = indexcut==j;
        indexcut(id) = [];
    else
        indexcut = indexcut(1:end-1);
    end
    for i = 1:Msyn
        ConMat(indexcut(i),j) = 1;
    end
end

% discretization in time
t_all = 500; % total simulation time
dt = 0.05; % time step for the numerical method
totalpts = floor(t_all/dt);
tvec = dt:dt:(dt*totalpts);

% single cell parameters
c = 1;
gna = 35;
gk = 9;
gl = 0.1;
ena = 55;
ek = -90;
el = -65;
phi = 5;

% network parameters
alpha = 12;
beta = 0.1;
thre = 0;
vsyn_in = -75;
k = 2;
gii = 0.1/Msyn;
Iapp = Idrive+randn(1,N)*Isigma; % applied current to each cell, row vector 

%% Initialization
% initialize membrane voltages and gating variables
% WB 1996 used v from -70 to -50, h and n at rest
% rest values: v=-64.0228, h=0.7809, n=0.0890
v = -60 + 20*randn(1,N);
h = 0.85 + 0.1*randn(1,N);
n = 0.05 + 0.1*randn(1,N);
s = 0.05 + 0.1*randn(1,N);
ind = 1;

spktime_size_pre = floor((t_all*N*250)/1000);
spktime = zeros(1, spktime_size_pre);
spkcell = zeros(1, spktime_size_pre);
spk = zeros(1,N);
Isyn_t = zeros(1,totalpts);
gsyn_t = zeros(1,totalpts);
vvec = zeros(N,totalpts);
svec = zeros(N,totalpts);
vt2 = v;

%% solve the main ODEs by RK4
for i = 1:totalpts
    t = dt*(i-1);
    vt1 = vt2;
    vt2 = v;
    [vk1, Isyn_vec1, gsyn_vec1] = dvdt(v, h, n, s);
    hk1 = dhdt(v, h);
    nk1 = dndt(v, n);
    sk1 = dsdt(v, s);
    
    [vk2, Isyn_vec2, gsyn_vec2] = dvdt(v+vk1*dt/2, h+hk1*dt/2, ...
        n+nk1*dt/2, s+sk1*dt/2);
    hk2 = dhdt(v+vk1*dt/2, h+hk1*dt/2);
    nk2 = dndt(v+vk1*dt/2, n+nk1*dt/2);
    sk2 = dsdt(v+vk1*dt/2, s+sk1*dt/2);
    
    [vk3, Isyn_vec3, gsyn_vec3] = dvdt(v+vk2*dt/2, h+hk2*dt/2, ...
        n+nk2*dt/2, s+sk2*dt/2);
    hk3 = dhdt(v+vk2*dt/2, h+hk2*dt/2);
    nk3 = dndt(v+vk2*dt/2, n+nk2*dt/2);
    sk3 = dsdt(v+vk2*dt/2, s+sk2*dt/2);
    
    [vk4, Isyn_vec4, gsyn_vec4] = dvdt(v+vk3*dt, h+hk3*dt, ...
        n+nk3*dt, s+sk3*dt);
    hk4 = dhdt(v+vk3*dt, h+hk3*dt);
    nk4 = dndt(v+vk3*dt, n+nk3*dt);
    sk4 = dsdt(v+vk3*dt, s+sk3*dt);
    
    v = v + (dt/6)*(vk1+2*vk2+2*vk3+vk4);
    h = h + (dt/6)*(hk1+2*hk2+2*hk3+hk4);
    n = n + (dt/6)*(nk1+2*nk2+2*nk3+nk4);
    s = s + (dt/6)*(sk1+2*sk2+2*sk3+sk4);

    Isyn_vec = (Isyn_vec1+2*Isyn_vec2+2*Isyn_vec3+Isyn_vec4)/6;
    gsyn_vec = (gsyn_vec1+2*gsyn_vec2+2*gsyn_vec3+gsyn_vec4)/6;
    vt3 = v;

    Isyn_t(:,i) = Isyn_vec;
    gsyn_t(:,i) = gsyn_vec;
    
    vvec(:,i) = v';
    svec(:,i) = s';

    %save the time and cell index for each spike
    for j = 1:N
        if vt2(j)>=0 && vt2(j)>=vt1(j) && vt2(j)>=vt3(j)
            spktemp = t+0.001;
            spk(j) = spktemp;
            spktime(ind) = spktemp;
            spkcell(ind) = j;
            ind = ind+1;
        end
    end
end

spktime(spktime==0) = [];
spkcell(spkcell==0) = [];

%% generate raster plot 
xlim_bd1 = t_all-100;
xlim_bd2 = t_all;
figure
for s = 1:N
    loc = (spkcell==s);
    ti = spktime(loc);
    scatter(ti, s*ones(1,length(ti)), 'k.');
    hold on
end
hold off
xlabel('t (ms)')
ylabel('Cell ID')
set(gcf,'unit','normalized','position',[0,0.1,0.3,0.3])
xlim([xlim_bd1,xlim_bd2])
set(gca,'Fontsize',22)
set(gca,'LineWidth',2)
box on
title(['I_\mu = ', num2str(Idrive), ', I_\sigma = ', num2str(Isigma)])

% save the data and raster plot
if save_ON == 1
    filename = ['I=', num2str(Idrive), '_std=', num2str(Isigma)]; 
    savefig([filename,'_rs.fig'])
    print('-dpng',[filename,'_rs.png'])
    save([filename,'.mat'],'spkcell','spktime','Isyn_t','gsyn_t',...
        'tvec','vvec','svec')
    close
end

%% functions
% functions for single cell    
    % gating variable m
    function alpha_m = alpham(v)
        alpha_m = -0.1*(v+35)./(exp(-0.1.*(v+35))-1);
    end        
    
    function beta_m = betam(v)
        beta_m = 4.*exp(-(v+60)./18);
    end 
    
    function m_infty = minf(v)
        m_infty = alpham(v)./(alpham(v)+betam(v));
    end
    
    % gating variable h
    function alpha_h = alphah(v)
        alpha_h = 0.07.*exp(-(v+58)./20);
    end        
    
    function beta_h = betah(v)
        beta_h = 1./(exp(-0.1*(v+28))+1);
    end 
    
    function dh = dhdt(v,h)
        dh = phi.*(alphah(v).*(1-h)-betah(v).*h);
    end
    
    % gating variable n
    function alpha_n = alphan(v)
        alpha_n = -0.01.*(v+34)./(exp(-0.1.*(v+34))-1);
    end        
    
    function beta_n = betan(v)
        beta_n = 0.125.*exp(-(v+44)./80);
    end 
    
    function dn = dndt(v,n)
        dn = phi.*(alphan(v).*(1-n)-betan(v).*n);
    end
    
% functions for synapses
    function qv = q(v) % sigmoid
        qv = 1./(1+exp(-(v-thre)./k));
    end
    
    function ds = dsdt(v,s)
        ds = alpha.*q(v).*(1-s)-beta.*s;
    end
    
    function [Isyn,Isyn_vec,wsyn_vec] = I_syn(v,s)
        wsyn = gii*s*ConMat; 
        Isyn = wsyn.*(v-vsyn_in); 
        Isyn_vec = mean(Isyn);
        wsyn_vec = mean(wsyn);     
    end
    
% main ODE
    function [dv, Isyn_vec, wsyn_vec] = dvdt(v, h, n, s)
        [Isyn, Isyn_vec, wsyn_vec] = I_syn(v, s);
        Ina = gna*(minf(v).^3).*h.*(v-ena);
        Ik = gk*(n.^4).*(v-ek);
        Il = gl*(v-el);
        dv = (Iapp-Isyn-Ina-Ik-Il) / c;  
    end

end