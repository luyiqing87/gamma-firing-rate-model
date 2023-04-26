function [spktime,spkcell,tvec,vvec_t,svec_t] = ...
    EI_network(gee,gii,gei,gie,I_exc,I_inh,x_ini)
% Examples: 
% (1) PING (Fig 1A): 
% [spktime,spkcell,tvec,svec_t,vvec_t] = EI_network(0.004,0.016,0.002,0.004,0.4,-1);
% (2) ING (Fig 1B):
% [spktime,spkcell,tvec,svec_t,vvec_t] = EI_network(0.004,0.016,0.002,0.004,0.4,1);
% Input: gee, gii, gei, gie -- maximal synaptic conductances
%        I_exc -- mean applied current to each E cell
%        I_inh -- mean applied current to each I cell
%        x_ini -- 3*N matrix, optional initial conditions 
% Output: spkcell -- row vector, the indices of cells that fire
%         spktime -- row vector, the corresponding spike time
%         tvec -- row vector, discretization of time in simulation 
%         vvec_t -- 2 * length(tvec) matrix, each row stores the time 
%         course of the mean membrane potential of E or I cells
%         svec_t -- 2 * length(tvec) matrix, each row stores the time 
%         course of the mean synaptic gating variable of E or I cells

%% plot settings
plot_rasters = 1; % set to 1 to generate a raster plot
plot_gsyn = 1; % set to 1 to plot the time courses of synaptic conductances
plot_Isyn = 1; % set to 1 to plot the time courses of synaptic currents
red = '#bc3333';
blue = '#142896';

%% set parameters
% network parameters
N = 1000; % total number of cells
num_ex = N*0.8; % number of E cells
num_in = N*0.2; % number of I cells
ConMat = ConEI(num_ex,num_in,0.05,0.3,0.3,0.3);
W = [gee*ones(num_ex),gei*ones(num_ex,num_in);...
      gie*ones(num_in,num_ex),gii*ones(num_in)];     
ConMatW = ConMat.*W; % connectivity matrix with weight

% discrete time steps
t_all = 500; % total simulation time
dt = 0.05; % time step for numerical method
totalpts = t_all/dt;
tvec = dt:dt:(dt*totalpts);

% applied currents
Iapp_ex = I_exc - abs(I_exc)*0.05 + abs(I_exc)*0.1*rand(1,num_ex); 
Iapp_in = I_inh - abs(I_inh)*0.05 + abs(I_inh)*0.1*rand(1,num_in);
Iapp = [Iapp_ex, Iapp_in];

% parameters for single neuron model
c = 1;
gna = 24;
gkdr = 3;
gl = 0.02;
ena = 55;
ek = -90;
el = -60;

% parameters for synapses
alpha = [ones(1,num_ex)*14/3,ones(1,num_in)*53/11];
beta = [ones(1,num_ex)*1/3,ones(1,num_in)*2/11];
tq = [ones(1,num_ex)*0.5,ones(1,num_in)*0.6];
vsyn_ex = 0;
vsyn_in = -75;

%% Initialization
if nargin == 6 % if no initial condition is given, generate it randomly
    v = -70+90*rand(1,N);
    h = 0.2+0.6*rand(1,N);
    n = 0.2+0.6*rand(1,N);
elseif nargin == 7
    v = x_ini(1,:);
    h = x_ini(2,:);
    n = x_ini(3,:);
end
s = zeros(1,N);
ind = 1;
spktime_size_pre = (t_all*N*250)/1000;
spktime = zeros(1,spktime_size_pre);
spkcell = zeros(1,spktime_size_pre);
spk = zeros(1,N);
vt2 = v;

Isyn_t = zeros(4,totalpts);
gsyn_t = zeros(4,totalpts);
vvec_t = zeros(2,totalpts);
svec_t = zeros(2,totalpts);

%% solve by RK4
for i = 1:totalpts
    t = dt*(i-1);
    vt1 = vt2;
    vt2 = v;
    [vk1,Isyn_vec1,gsyn_vec1] = dvdt(v,h,n,s);
    hk1 = dhdt(v,h);
    nk1 = dndt(v,n);
    sk1 = dsdt(spk,s,t);
    
    [vk2,Isyn_vec2,gsyn_vec2] = dvdt(v+vk1*dt/2,h+hk1*dt/2,...
        n+nk1*dt/2,s+sk1*dt/2);
    hk2 = dhdt(v+vk1*dt/2,h+hk1*dt/2);
    nk2 = dndt(v+vk1*dt/2,n+nk1*dt/2);
    sk2 = dsdt(spk,s+sk1*dt/2,t+dt/2);
    
    [vk3,Isyn_vec3,gsyn_vec3] = dvdt(v+vk2*dt/2,h+hk2*dt/2,...
        n+nk2*dt/2,s+sk2*dt/2);
    hk3 = dhdt(v+vk2*dt/2,h+hk2*dt/2);
    nk3 = dndt(v+vk2*dt/2,n+nk2*dt/2);
    sk3 = dsdt(spk,s+sk2*dt/2,t+dt/2);
     
    [vk4,Isyn_vec4,gsyn_vec4] = dvdt(v+vk3*dt,h+hk3*dt,n+nk3*dt,s+sk3*dt);
    hk4 = dhdt(v+vk3*dt,h+hk3*dt);
    nk4 = dndt(v+vk3*dt,n+nk3*dt);
    sk4 = dsdt(spk,s+sk3*dt,t+dt);
    
    v = v+(dt/6)*(vk1+2*vk2+2*vk3+vk4);
    h = h+(dt/6)*(hk1+2*hk2+2*hk3+hk4);
    n = n+(dt/6)*(nk1+2*nk2+2*nk3+nk4);
    s = s+(dt/6)*(sk1+2*sk2+2*sk3+sk4);
    Isyn_vec = (Isyn_vec1+2*Isyn_vec2+2*Isyn_vec3+Isyn_vec4)/6;
    gsyn_vec = (gsyn_vec1+2*gsyn_vec2+2*gsyn_vec3+gsyn_vec4)/6;

    vt3=v;
    Isyn_t(:,i) = Isyn_vec;
    gsyn_t(:,i) = gsyn_vec;
    vvec_t(:,i) = [mean(v(1:800));mean(v(801:1000))];
    svec_t(:,i) = [mean(s(1:800));mean(s(801:1000))];
    
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

spktime(spktime==0)=[];
spkcell(spkcell==0)=[];

%% generate plots
[~, Index_ex] = sort(Iapp_ex,'descend');
if plot_rasters == 1
    figure(1)
    for s = 1:num_ex
        loc = spkcell==Index_ex(s);
        ti = spktime(loc);
        plot(ti,s*ones(1,length(ti)),'.','MarkerSize',6,'Color',red);
        hold on
    end
    for s = (num_ex+1):N
        loc = spkcell==s;
        ti = spktime(loc);
        plot(ti,s*ones(1,length(ti)),'.','MarkerSize',6,'Color',blue);
        hold on
    end
    hold off
    xlabel('t (ms)')
    ylabel('Cell Index')
    xlim([t_all-100, t_all])
    xticks([400,450,500])
    yticks([0,400,800,1000])
    ylim([0,1000])
    title(['I_I=',num2str(I_inh)])
    set(gca,'FontSize',18)
    set(gca,'LineWidth',2)
    box on
    set(gcf,'unit','centimeters','position',[0,0,18,10])
end

if plot_Isyn == 1
    figure(2)
    plot(tvec,-Isyn_t(1,:),'LineWidth',2,'Color',red,'LineStyle','-')
    hold on
    plot(tvec,-Isyn_t(2,:),'LineWidth',2,'Color',blue,'LineStyle','-')
    hold on
    plot(tvec,-Isyn_t(3,:),'LineWidth',2,'Color',red,'LineStyle','-.')
    hold on
    plot(tvec,-Isyn_t(4,:),'LineWidth',2,'Color',blue,'LineStyle','-.')
    set(gca,'Fontsize',18)
    set(gca,'LineWidth',2)
    box on
    xlabel('t (ms)')
    ylabel('Isyn')
    legend('from E to E','from I to I', 'from E to I','from I to E')
    set(gcf,'unit','centimeters','position',[0,10,18,10])
    xlim([t_all-100,t_all])
    xticks([400,450,500])
    ylim([-50,20])

end

if plot_gsyn == 1
    figure(3)
    plot(tvec,gsyn_t(1,:),'LineWidth',2,'Color',red,'LineStyle','-')
    hold on
    plot(tvec,gsyn_t(2,:),'LineWidth',2,'Color',blue,'LineStyle','-')
    hold on
    plot(tvec,gsyn_t(3,:),'LineWidth',2,'Color',red,'LineStyle','-.')
    hold on
    plot(tvec,gsyn_t(4,:),'LineWidth',2,'Color',blue,'LineStyle','-.')    
    ylim([-0.05,0.85])
    set(gca,'Fontsize',18)
    set(gca,'LineWidth',2)
    box on
    xlabel('t (ms)')
    ylabel('gsyn')
    legend('from E to E','from I to I', 'from E to I','from I to E')
    set(gcf,'unit','centimeters','position',[0,10,18,10])
    xlim([t_all-100,t_all])
    xticks([400,450,500])
end

%% functions 
% functions for single neuron model
    function h_inf = hinf(x)
        h_inf = 1./(1+exp((x+53)./7));
    end        
    
    function tau_h = tauh(x)
	    tau_h = 0.37+2.78./(1+exp((x+40.5)./6));
    end
    
    function dh = dhdt(v,h)
        dh = (hinf(v)-h)./(tauh(v));
    end
    
    function n_inf = ninf(x)
	    n_inf = 1./(1+exp((-x-30)./10));
    end
    
    function tau_n = taun(x)
	    tau_n = 0.37+1.85./(1+exp((x+27)./15));
    end
    
    function dn = dndt(v,n)
        dn=(ninf(v)-n)./(taun(v));
    end

    function m_inf = minf(x)
	    m_inf = 1./(1+exp((-x-30)./9.5));
    end

% functions for synaptic current
    function ds = dsdt(spk,s,t)
        q = double((t-spk)<=tq);
	    ds = q.*alpha.*(1-s)-beta.*s;
    end
    
    function [I_syn,Isyn_vec,gsyn_vec] = Isyn(v,s)
        vE = v(1:num_ex);
        vI = v(num_ex+1:N);
        sE = s(1:num_ex);
        sI = s(num_ex+1:N);
        mat_WEE = ConMatW(1:num_ex,1:num_ex);
        mat_WEI = ConMatW(1:num_ex,num_ex+1:end);
        mat_WIE = ConMatW(num_ex+1:end,1:num_ex);
        mat_WII = ConMatW(num_ex+1:end,num_ex+1:end);
        gsyn_EE = sE*mat_WEE;
        gsyn_EI = sE*mat_WEI;
        gsyn_IE = sI*mat_WIE;
        gsyn_II = sI*mat_WII;
        Isyn_EE = gsyn_EE.*(vE-vsyn_ex);
        Isyn_EI = gsyn_EI.*(vI-vsyn_ex);
        Isyn_IE = gsyn_IE.*(vE-vsyn_in);
        Isyn_II = gsyn_II.*(vI-vsyn_in);
        Isyn_E = Isyn_EE+Isyn_IE;
        Isyn_I = Isyn_EI+Isyn_II;
        I_syn = [Isyn_E,Isyn_I];   
        Isyn_vec = [mean(Isyn_EE);mean(Isyn_II);mean(Isyn_EI);mean(Isyn_IE)];
        gsyn_vec = [mean(gsyn_EE);mean(gsyn_II);mean(gsyn_EI);mean(gsyn_IE)];
    end

% main ODEs
    function [dv,Isyn_vec,gsyn_vec] = dvdt(v,h,n,s)
        [I_syn,Isyn_vec,gsyn_vec] = Isyn(v,s);
        dv = (Iapp-I_syn-gna*(minf(v).^3).*h.*(v-ena)-...
              gkdr*(n.^4).*(v-ek)-gl*(v-el))/c;    
    end

end