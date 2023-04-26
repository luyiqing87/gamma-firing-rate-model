% The code plots the mean-field r,s,v time courses of Wang-Buszaki network

%% run the simulation
Idrive = 2;
Isigma = 0.02;
[spktime, spkcell, tvec, vvec, svec] = WB_II(Idrive, Isigma, 0);

%% calculate instantaneous firing rate
N = 200;
t_all = tvec(end);
window_size = 1;
window_slide = 0.05;
window_num=floor((t_all-window_size)/window_slide);
t_inst = zeros(1,window_num);
f_inst = zeros(1,window_num);
for j = 1:window_num
    bd1 = (j-1)*window_slide;
    bd2 = j*window_slide+window_size;
    t_inst(j) = (bd1+bd2)/2;
    f_inst(j) = length(find(spktime>=bd1 & spktime<bd2))...
        *(1000/window_size)/N;
end

%% plot the analog of r, s, v time courses
% colors
blue1='#142896'; 
blue2='#5a68b1';
blue3='#a0a7cd';
% generate plot
figure
hold on
vscaled_t = ((mean(vvec)+75)*1.1/75)-0.1;
rscaled_t = f_inst/max(f_inst);
plot(t_inst,rscaled_t,'Color',blue1,'LineStyle','-','LineWidth',2)
plot(tvec,mean(svec),'Color',blue2,'LineStyle','--','LineWidth',2)
plot(tvec,vscaled_t,'Color',blue3,'LineStyle','-.','LineWidth',2)
legend('r', 's', 'v')
ylim([-0.1, 1.1])
xlim([t_all-100, t_all])
xlabel('t (ms)')
set(gca,'Fontsize',22)
set(gca,'LineWidth',2)
box on
title(['I_\mu = ', num2str(Idrive), ', I_\sigma = ', num2str(Isigma)])
set(gcf,'unit','normalized','position',[0,0.1,0.3,0.3])

