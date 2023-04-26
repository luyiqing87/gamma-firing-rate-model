% The code generates the analog of r,s,v time courses for E and I
% populations of an E-I spiking network

%% run the spiking network simulation
I_E = 0.4;
I_I = -1;
[spktime,spkcell,tvec,vvec,svec] = ...
    EI_network(0.004, 0.016, 0.002, 0.004, I_E, I_I);

%% calculate the average firing rate of each population
T_all = tvec(end);
[~,indEX] = find(spkcell<=800);
[~,indIN] = find(spkcell>800);
freqEX = length(indEX)*1000/(T_all*800);
freqIN = length(indIN)*1000/(T_all*200);
spktimeEX = spktime(indEX);
spktimeIN = spktime(indIN);
disp(['mean E rate = ',num2str(freqEX),' Hz'])
disp(['mean I rate = ',num2str(freqIN),' Hz'])

%% calculate the instantaneous firing rate
win_size = 1;
dt_slide = 0.1;
win_num = floor((T_all-win_size)/dt_slide);
freqEX_t = zeros(1,win_num)*NaN;
freqIN_t = zeros(1,win_num)*NaN;
tvec2 = zeros(1,win_num)*NaN;

for jj = 1:win_num
    t1 = dt_slide * (jj-1);
    t2 = dt_slide * (jj-1) + win_size;
    EX_inwin = find(spktimeEX>=t1 & spktimeEX<t2);
    IN_inwin = find(spktimeIN>=t1 & spktimeIN<t2);
    freqEX_t(jj) = length(EX_inwin)*1000/(win_size*800);
    freqIN_t(jj) = length(IN_inwin)*1000/(win_size*200);
    tvec2(jj) = (t1+t2)/2;
end

%% plot the r,s,v time courses of E and I populations
% colors
red1 = '#bc3333';
red2 = '#ca6f6f';
red3 = '#d8abac';
blue3 = '#a0a7cd';
blue2 = '#5a68b1';
blue1 = '#142896'; 
sE = svec(1,:);
sI = svec(2,:);
vE_u = vvec(1,:);
vI_u = vvec(2,:);
% normalize
vE = 1+1.1*vE_u/75;
vI = 1+1.1*vI_u/75;
rE = freqEX_t/1000;
rI = freqIN_t/1000;

figure
set(gcf,'unit','centimeters','position',[0,0,18,12])
pos1 = [0.15 0.6 0.75 0.35];
pos2 = [0.15 0.2 0.75 0.35];
wd = 3;
subplot('Position', pos1)
hold on
p1 = plot(tvec2,rE,'LineWidth',wd,'Color',red1);
p2 = plot(tvec,sE,'-.','LineWidth',wd,'Color',red2);
p3 = plot(tvec,vE,':','LineWidth',wd,'Color',red3);
set(gca,'Fontsize',16)
set(gca,'LineWidth',2)
xlim([350,450])
xticks(350:50:450)
xticklabels({'0','50','100'})
ylim([-0.1,1.1])
box on
ylabel('EX')
xticklabels([])
legend([p1 p2 p3],{'r_{E}','s_{E}','v_{E}'},'Location','northeast','NumColumns',3)

subplot('Position', pos2)
hold on
p4 = plot(tvec2,rI,'LineWidth',wd,'Color',blue1);
p5 = plot(tvec,sI,'-.','LineWidth',wd,'Color',blue2);
p6 = plot(tvec,vI,':','LineWidth',wd,'Color',blue3);
set(gca,'Fontsize',16)
set(gca,'LineWidth',2)
xlim([350,450])
xticks(350:50:450)
xticklabels({'0','50','100'})
ylim([-0.1,1.1])
box on
xlabel('t (ms)')
ylabel('IN')
legend([p4 p5 p6],{'r_{I}','s_{I}','v_{I}'},'Location','northeast','NumColumns',3)