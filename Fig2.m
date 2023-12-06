%% Fig2

%% stimulus phase coherence (0.125s)

clear
clc

fs = 100;

% experiement 1--- 125ms

load(".\data_availability\Stimuli\125_exp1_env.mat");
load(".\data_availability\Stimuli\anechoic_env.mat");

sig_env_echo = [];
sig_env = [];
for i = 1:size(sig_speech_echo_sum,2)
    sig_env_echo = [sig_env_echo;sig_speech_echo_sum{i}]; % echoic speech after depression
    sig_env = [sig_env;sig_speech_sum{i}];% original anechoic speech
end

sig1 = mean(sig_env,2);
sig2 = mean(sig_env_echo,2);

for ind = 1:1000
    try  x = sig1(fs*2*(ind-1)+1:fs*2*ind);
    catch,disp([num2str(ind)]);break;end
    sig_depn(:,ind) = x;
    try  y = sig2(fs*2*(ind-1)+1:fs*2*ind);
    catch,break;end
    sig_envn(:,ind) = y;
end
fre_dep = fft(sig_depn);
fre_env = fft(sig_envn);
phy = angle(fre_dep)-angle(fre_env);
phy_coh = mean(cos(phy),2).^2+mean(sin(phy),2).^2;
figure
f = 1:size(phy_coh,1);f=f-1;f=f/size(phy_coh,1)*fs;
p1 = plot(f,mean(phy_coh,2),'color','k','linewidth',2)
hold on
plot([4 4],[0 1],'color',[176 31 36]./255,'linestyle','--','linewidth',2);

xlim([0.8 14])
ylim([0 1])
xlabel('Frequency(Hz)')
ylabel('Phase coherence')
title('0.125-s echoic speech','FontSize',15,'FontWeight','bold')
set(gcf,'unit','centimeters','position',[1.50 1.5 15 10])
set(gca,'fontname','calibri');
set(gca,'FontSize',14,'Linewidth',0.8,'GridAlpha',.15)

% experiement 1--- 250ms
clear
clc
load(".\data_availability\Stimuli\250_exp1_env.mat");
load(".\data_availability\Stimuli\anechoic_env.mat");

sig_env_echo = [];
sig_env = [];
for i = 1:size(sig_speech_echo_sum2,2)
    sig_env_echo = [sig_env_echo;sig_speech_echo_sum2{i}]; % echoic speech after depression
    sig_env = [sig_env;sig_speech_sum{i}];% original anechoic speech
end

sig1 = mean(sig_env,2);
sig2 = mean(sig_env_echo,2);

for ind = 1:1000
    try  x = sig1(fs*2*(ind-1)+1:fs*2*ind);
    catch,disp([num2str(ind)]);break;end
    sig_depn(:,ind) = x;
    try  y = sig2(fs*2*(ind-1)+1:fs*2*ind);
    catch,break;end
    sig_envn(:,ind) = y;
end
fre_dep = fft(sig_depn);
fre_env = fft(sig_envn);
phy = angle(fre_dep)-angle(fre_env);
phy_coh = mean(cos(phy),2).^2+mean(sin(phy),2).^2;
figure
f = 1:size(phy_coh,1);f=f-1;f=f/size(phy_coh,1)*fs;
p1 = plot(f,mean(phy_coh,2),'color','k','linewidth',2)
hold on
plot([2 2],[0 1],'color',[176 31 36]./255,'linestyle','--','linewidth',2);
plot([6 6],[0 1],'color',[176 31 36]./255,'linestyle','--','linewidth',2);

xlim([0.8 14])
ylim([0 1])
xlabel('Frequency(Hz)')
ylabel('Phase coherence')
title('0.25-s echoic speech','FontSize',15,'FontWeight','bold')
% set(gca,'fontsize',15);
set(gcf,'unit','centimeters','position',[1.50 1.5 15 10])
set(gca,'fontname','calibri');
set(gca,'FontSize',14,'Linewidth',0.8,'GridAlpha',.15)

%% phase coherence (exp 1)
clear
clc
load('.\data_availability\phase_coherence\Exp1\phy_coh_clean.mat');
color_list = {'#0072BD','#D95319','#EDB120'};
fs = 100;
cili = 1:3:306;
tidu = setdiff(1:306,cili);

f = 1:size(phy_coh_clean,1);f=f-1;f=f/size(phy_coh_clean,1)*fs;
sub = [1:15];
condition = {'Anechoic speech','0.125-s echoic speech','0.25-s echoic speech'};
figure
for con = 1:2
%     figure
    pc_2 = squeeze(mean(phy_coh_clean(:,tidu,con,sub),2));
    err_2 = std(pc_2,0,2)./sqrt(length(sub));
    f2 = fill([f,fliplr(f)],[mean(pc_2,2)+err_2;flipud(mean(pc_2,2)-err_2)],[195 195 195]./255);
    f2.EdgeColor = 'none';
    f2.FaceColor = color_list{con};
    f2.FaceAlpha = 0.3; 
    hold on
    p2 = plot(f,mean(mean(phy_coh_clean(:,tidu,con,sub),2),4),'linewidth',2,'color',color_list{con});

    xlim([0.8 14])
    ylabel('phase coherence')
    ylim([0.002 0.012])
    title('Experiment 1 (0.125)','FontSize',15,'FontWeight','bold')
    xlabel('frequency (Hz)')
    set(gca,'fontname','calibri');
    set(gcf,'unit','centimeters','position',[4 4 15 10])
    set(gca,'FontSize',14,'Linewidth',0.8,'GridAlpha',.15)

    p3 = plot([4 4],[0 0.015],'color',color_list{2},'linewidth',2,'LineStyle','--');
load('.\data_availability\phase_coherence\Exp1\phy_coh_shuffle.mat');
phy_coh_shuffle = phy_sum;
p4 = plot(f,squeeze(mean(mean(phy_coh_shuffle(:,:,:),3),2)),'color','k','LineWidth',1.5,'LineStyle',':');

for i = 1:21
    a = mean(mean(phy_coh_clean(i,tidu,con,sub),2),4);
    b = squeeze(mean(phy_coh_shuffle(i,:,:),3));
    M_value(i) = (length(find(a<b))+1)/(size(phy_coh_shuffle,2)+1);
end
% 
[h,q,p_new] = fdr_bh(M_value);
for i = 1:21
    p_value = p_new(i);
    clear a b
    if p_value<0.05
        f4 = fill([[f(i)-0.25 f(i)+0.25],[f(i)+0.25 f(i)-0.25]],...
            [0.0102+0.0005*con+0.00015,0.0102+0.0005*con+0.00015,0.0102+0.0005*con-0.00015,0.0102+0.0005*con-0.00015],'k');
        f4.EdgeColor = 'none';
        f4.FaceColor = color_list{con};
        f4.FaceAlpha = 0.3; 
    end
end
end

figure
for con = [1 3]
%     figure
    pc_2 = squeeze(mean(phy_coh_clean(:,tidu,con,sub),2));
    err_2 = std(pc_2,0,2)./sqrt(length(sub));
    f2 = fill([f,fliplr(f)],[mean(pc_2,2)+err_2;flipud(mean(pc_2,2)-err_2)],[195 195 195]./255);
    f2.EdgeColor = 'none';
    f2.FaceColor = color_list{con};
    f2.FaceAlpha = 0.3; 
    hold on
    p2 = plot(f,mean(mean(phy_coh_clean(:,tidu,con,sub),2),4),'linewidth',2,'color',color_list{con});

    xlim([0.8 14])
    ylabel('phase coherence')
    ylim([0.002 0.012])
    title('Experiment 1 (0.25)','FontSize',15,'FontWeight','bold')
    xlabel('frequency (Hz)')
    set(gca,'fontname','calibri');
    set(gcf,'unit','centimeters','position',[4 4 15 10])
    set(gca,'FontSize',14,'Linewidth',0.8,'GridAlpha',.15)

    p3 = plot([2 2],[0 0.015],'color',color_list{3},'linewidth',2,'LineStyle','--');
    p3 = plot([6 6],[0 0.015],'color',color_list{3},'linewidth',2,'LineStyle','--');
load('F:\meg_data_copy\code_availability\data_availability\phase_coherence\Exp1\phy_coh_shuffle.mat');
phy_coh_shuffle = phy_sum;
p4 = plot(f,squeeze(mean(mean(phy_coh_shuffle(:,:,:),3),2)),'color','k','LineWidth',1.5,'LineStyle',':');

for i = 1:21
    a = mean(mean(phy_coh_clean(i,tidu,con,sub),2),4);
    b = squeeze(mean(phy_coh_shuffle(i,:,:),3));
    M_value(i) = (length(find(a<b))+1)/(size(phy_coh_shuffle,2)+1);
%     p_value = length(find(a<b))/size(phy_coh_shuffle,2);
end
% 
[h,q,p_new] = fdr_bh(M_value);
for i = 1:21
    p_value = p_new(i);
    clear a b
    if p_value<0.05
        f4 = fill([[f(i)-0.25 f(i)+0.25],[f(i)+0.25 f(i)-0.25]],...
            [0.0102+0.0005*con+0.00015,0.0102+0.0005*con+0.00015,0.0102+0.0005*con-0.00015,0.0102+0.0005*con-0.00015],'k');
        f4.EdgeColor = 'none';
        f4.FaceColor = color_list{con};
        f4.FaceAlpha = 0.3; 
    end
end
end
%% topography at echo-sensitive frequencies
clear
clc
close all
load('.\data_availability\phase_coherence\Exp1\phy_coh_clean.mat');

Color_temp = [[46 75 160]./255;[250 164 25]./255];
addpath '.\tools'
cili = 1:3:306;
fs = 100;
f = 1:size(phy_coh_clean,1);f=f-1;f=f/size(phy_coh_clean,1)*fs;
phy_coh_tidu = (phy_coh_clean(:,2:3:306,:,:)+phy_coh_clean(:,3:3:306,:,:))/2;
condition = {'Anechoic','0.125s echoic','0.25s echoic'};
for fre = [2 4 6]
for con = 1:3
    figure
    
    PKUtopoGradOutline(mean(phy_coh_tidu(find(f == fre),:,con,:),4)./ ...
        max(mean(phy_coh_tidu(find(f == fre),:,con,:),4)),3);
    colormap(flipud(othercolor('RdBu4')))
    set(gcf,'unit','centimeters','position',[1.5 1.5 6 4.8])
    colorbar
    title([condition{con} '-' num2str(fre) 'Hz'])
    set(gca,'fontname','arial')
    set(gca,'fontsize',13);

    clim([0 1])
end
end

%% MEG bar plot
clear
clc

% load('..\pilot_ProData\phy_coh\phy_coh_clean_0.8-10_1.875.mat');
load('.\data_availability\phase_coherence\Exp1\phy_coh_clean.mat');
phy_coh1 = phy_coh_clean;
load('.\data_availability\phase_coherence\Exp2\phy_coh_clean.mat');
phy_coh2 = phy_coh_clean;
load('.\data_availability\phase_coherence\Exp1\phy_coh_shuffle.mat');
phy_coh_shuffle1 = mean(mean(phy_sum,3),2);
load('.\data_availability\phase_coherence\Exp2\phy_coh_shuffle.mat');
phy_coh_shuffle2 = mean(mean(phy_sum,3),2);


fs = 100;
f = 1:size(phy_coh1,1);f=f-1;f=f/size(phy_coh1,1)*fs;

i = 1;
s_sum1 = [];
err1_sum = [];
s_sum2 = [];
err2_sum = [];
fre_range = [2 4 6];


for frebin = fre_range
fp = find(f == frebin);
tidu = setdiff(1:306,1:3:306);

pc_sum1 = squeeze(mean(phy_coh1(fp,tidu,[2 1 3],:),2));
pc_sum2 = squeeze(mean(phy_coh2(fp,tidu,[2 1 3],:),2));

err1 = std(pc_sum1,0,2)./sqrt(15);
err2 = std(pc_sum2,0,2)./sqrt(12);
pc_sum1_sum(:,:,i) = pc_sum1;
pc_sum2_sum(:,:,i) = pc_sum2;
s_sum1 = [s_sum1 squeeze(mean(pc_sum1,2))];
err1_sum = [err1_sum err1];
s_sum2 = [s_sum2 squeeze(mean(pc_sum2,2))];
err2_sum = [err2_sum err2];
i = i+1;
end

s_sum1 = s_sum1';
s_sum2 = s_sum2';

% experiment 1
figure;
b = bar(s_sum1,'FaceColor','flat','LineStyle','-','LineWidth',1,'BarWidth',0.8);

for bbb = 1:size(b,2)
    x_pos(bbb,:) = b(:,bbb).XEndPoints;
end
hold on

% individual dots lines
for i = 1:length(fre_range)
plot(x_pos(:,i),pc_sum1_sum(:,:,i),'color',[195 195 195]./255,'markersize',10)
end

% chance level
i = 1;
for frebin = fre_range
fp = find(f == frebin);
plot(x_pos(:,i),ones(3,1)*phy_coh_shuffle1(fp),'color','k','LineStyle',':','LineWidth',1)
i = i+1;
end

% individual dots
for i = 1:length(fre_range)
scatter(x_pos(:,i),pc_sum1_sum(:,:,i),20,'filled','MarkerFaceColor',[130 130 130]./255,'MarkerEdgeColor',[130 130 130]./255,'LineWidth',0.5);
end

% errorbar
e = errorbar(x_pos,s_sum1',err1_sum,'linestyle','none','LineWidth',1,'CapSize',10,'Color','k');
%
ylabel('coherence')
title(['Experiment 1'])
set(gca,'xticklabel',{'2 Hz','4 Hz','6 Hz'})
set(gca,'fontsize',13)

box off
ylim([0 0.015])
set(gcf,'unit','centimeters','position',[3 3 18 12])


%% exp 1 (significant test)
clear
clc
fs = 100;

load('.\data_availability\phase_coherence\Exp1\phy_coh_shuffle.mat');
phy_coh_shuffle1 = phy_sum;
load('.\data_availability\phase_coherence\Exp1\phy_coh_clean.mat');
phy_coh1 = phy_coh_clean;

tidu = setdiff(1:306,1:3:306);
f = 1:size(phy_coh1,1);f=f-1;f=f/size(phy_coh1,1)*fs;

% permutation
randddd = dec2bin(0:(2^15-1), 15);
rd = double(randddd)-48;
for i = 1:size(rd,1)
%     tic
    label = rd(i,:);
    pos = find(label == 1);
    phy_coh1_chance = squeeze(mean(phy_coh1(:,tidu,[1 2],:),2));
    phy_coh1_chance(:,:,pos) = squeeze(mean(phy_coh1(:,tidu,[2 1],pos),2));
    pc1_chance_125_sum(:,:,:,i) = phy_coh1_chance;
    clear phy_coh1_chance
    
%     toc
end
pc1_chance_125_sum = single(pc1_chance_125_sum);
% save('pilot_ProData\phy_coh_chance\pc1_chance_125.mat','pc1_chance_125_sum');


for i = 1:size(rd,1)
%     tic
    label = rd(i,:);
    pos = find(label == 1);
    phy_coh1_chance = squeeze(mean(phy_coh1(:,tidu,[1 3],:),2));
    phy_coh1_chance(:,:,pos) = squeeze(mean(phy_coh1(:,tidu,[3 1],pos),2));
    pc1_chance_250_sum(:,:,:,i) = phy_coh1_chance;
    clear phy_coh1_chance
%     toc
end
pc1_chance_250_sum = single(pc1_chance_250_sum);
% save('pilot_ProData\phy_coh_chance\pc1_chance_250.mat','pc1_chance_250_sum');

for i = 1:size(rd,1)
%     tic
    label = rd(i,:);
    pos = find(label == 1);
    phy_coh1_chance = squeeze(mean(phy_coh1(:,tidu,[2 3],:),2));
    phy_coh1_chance(:,:,pos) = squeeze(mean(phy_coh1(:,tidu,[3 2],pos),2));
    pc1_chance_125_250_sum(:,:,:,i) = phy_coh1_chance;
    clear phy_coh1_chance
    
%     toc
end
pc1_chance_125_250_sum = single(pc1_chance_125_250_sum);
% save('pilot_ProData\phy_coh_chance\pc1_chance_125_250.mat','pc1_chance_125_250_sum');


fff = 1;
for f_hz = [2 4 6]
for con = 1:3
    a = mean(mean(phy_coh1(find(f == f_hz),tidu,con,:),2),4);
    b = squeeze(mean(phy_coh_shuffle1(find(f == f_hz),:,con),3));
    M1_value(con) = (length(find(a<b))+1)/(size(phy_coh_shuffle1,2)+1);
end

% 125 condition
cc1 = squeeze(pc1_chance_125_sum(find(f == f_hz),1,:,:));
cc2 = squeeze(pc1_chance_125_sum(find(f == f_hz),2,:,:));
cc_con1(:,1) = mean(cc2-cc1,1);
real_con1(1) = mean(mean(phy_coh1(find(f == f_hz),tidu,2,:),2)-mean(phy_coh1(find(f == f_hz),tidu,1,:),2),4);
p1(fff,1) = (min(length(find(cc_con1(:,1)>real_con1(1))),length(find(cc_con1(:,1)<=real_con1(1))))+1)./(size(cc_con1,1)+1);

% 250 condition
cc1 = squeeze(pc1_chance_250_sum(find(f == f_hz),1,:,:));
cc2 = squeeze(pc1_chance_250_sum(find(f == f_hz),2,:,:));
cc_con1(:,2) = mean((cc2-cc1),1);
real_con1(2) = mean(mean(phy_coh1(find(f == f_hz),tidu,3,:),2)-mean(phy_coh1(find(f == f_hz),tidu,1,:),2),4);
p1(fff,2) = (min(length(find(cc_con1(:,2)>real_con1(2))),length(find(cc_con1(:,2)<=real_con1(2))))+1)./(size(cc_con1,1)+1);

% 125 & 250 condition
cc1 = squeeze(pc1_chance_125_250_sum(find(f == f_hz),1,:,:));
cc2 = squeeze(pc1_chance_125_250_sum(find(f == f_hz),2,:,:));
cc_con1(:,3) = mean((cc2-cc1),1);
real_con1(3) = mean(mean(phy_coh1(find(f == f_hz),tidu,3,:),2)-mean(phy_coh1(find(f == f_hz),tidu,2,:),2),4);
p1(fff,3) = (min(length(find(cc_con1(:,3)>real_con1(3))),length(find(cc_con1(:,3)<=real_con1(3))))+1)./(size(cc_con1,1)+1);
p1(fff,4:6) = M1_value;
fff = fff+1;
end
[~,~,p1_new_new] = fdr_bh(p1);