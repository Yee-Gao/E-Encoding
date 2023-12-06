%% MEG bar plot (1 Hz)  
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
f = 1:size(phy_coh2,1);f=f-1;f=f/size(phy_coh1,1)*fs;

i = 1;
for frebin = [1]
fp = find(f == frebin);
tidu = setdiff(1:306,1:3:306);

pc_sum1 = squeeze(mean(phy_coh1(fp,tidu,[2 1 3],:),2));
pc_sum2 = squeeze(mean(phy_coh2(fp,tidu,[2 1 3],:),2));

err1 = std(pc_sum1,0,2)./sqrt(15);
err2 = std(pc_sum2,0,2)./sqrt(12);

s_sum =cat(2,squeeze(mean(pc_sum1,2)),squeeze(mean(pc_sum2,2)));
s_sum = s_sum';

err = [err1,err2];

figure
b = bar(s_sum,'FaceColor','flat','LineStyle','-','LineWidth',1,'BarWidth',0.8);
% b = bar(s_sum','FaceColor','flat','LineStyle','-','LineWidth',1,'BarWidth',0.8);

for bbb = 1:size(b,2)
    x_pos(bbb,:) = b(:,bbb).XEndPoints;
end
hold on

plot(x_pos(:,1),pc_sum1,'color',[195 195 195]./255,'markersize',10)
plot(x_pos(:,2),pc_sum2,'color',[195 195 195]./255,'markersize',10)
plot(x_pos(:,1),ones(3,1)*phy_coh_shuffle1(fp),'color','k','LineStyle','-','LineWidth',1)
plot(x_pos(:,2),ones(3,1)*phy_coh_shuffle2(fp),'color','k','LineStyle','-','LineWidth',1)

scatter(x_pos(:,1),pc_sum1,20,'filled','MarkerFaceColor',[130 130 130]./255,'MarkerEdgeColor',[130 130 130]./255,'LineWidth',0.5);
scatter(x_pos(:,2),pc_sum2,20,'filled','MarkerFaceColor',[130 130 130]./255,'MarkerEdgeColor',[130 130 130]./255,'LineWidth',0.5);

e = errorbar(x_pos,s_sum',err,'linestyle','none','LineWidth',1,'CapSize',10,'Color','k');
ylabel('coherence')
title([num2str(frebin) 'Hz'])
set(gca,'xticklabel',{'Exp 1','Exp 2'})
set(gca,'fontsize',13)
if i == 4
legend('0.125-s echoic speech','anechoic speech','0.25-s echoic speech')
end
box off
ylim([0 0.023])
% grid on
i = i+1;

end
set(gcf,'unit','centimeters','position',[3 3 15 10])
