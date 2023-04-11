%%  phase coherence
clear
clc
% load('..\pilot_ProData\phy_coh\phy_coh_clean_0.8-10_1.875.mat');
load('.\data_availability\phase_coherence\Exp1\phy_coh_clean.mat');

fs = 100;
cili = 1:3:306;
tidu = setdiff(1:306,cili);

f = 1:size(phy_coh_clean,1);f=f-1;f=f/size(phy_coh_clean,1)*fs;
sub = [1:15];
condition = {'Anechoic speech','0.125-s echoic speech','0.25-s echoic speech'};
for con = 1:3
    figure
    pc_2 = squeeze(mean(phy_coh_clean(:,tidu,1,sub),2));
    err_2 = std(pc_2,0,2)./sqrt(length(sub));
    f2 = fill([f,fliplr(f)],[mean(pc_2,2)+err_2;flipud(mean(pc_2,2)-err_2)],[195 195 195]./255);
    f2.EdgeColor = 'none';
    f2.FaceColor = 'k';
    f2.FaceAlpha = 0.3; 
    hold on
    p2 = plot(f,mean(mean(phy_coh_clean(:,tidu,1,sub),2),4),'linewidth',2,'color','k');

    if con ~= 1
        pc_1 = squeeze(mean(phy_coh_clean(:,tidu,con,sub),2));
        err_1 = std(pc_1,0,2)./sqrt(length(sub));
        f1 = fill([f,fliplr(f)],[mean(pc_1,2)+err_1;flipud(mean(pc_1,2)-err_1)],[195 195 195]./255);
        f1.EdgeColor = 'none';
        f1.FaceColor = [176 31 36]./255;
        f1.FaceAlpha = 0.3; 
        hold on
        p1 = plot(f,mean(mean(phy_coh_clean(:,tidu,con,sub),2),4),'linewidth',2,'color',[176 31 36]./255);
    end

    xlim([1 14])
    xlabel('Frequency(Hz)')
    ylabel('Phase coherence')
    set(gca,'fontname','arial')
    set(gca,'fontsize',13);
    set(gcf,'unit','centimeters','position',[1.50 1.5 10.5 8])
    title(condition{con})
    ylim([0.002 0.012])


if con == 2
    p3 = plot([4 4],[0 0.015],'color',[176 31 36]./255,'linestyle','--','linewidth',2);
elseif con == 3
    p3 = plot([2 2],[0 0.015],'color',[176 31 36]./255,'linestyle','--','linewidth',2);
    p3 = plot([6 6],[0 0.015],'color',[176 31 36]./255,'linestyle','--','linewidth',2);
end

load('.\data_availability\phase_coherence\Exp1\phy_coh_shuffle.mat');
phy_coh_shuffle = phy_sum;
p4 = plot(f,squeeze(mean(mean(phy_coh_shuffle(:,:,:),3),2)),'color','k','LineWidth',1.5,'LineStyle','--');

for i = 1:21
    a = mean(mean(phy_coh_clean(i,tidu,con,sub),2),4);
    b = squeeze(mean(phy_coh_shuffle(i,:,:),3));
    M_value(i) = (length(find(a<b))+1)/(size(phy_coh_shuffle,2)+1);

end

[h,q,p_new] = fdr_bh(M_value);
for i = 1:21
    p_value = p_new(i);
    clear a b
    if p_value<0.05
        f4 = fill([[f(i)-0.25 f(i)+0.25],[f(i)+0.25 f(i)-0.25]],...
            [0.0115+0.00015,0.0115+0.00015,0.0115-0.00015,0.0115-0.00015],'k');
        f4.EdgeColor = 'none';
        f4.FaceColor = [176 31 36]./255;
        f4.FaceAlpha = 0.3; 
    end
end
end


%% topograph at echo-sensitive frequencies
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
    set(gcf,'unit','centimeters','position',[1.5 1.5 10 8])
    colorbar
    title([condition{con} '-' num2str(fre) 'Hz'])
    set(gca,'fontname','arial')
    set(gca,'fontsize',13);

    clim([0 1])
end
end