%% Experiment 1
% 0.8-5
clear
clc
% close all
fs = 100;

mm = 1;

i = 1;
for sub = [1:15]
    load(['.\data_availability\TRF\Narrow_TRF\Exp1\0.5-5\idealized_model\wid_trf_' num2str(sub) '.mat']);
    pp_sum0(:,i) = pp_new2_2;
    h_sum0(:,:,:,i) = h_new2_2;
    load(['.\data_availability\TRF\Narrow_TRF\Exp1\0.5-5\mixture_model\wid_trf_' num2str(sub) '.mat']);
    pp_sum1(:,i) = pp_new2_3;
    h_sum1(:,:,:,i) = h_new2_3;
    load(['.\data_availability\TRF\Narrow_TRF\Exp1\0.5-5\streaming_model\wid_trf_' num2str(sub) '.mat']);
    pp_sum2(:,i) = pp_new2_4;
    h_sum2(:,:,:,i) = h_new2_4;
    i = i+1;
end
figure
tidu = setdiff(1:306,1:3:306);
pp_sum = squeeze(cat(3,squeeze(mean(pp_sum1(tidu,:),1)),squeeze(mean(pp_sum2(tidu,:),1)),squeeze(mean(pp_sum0(tidu,:),1))));
pp_sum_sum(:,:,mm) = pp_sum;
mm = mm+1;


err = squeeze(std(pp_sum,0,1)/sqrt(size(pp_sum,1)));
b = bar(squeeze(mean(pp_sum,1)));
hold on
for i = 1:size(err,1)
    xx(:,i) = b(1,i).XEndPoints';
end
b(1).FaceColor = 'none';
b(1).LineWidth = 2;
set(gca,'fontname','arial')
set(gca,'fontsize',15);
set(gcf,'unit','centimeters','position',[1.50 1.5 15 10])
ylabel('predictive power')
title(['Exp1:[0.8-5 Hz]'])
ylim([0 0.1])
set(gca,'ytick',0:0.02:0.4)
set(gca,'xticklabel',{'base','stream','idealized'})
plot(pp_sum','color',[195 195 195]./255,'linewidth',1)
s1 = scatter(ones(15,1),pp_sum(:,1),'MarkerFaceColor',[250 164 25]./255,'MarkerEdgeColor','none','SizeData',40)
s2 = scatter(ones(15,1)*2,pp_sum(:,2),'MarkerFaceColor','#228B22','MarkerEdgeColor','none','SizeData',40)
s3 = scatter(ones(15,1)*3,pp_sum(:,3),'MarkerFaceColor',[46 75 160]./255,'MarkerEdgeColor','none','SizeData',40)
errorbar(xx,squeeze(mean(pp_sum,1)),err,'linestyle','none','linewidth',2,'color','k','CapSize',10)
box off

% 5-10
clear
clc
% close all
fs = 100;

mm = 1;

i = 1;
for sub = [1:15]
    load(['.\data_availability\TRF\Narrow_TRF\Exp1\5-10\idealized_model\wid_trf_' num2str(sub) '.mat']);
    pp_sum0(:,i) = pp_new2_2;
    h_sum0(:,:,:,i) = h_new2_2;
    load(['.\data_availability\TRF\Narrow_TRF\Exp1\5-10\mixture_model\wid_trf_' num2str(sub) '.mat']);
    pp_sum1(:,i) = pp_new2_3;
    h_sum1(:,:,:,i) = h_new2_3;
    load(['.\data_availability\TRF\Narrow_TRF\Exp1\5-10\streaming_model\wid_trf_' num2str(sub) '.mat']);
    pp_sum2(:,i) = pp_new2_4;
    h_sum2(:,:,:,i) = h_new2_4;
    i = i+1;
end
figure
tidu = setdiff(1:306,1:3:306);
pp_sum = squeeze(cat(3,squeeze(mean(pp_sum1(tidu,:),1)),squeeze(mean(pp_sum2(tidu,:),1)),squeeze(mean(pp_sum0(tidu,:),1))));
pp_sum_sum(:,:,mm) = pp_sum;
mm = mm+1;


err = squeeze(std(pp_sum,0,1)/sqrt(size(pp_sum,1)));
b = bar(squeeze(mean(pp_sum,1)));
hold on
for i = 1:size(err,1)
    xx(:,i) = b(1,i).XEndPoints';
end
b(1).FaceColor = 'none';
b(1).LineWidth = 2;
set(gca,'fontname','arial')
set(gca,'fontsize',15);
set(gcf,'unit','centimeters','position',[1.50 1.5 15 10])
ylabel('predictive power')
title(['Exp1:[5-10 Hz]'])
ylim([0 0.045])
set(gca,'ytick',0:0.02:0.4)
set(gca,'xticklabel',{'base','stream','idealized'})
plot(pp_sum','color',[195 195 195]./255,'linewidth',1)
s1 = scatter(ones(15,1),pp_sum(:,1),'MarkerFaceColor',[250 164 25]./255,'MarkerEdgeColor','none','SizeData',40)
s2 = scatter(ones(15,1)*2,pp_sum(:,2),'MarkerFaceColor','#228B22','MarkerEdgeColor','none','SizeData',40)
s3 = scatter(ones(15,1)*3,pp_sum(:,3),'MarkerFaceColor',[46 75 160]./255,'MarkerEdgeColor','none','SizeData',40)
errorbar(xx,squeeze(mean(pp_sum,1)),err,'linestyle','none','linewidth',2,'color','k','CapSize',10)
box off
%% Experiment 2

% 0.8-5
clear
clc
% close all
fs = 100;
mm = 1;

i = 1;
for sub = [1:12]
    load(['.\data_availability\TRF\Narrow_TRF\Exp2\0.5-5\idealized_model\wid_trf_' num2str(sub) '.mat']);
    pp_sum0(:,i) = pp_new2_2;
    h_sum0(:,:,:,i) = h_new2_2;
    load(['.\data_availability\TRF\Narrow_TRF\Exp2\0.5-5\mixture_model\wid_trf_' num2str(sub) '.mat']);
    pp_sum1(:,i) = pp_new2_3;
    h_sum1(:,:,:,i) = h_new2_3;
    load(['.\data_availability\TRF\Narrow_TRF\Exp2\0.5-5\streaming_model\wid_trf_' num2str(sub) '.mat']);
    pp_sum2(:,i) = pp_new2_4;
    h_sum2(:,:,:,i) = h_new2_4;
    i = i+1;
end
figure
tidu = setdiff(1:306,1:3:306);
pp_sum = squeeze(cat(3,squeeze(mean(pp_sum1(tidu,:),1)),squeeze(mean(pp_sum2(tidu,:),1)),squeeze(mean(pp_sum0(tidu,:),1))));
pp_sum_sum(:,:,mm) = pp_sum;

err = squeeze(std(pp_sum,0,1)/sqrt(size(pp_sum,1)));
b = bar(squeeze(mean(pp_sum,1)));
hold on
for i = 1:size(err,1)
    xx(:,i) = b(1,i).XEndPoints';
end
b(1).FaceColor = 'none';
b(1).LineWidth = 2;
set(gca,'fontname','arial')
set(gca,'fontsize',15);
set(gcf,'unit','centimeters','position',[1.50 1.5 15 10])
ylabel('predictive power')
ylim([0 0.12])
title(['Exp2:[0.8-5 Hz]'])
set(gca,'ytick',[-0.02:0.02:0.12])
set(gca,'xticklabel',{'base','stream','idealized'})
plot(pp_sum','color',[195 195 195]./255,'linewidth',1)
s1 = scatter(ones(12,1),pp_sum(:,1),'MarkerFaceColor',[250 164 25]./255,'MarkerEdgeColor','none','SizeData',40)
s2 = scatter(ones(12,1)*2,pp_sum(:,2),'MarkerFaceColor','#228B22','MarkerEdgeColor','none','SizeData',40)
s3 = scatter(ones(12,1)*3,pp_sum(:,3),'MarkerFaceColor',[46 75 160]./255,'MarkerEdgeColor','none','SizeData',40)
errorbar(xx,squeeze(mean(pp_sum,1)),err,'linestyle','none','linewidth',2,'color','k','CapSize',10)
box off


% 5-10
clear
clc
fs = 100;
mm = 1;

i = 1;
for sub = [1:12]
    load(['.\data_availability\TRF\Narrow_TRF\Exp2\5-10\idealized_model\wid_trf_' num2str(sub) '.mat']);
    pp_sum0(:,i) = pp_new2_2;
    h_sum0(:,:,:,i) = h_new2_2;
    load(['.\data_availability\TRF\Narrow_TRF\Exp2\5-10\mixture_model\wid_trf_' num2str(sub) '.mat']);
    pp_sum1(:,i) = pp_new2_3;
    h_sum1(:,:,:,i) = h_new2_3;
    load(['.\data_availability\TRF\Narrow_TRF\Exp2\5-10\streaming_model\wid_trf_' num2str(sub) '.mat']);
    pp_sum2(:,i) = pp_new2_4;
    h_sum2(:,:,:,i) = h_new2_4;
    i = i+1;
end
figure
tidu = setdiff(1:306,1:3:306);
pp_sum = squeeze(cat(3,squeeze(mean(pp_sum1(tidu,:),1)),squeeze(mean(pp_sum2(tidu,:),1)),squeeze(mean(pp_sum0(tidu,:),1))));
pp_sum_sum(:,:,mm) = pp_sum;

% pp_sum = squeeze(cat(3,squeeze(mean(pp_sum1(tidu,:),1)),squeeze(mean(pp_sum2(tidu,:),1))));
err = squeeze(std(pp_sum,0,1)/sqrt(size(pp_sum,1)));
b = bar(squeeze(mean(pp_sum,1)));
hold on
for i = 1:size(err,1)
    xx(:,i) = b(1,i).XEndPoints';
end
b(1).FaceColor = 'none';
b(1).LineWidth = 2;
set(gca,'fontname','arial')
set(gca,'fontsize',15);
set(gcf,'unit','centimeters','position',[1.50 1.5 15 10])
ylabel('predictive power')
ylim([0 0.04])
title(['Exp2:[5-10 Hz]'])
set(gca,'ytick',[-0.02:0.02:0.12])
set(gca,'xticklabel',{'base','stream','idealized'})
plot(pp_sum','color',[195 195 195]./255,'linewidth',1)
s1 = scatter(ones(12,1),pp_sum(:,1),'MarkerFaceColor',[250 164 25]./255,'MarkerEdgeColor','none','SizeData',40)
s2 = scatter(ones(12,1)*2,pp_sum(:,2),'MarkerFaceColor','#228B22','MarkerEdgeColor','none','SizeData',40)
s3 = scatter(ones(12,1)*3,pp_sum(:,3),'MarkerFaceColor',[46 75 160]./255,'MarkerEdgeColor','none','SizeData',40)
errorbar(xx,squeeze(mean(pp_sum,1)),err,'linestyle','none','linewidth',2,'color','k','CapSize',10)
box off

%% Experiment 3
% 0.8-5
clear
clc
exp_con = {'attend','unattend','vocode'};
fs = 100;

mm=1;

i = 1;
for sub = [1:9 12 14:24]

    load(['.\data_availability\TRF\Narrow_TRF\Exp3\0.5-5\idealized_model\wid_trf_' num2str(sub) '.mat']);
    pp_sum0(:,:,i) = pp_new2_2;
    hh_sum0(:,:,:,:,i) = h_new2_2;
    load(['.\data_availability\TRF\Narrow_TRF\Exp3\0.5-5\mixture_model\wid_trf_' num2str(sub) '.mat']);
    pp_sum1(:,:,i) = pp_new2_3;
    hh_sum1(:,:,:,:,i) = h_new2_3;
    load(['.\data_availability\TRF\Narrow_TRF\Exp3\0.5-5\streaming_model\wid_trf_' num2str(sub) '.mat']);
    pp_sum2(:,:,i) = pp_new2_4;
    hh_sum2(:,:,:,:,i) = h_new2_4;
    i = i+1;
end

tidu = setdiff(1:306,1:3:306);

pp_sum = cat(3,squeeze(mean(pp_sum1(tidu,:,:),1)),squeeze(mean(pp_sum2(tidu,:,:),1)),squeeze(mean(pp_sum0(tidu,:,:),1)));
pp_sum_sum(:,:,:,mm) = pp_sum;
mm = mm+1;

err = squeeze(std(pp_sum,0,2)/sqrt(size(pp_sum,2)));

for con = 1:3
figure
b = bar(squeeze(mean(pp_sum(con,:,:),2)));
hold on


b(1).FaceColor = 'none';
b(1).LineWidth = 2;
b(1).BarWidth = 0.8;
plot([1 2 3],squeeze(pp_sum(con,:,:))','color',[195 195 195]./255)

for i = 1:21
s1 = scatter(1,pp_sum(con,i,1),'MarkerFaceColor',[250 164 25]./255,'MarkerEdgeColor','none','SizeData',40);
hold on
s2 = scatter(2,pp_sum(con,i,2),'MarkerFaceColor','#228B22','MarkerEdgeColor','none','SizeData',40);
s3 = scatter(3,pp_sum(con,i,3),'MarkerFaceColor',[46 75 160]./255,'MarkerEdgeColor','none','SizeData',40);
end
e = errorbar([1 2 3],squeeze(mean(pp_sum(con,:,:),2)),err(con,:),'linestyle','none','linewidth',2,'color','k')
e.CapSize = 10;

set(gca,'xticklabel',{'baseline','streaming','idealized'})
set(gca,'fontname','arial')
set(gca,'fontsize',15);
set(gcf,'unit','centimeters','position',[1.50 1.5 15 10])
ylabel('predictive power')
set(gca,'ytick',[0:0.02:0.4])
ylim([0 0.1])
title([exp_con{con} ', Exp3:[0.8-5 Hz]'])

box off
% saveas(gcf,['.\figure\figure_supplement\pp_tidu_con_' num2str(con) '-' fre_band '.svg'])
end

% 5-10
clear
clc
exp_con = {'attend','unattend','vocode'};
fs = 100;
mm=1;

i = 1;
for sub = [1:9 12 14:24]

    load(['.\data_availability\TRF\Narrow_TRF\Exp3\5-10\idealized_model\wid_trf_' num2str(sub) '.mat']);
    pp_sum0(:,:,i) = pp_new2_2;
    hh_sum0(:,:,:,:,i) = h_new2_2;
    load(['.\data_availability\TRF\Narrow_TRF\Exp3\5-10\mixture_model\wid_trf_' num2str(sub) '.mat']);
    pp_sum1(:,:,i) = pp_new2_3;
    hh_sum1(:,:,:,:,i) = h_new2_3;
    load(['.\data_availability\TRF\Narrow_TRF\Exp3\5-10\streaming_model\wid_trf_' num2str(sub) '.mat']);
    pp_sum2(:,:,i) = pp_new2_4;
    hh_sum2(:,:,:,:,i) = h_new2_4;
    i = i+1;
end

tidu = setdiff(1:306,1:3:306);
pp_sum = cat(3,squeeze(mean(pp_sum1(tidu,:,:),1)),squeeze(mean(pp_sum2(tidu,:,:),1)),squeeze(mean(pp_sum0(tidu,:,:),1)));
pp_sum_sum(:,:,:,mm) = pp_sum;
mm = mm+1;

err = squeeze(std(pp_sum,0,2)/sqrt(size(pp_sum,2)));

for con = 1:3
figure
b = bar(squeeze(mean(pp_sum(con,:,:),2)));
hold on


b(1).FaceColor = 'none';
b(1).LineWidth = 2;
b(1).BarWidth = 0.8;
plot([1 2 3],squeeze(pp_sum(con,:,:))','color',[195 195 195]./255)

for i = 1:21
s1 = scatter(1,pp_sum(con,i,1),'MarkerFaceColor',[250 164 25]./255,'MarkerEdgeColor','none','SizeData',40);
hold on
s2 = scatter(2,pp_sum(con,i,2),'MarkerFaceColor','#228B22','MarkerEdgeColor','none','SizeData',40);
s3 = scatter(3,pp_sum(con,i,3),'MarkerFaceColor',[46 75 160]./255,'MarkerEdgeColor','none','SizeData',40);
end
e = errorbar([1 2 3],squeeze(mean(pp_sum(con,:,:),2)),err(con,:),'linestyle','none','linewidth',2,'color','k')
e.CapSize = 10;

set(gca,'xticklabel',{'baseline','streaming','idealized'})
set(gca,'fontname','arial')
set(gca,'fontsize',15);
set(gcf,'unit','centimeters','position',[1.50 1.5 15 10])
ylabel('predictive power')
set(gca,'ytick',[0:0.02:0.4])
ylim([0 0.04])
title([exp_con{con} ', Exp3:[5-10 Hz]'])

box off
% saveas(gcf,['.\figure\figure_supplement\pp_tidu_con_' num2str(con) '-' fre_band '.svg'])
end