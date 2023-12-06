%% Experiment 1
% -----------------------------------------------------------
% nodepression
clear
clc
% close all
fs = 100;
tik = 0.01;
d = fs*1;
lag = 'all';
shift = 0.2;
i = 1;
tidu = setdiff(1:306,1:3:306);

load(['.\data_availability\TRF\Exp1\idealized_model\wid_trf.mat']);
pp_sum0 = squeeze(pp_new2_2_sum);
h_sum0 = h_new2_2_sum;
load(['.\data_availability\TRF\Exp1\mixture_model\wid_trf.mat']);
pp_sum1 = squeeze(pp_new2_3_sum);
h_sum1 = h_new2_3_sum;
load(['.\data_availability\TRF\Exp1\streaming_model\wid_trf.mat']);
pp_sum2 = squeeze(pp_new2_4_sum);
h_sum2 = h_new2_4_sum;
pp_sum = squeeze(cat(3,squeeze(mean(pp_sum1(tidu,:),1)),squeeze(mean(pp_sum2(tidu,:),1)),squeeze(mean(pp_sum0(tidu,:),1))));

pp_nodep = pp_sum;

% depression

load(['.\data_availability\TRF\Ada+TRF\Exp1\ada_idealized_model\wid_trf.mat']);
pp_sum0 = squeeze(pp_new2_2_sum);
load(['.\data_availability\TRF\Ada+TRF\Exp1\ada_mixture_model\wid_trf.mat']);
pp_sum1 = squeeze(pp_new2_3_sum);
load(['.\data_availability\TRF\Ada+TRF\Exp1\ada_streaming_model\wid_trf.mat']);
pp_sum2 = squeeze(pp_new2_4_sum);
pp_sum = squeeze(cat(3,squeeze(mean(pp_sum1(tidu,:),1)),squeeze(mean(pp_sum2(tidu,:),1)),squeeze(mean(pp_sum0(tidu,:),1))));
pp_dep = pp_sum;


figure
tidu = setdiff(1:306,1:3:306);
pp_sum = pp_dep-pp_nodep;
err = squeeze(std(pp_sum,0,1)/sqrt(size(pp_sum,1)));
e.CapSize = 10;
b = bar(squeeze(mean(pp_sum,1)));
hold on
for i = 1:size(err,1)
    xx(:,i) = b(1,i).XEndPoints';
end
b(1).FaceColor = 'none';
b(1).LineWidth = 2;
set(gca,'fontname','arial')
set(gca,'fontsize',15);
ylabel('Correlation')
s1 = scatter(ones(15,1)+0.1*randn(15,1),pp_sum(:,1),'MarkerFaceColor',[130 130 130]./255,'MarkerEdgeColor','none','SizeData',40)
s2 = scatter(ones(15,1)*2+0.1*randn(15,1),pp_sum(:,2),'MarkerFaceColor',[130 130 130]./255,'MarkerEdgeColor','none','SizeData',40)
s3 = scatter(ones(15,1)*3+0.1*randn(15,1),pp_sum(:,3),'MarkerFaceColor',[130 130 130]./255,'MarkerEdgeColor','none','SizeData',40)
errorbar(xx,squeeze(mean(pp_sum,1)),err,'linestyle','none','linewidth',2,'color','k','CapSize',10)
box off
set(gcf,'Units','centimeters','Position',[1.5 1.5 12 10]);
set(gca,'xticklabel',{'Mixture','Stream','Ideal'})
ylim([-1e-3 6.5e-3])
title("Predictive power difference(exp1)")

%% Experiment 2
% -------------------------------------------------------

% nodep
clear
clc
% close all
fs = 100;
d = fs*1;
lag = 'all';
shift = 0.2;


load(['.\data_availability\TRF\Exp2\idealized_model\wid_trf.mat']);
pp_sum0 = pp_new2_2_sum(:,1:12);
h_sum0 = h_new2_2_sum(:,:,:,1:12);
load(['.\data_availability\TRF\Exp2\mixture_model\wid_trf.mat']);
pp_sum1 = pp_new2_3_sum(:,1:12);
h_sum1 = h_new2_3_sum(:,:,:,1:12);
load(['.\data_availability\TRF\Exp2\streaming_model\wid_trf.mat']);
pp_sum2 = pp_new2_4_sum(:,1:12);
h_sum2 = h_new2_4_sum(:,:,:,1:12);

tidu = setdiff(1:306,1:3:306);
pp_sum = squeeze(cat(3,squeeze(mean(pp_sum1(tidu,:),1)),squeeze(mean(pp_sum2(tidu,:),1)),squeeze(mean(pp_sum0(tidu,:),1))));
pp_nodep = pp_sum;


% dep

load(['.\data_availability\TRF\Ada+TRF\Exp2\ada_idealized_model\wid_trf.mat']);
pp_sum0 = pp_new2_2_sum(:,1:12);
load(['.\data_availability\TRF\Ada+TRF\Exp2\ada_mixture_model\wid_trf.mat']);
pp_sum1 = pp_new2_3_sum(:,1:12);
load(['.\data_availability\TRF\Ada+TRF\Exp2\ada_streaming_model\wid_trf.mat']);
pp_sum2 = pp_new2_4_sum(:,1:12);

tidu = setdiff(1:306,1:3:306);
pp_sum = squeeze(cat(3,squeeze(mean(pp_sum1(tidu,:),1)),squeeze(mean(pp_sum2(tidu,:),1)),squeeze(mean(pp_sum0(tidu,:),1))));
pp_dep = pp_sum;
figure
tidu = setdiff(1:306,1:3:306);
pp_sum = pp_dep-pp_nodep;
err = squeeze(std(pp_sum,0,1)/sqrt(size(pp_sum,1)));
e.CapSize = 10;
b = bar(squeeze(mean(pp_sum,1)));
hold on
for i = 1:size(err,1)
    xx(:,i) = b(1,i).XEndPoints';
end
b(1).FaceColor = 'none';
b(1).LineWidth = 2;
set(gca,'fontname','arial')
set(gca,'fontsize',15);
ylabel('Correlation')
s1 = scatter(ones(12,1)+0.1*randn(12,1),pp_sum(:,1),'MarkerFaceColor',[130 130 130]./255,'MarkerEdgeColor','none','SizeData',40)
s2 = scatter(ones(12,1)*2+0.1*randn(12,1),pp_sum(:,2),'MarkerFaceColor',[130 130 130]./255,'MarkerEdgeColor','none','SizeData',40)
s3 = scatter(ones(12,1)*3+0.1*randn(12,1),pp_sum(:,3),'MarkerFaceColor',[130 130 130]./255,'MarkerEdgeColor','none','SizeData',40)
errorbar(xx,squeeze(mean(pp_sum,1)),err,'linestyle','none','linewidth',2,'color','k','CapSize',10)
box off
set(gcf,'Units','centimeters','Position',[1.5 1.5 12 10]);
set(gca,'xticklabel',{'Mixture','Stream','Ideal'})
ylim([-1e-3 6e-3])
title("Predictive power difference(exp2)")

%% Experiment 3
% -------------------------------------------------------

% nodepression
clear
clc
fs = 100;
tik = 0.1;
d = fs*1;
lag = 'all';
i = 1;

load(['.\data_availability\TRF\Exp3\idealized_model\wid_trf.mat']);
pp_sum0 = pp_new2_2_sum;
hh_sum0 = h_new2_2_sum;
load(['.\data_availability\TRF\Exp3\mixture_model\wid_trf.mat']);
pp_sum1 = pp_new2_3_sum;
hh_sum1 = h_new2_3_sum;

load(['.\data_availability\TRF\Exp3\streaming_model\wid_trf.mat']);
pp_sum2 = pp_new2_4_sum;
hh_sum2 = h_new2_4_sum;


tidu = setdiff(1:306,1:3:306);
% tidu = 2:3:306;
pp_sum = cat(3,squeeze(mean(pp_sum1(tidu,:,:),1)),squeeze(mean(pp_sum2(tidu,:,:),1)),squeeze(mean(pp_sum0(tidu,:,:),1)));
pp_sum_nodep = pp_sum;


% depression

load(['.\data_availability\TRF\Ada+TRF\Exp3\ada_idealized_model\att\wid_trf.mat']);
pp_sum0(:,1,:) = pp_new2_2_sum(:,1,:);
load(['.\data_availability\TRF\Ada+TRF\Exp3\ada_idealized_model\unatt\wid_trf.mat']);
pp_sum0(:,2,:) = pp_new2_2_sum(:,2,:);
load(['.\data_availability\TRF\Ada+TRF\Exp3\ada_idealized_model\vo\wid_trf.mat']);
pp_sum0(:,3,:) = pp_new2_2_sum(:,3,:);



load(['.\data_availability\TRF\Ada+TRF\Exp3\ada_mixture_model\att\wid_trf.mat']);
pp_sum1(:,1,:) = pp_new2_3_sum(:,1,:);
load(['.\data_availability\TRF\Ada+TRF\Exp3\ada_mixture_model\unatt\wid_trf.mat']);
pp_sum1(:,2,:) = pp_new2_3_sum(:,2,:);
load(['.\data_availability\TRF\Ada+TRF\Exp3\ada_mixture_model\vo\wid_trf.mat']);
pp_sum1(:,3,:) = pp_new2_3_sum(:,3,:);


load(['.\data_availability\TRF\Ada+TRF\Exp3\ada_streaming_model\att\wid_trf.mat']);
pp_sum2(:,1,:) = pp_new2_4_sum(:,1,:);
load(['.\data_availability\TRF\Ada+TRF\Exp3\ada_streaming_model\unatt\wid_trf.mat']);
pp_sum2(:,2,:) = pp_new2_4_sum(:,2,:);
load(['.\data_availability\TRF\Ada+TRF\Exp3\ada_streaming_model\vo\wid_trf.mat']);
pp_sum2(:,3,:) = pp_new2_4_sum(:,3,:);

tidu = setdiff(1:306,1:3:306);
% tidu = 2:3:306;
pp_sum = cat(3,squeeze(mean(pp_sum1(tidu,:,:),1)),squeeze(mean(pp_sum2(tidu,:,:),1)),squeeze(mean(pp_sum0(tidu,:,:),1)));
pp_sum_dep = pp_sum;

pp_sum = pp_sum_dep-pp_sum_nodep;

err = squeeze(std(pp_sum,0,2)/sqrt(size(pp_sum,2)));

for con = 1:3
figure
b = bar(squeeze(mean(pp_sum(con,:,:),2)));
hold on


b(1).FaceColor = 'none';
b(1).LineWidth = 2;
b(1).BarWidth = 0.8;
% plot([1 2 3],squeeze(pp_sum(con,:,:))','color',[195 195 195]./255)

for i = 1:21
s1 = scatter(1+0.1*randn(1),pp_sum(con,i,1),'MarkerFaceColor',[130 130 130]./255,'MarkerEdgeColor','none','SizeData',40);
hold on
s2 = scatter(2+0.1*randn(1),pp_sum(con,i,2),'MarkerFaceColor',[130 130 130]./255,'MarkerEdgeColor','none','SizeData',40);
s3 = scatter(3+0.1*randn(1),pp_sum(con,i,3),'MarkerFaceColor',[130 130 130]./255,'MarkerEdgeColor','none','SizeData',40);
end
e = errorbar([1 2 3],squeeze(mean(pp_sum(con,:,:),2)),err(con,:),'linestyle','none','linewidth',2,'color','k')
e.CapSize = 10;

set(gca,'xticklabel',{'Mixture','Stream','Ideal'})
set(gca,'fontname','arial')
set(gca,'fontsize',15);
set(gcf,'unit','centimeters','position',[1.50 1.5 12 10])
ylabel('Correlation')
title('Predictive power difference')
% set(gca,'ytick',[-0.02:0.02:0.1])
% ylim([0 0.1])
box off
ylim([-1.8e-3 6e-3])
end