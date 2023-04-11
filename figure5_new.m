%% Exp1(predictive power)
clear
clc
close all
fs = 100;
tik = 0.1;
d = fs*1;
lag = 'all';
addpath '.\tools'
i = 1;
shift = 0.2;
for sub = [1:15]
    load(['.\data_availability\TRF\Exp1\baseline_model\tik' ...
        num2str(tik) 'd' num2str(d) '\wid_trf_' num2str(sub) '.mat']);

    pp_sum1(:,i) = pp_new2_3;
    h_sum1(:,:,:,i) = h_new2_3;
    i = i+1;
end



tik = 0.09;
d = fs*1;
lag = 'all';
i = 1;
shift = 0.2;
for sub = [1:15]
    load(['.\data_availability\TRF\Exp1\streaming_model\tik' ...
        num2str(tik) 'd' num2str(d) '\wid_trf_' num2str(sub) '.mat']);
    pp_sum2(:,i) = pp_new2_4;
    h_sum2(:,:,:,i) = h_new2_4;
    i = i+1;
end



tik = 0.35;
d = fs*1;
lag = 'all';
i = 1;
shift = 0.2;
for sub = [1:15]
    load(['.\data_availability\TRF\Exp1\idealized_model\tik' ...
        num2str(tik) 'd' num2str(d) '\wid_trf_' num2str(sub) '.mat']);
    pp_sum0(:,i) = pp_new2_2;
    h_sum0(:,:,:,i) = h_new2_2;
    i = i+1;
end


tidu = setdiff(1:306,1:3:306);
pp_sum = squeeze(cat(3,squeeze(mean(pp_sum1(tidu,:),1)),squeeze(mean(pp_sum2(tidu,:),1)),squeeze(mean(pp_sum0(tidu,:),1))));
err = squeeze(std(pp_sum,0,1)/sqrt(size(pp_sum,1)));
b = bar(squeeze(mean(pp_sum,1)));
hold on
for i = 1:size(err,1)
    xx(:,i) = b(1,i).XEndPoints';
end
b(1).FaceColor = 'none';
b(1).LineWidth = 2;
b(1).BarWidth = 0.8;
set(gca,'fontname','arial')
set(gca,'fontsize',13);
set(gcf,'unit','centimeters','position',[1.50 1.5 8 10])
ylabel('predictive power')
ylim([0 0.1])
set(gca,'xticklabel',{'echo','first+second'})
plot(pp_sum','color',[195 195 195]./255,'linewidth',1)
s1 = scatter(ones(15,1),pp_sum(:,1),'MarkerFaceColor',[250 164 25]./255,'MarkerEdgeColor','none','SizeData',40)
s2 = scatter(ones(15,1)*2,pp_sum(:,2),'MarkerFaceColor',[46 75 160]./255,'MarkerEdgeColor','none','SizeData',40)
s3 = scatter(ones(15,1)*3,pp_sum(:,3),'MarkerFaceColor',[46 75 160]./255,'MarkerEdgeColor','none','SizeData',40)

e = errorbar(xx,squeeze(mean(pp_sum,1)),err,'linestyle','none','linewidth',2,'color','k')
e.CapSize = 10;

set(gca,'xticklabel',{'baseline','streaming','idealized'})
set(gca,'fontname','arial')
set(gca,'fontsize',15);
set(gcf,'unit','centimeters','position',[1.50 1.5 15 10])
ylabel('predictive power')
set(gca,'ytick',[-0.02:0.02:0.1])
box off

addpath 'E:\echo_gjx\others\code\tools'
Bv = 10000;
[O,A]=bootstrap_for_vector(Bv,0.05,@compare_mean,pp_sum(:,[1 2])); 
p1(1)=2*[min(sum(A<=0),sum(A>=0))+1]/[Bv+1];

[O,A]=bootstrap_for_vector(Bv,0.05,@compare_mean,pp_sum(:,[3 2])); 
p1(2)=2*[min(sum(A<=0),sum(A>=0))+1]/[Bv+1];
[h,q,p_new] = fdr_bh(p1)


%% TRF

fs = 100;
t = -20:1:79;t=t/fs;

figure
% subplot(1,2,2)
h1 = squeeze(rms(h_sum2(:,1,tidu,:),3))-mean(squeeze(rms(h_sum2(1:21,1,tidu,:),3)));
err_1 = std(h1,0,2)./sqrt(size(h1,2));
h2 = squeeze(rms(h_sum2(:,2,tidu,:),3))-mean(squeeze(rms(h_sum2(1:21,2,tidu,:),3)));
err_2 = std(h1,0,2)./sqrt(size(h2,2));
f1 = fill([t fliplr(t)],[rms(h1,2)+err_1;flipud(rms(h1,2)-err_1)],[195 195 195]./255);
hold on
f1.EdgeColor = 'none';
f1.FaceColor = [250 164 25]./255;
f1.FaceAlpha = 0.3;
f2 = fill([t fliplr(t)],[rms(h2,2)+err_2;flipud(rms(h2,2)-err_2)],[195 195 195]./255);
hold on
f2.EdgeColor = 'none';
f2.FaceColor = [46 75 160]./255;
f2.FaceAlpha = 0.3;

p1 = plot(t,squeeze(rms(h1,2)),'linewidth',2,'color',[250 164 25]./255)
hold on
p2 = plot(t,squeeze(rms(h2,2)),'linewidth',2,'color',[46 75 160]./255)
legend([p1 p2],'direct sound','echo');
set(gca,'fontname','arial')
set(gca,'fontsize',13);
set(gca,'ytick',[0:0.001:0.003])
set(gcf,'unit','centimeters','position',[1.50 1.5 10.5 8])
xlabel('Time (s)')
xlim([-0.2 0.8])
ylim([0 2.1e-3])
title('Exp1 (TRF_{Streaming)}')
% sgtitle('Pilot experiment concat TRF')

%% topography
addpath '.\tools'
clear p1
% echo
pp1_tidu = sqrt((pp_sum1(2:3:306,:).^2+pp_sum1(3:3:306,:).^2)/2);
% first+second
pp2_tidu = sqrt((pp_sum2(2:3:306,:).^2+pp_sum2(3:3:306,:).^2)/2);
% bootstrap
Bv = 10000;

for chan = 1:102
    con = [squeeze(pp1_tidu(chan,:));squeeze(pp2_tidu(chan,:))];
    [O,A]=bootstrap_for_vector(Bv,0.05,@compare_mean,con'); 
    p1(chan)=2*[min(sum(A<=0),sum(A>=0))+1]/[Bv+1];
end

[h,q,p_new] = fdr_bh(p1);
P_value = 0.05;

figure
PKUtopoGradOutline(mean(pp2_tidu(:,:)-pp1_tidu(:,:),2),3);
colormap(flipud(othercolor('RdBu4')))
set(gcf,'unit','centimeters','position',[1.5 1.5 10 8])
colorbar
% title(condition{i})
hold on
PKUtopoChN(find(p_new<P_value),10);
set(gca,'fontsize',13);
set(gca,'fontname','Arial');
title(['Exp1 (topography)'])
clim([-0.02 0.02])


%% Exp2
clear 
clc
fs = 100;
tik = 0.18;
d = fs*1;
lag = 'all';
addpath '.\tools'
i = 1;
shift = 0.2;
for sub = [1:12]
    load(['.\data_availability\TRF\Exp2\idealized_model\tik' num2str(tik) 'd' num2str(d) '\wid_trf_' num2str(sub) '.mat']);

    pp_sum0(:,i) = pp_new2_2;
    h_sum0(:,:,:,i) = h_new2_2;
    i = i+1;
end

tik = 0.06;
i = 1;
shift = 0.2;
for sub = [1:12]
    load(['.\data_availability\TRF\Exp2\baseline_model\tik' num2str(tik) 'd' num2str(d) '\wid_trf_' num2str(sub) '.mat']);
    
    pp_sum1(:,i) = pp_new2_3;
    h_sum1(:,:,:,i) = h_new2_3;
    i = i+1;
end

tik = 0.09;
i = 1;
shift = 0.2;
for sub = [1:12]
    load(['.\data_availability\TRF\Exp2\streaming_model\tik' num2str(tik) 'd' num2str(d) '\wid_trf_' num2str(sub) '.mat']);

    pp_sum2(:,i) = pp_new2_4;
    h_sum2(:,:,:,i) = h_new2_4;
    i = i+1;
end


tidu = setdiff(1:306,1:3:306);
pp_sum = squeeze(cat(3,squeeze(mean(pp_sum1(tidu,:),1)),squeeze(mean(pp_sum2(tidu,:),1)),squeeze(mean(pp_sum0(tidu,:),1))));
err = squeeze(std(pp_sum,0,1)/sqrt(size(pp_sum,1)));
b = bar(squeeze(mean(pp_sum,1)));
hold on
for i = 1:size(err,1)
    xx(:,i) = b(1,i).XEndPoints';
end
b(1).FaceColor = 'none';
b(1).LineWidth = 2;
b(1).BarWidth = 0.8;
set(gca,'fontname','arial')
set(gca,'fontsize',13);
set(gcf,'unit','centimeters','position',[1.50 1.5 8 10])
ylabel('predictive power')
ylim([0 0.1])
set(gca,'xticklabel',{'echo','first+second'})
plot(pp_sum','color',[195 195 195]./255,'linewidth',1)
s1 = scatter(ones(12,1),pp_sum(:,1),'MarkerFaceColor',[250 164 25]./255,'MarkerEdgeColor','none','SizeData',40)
s2 = scatter(ones(12,1)*2,pp_sum(:,2),'MarkerFaceColor',[46 75 160]./255,'MarkerEdgeColor','none','SizeData',40)
s3 = scatter(ones(12,1)*3,pp_sum(:,3),'MarkerFaceColor',[46 75 160]./255,'MarkerEdgeColor','none','SizeData',40)


e = errorbar(xx,squeeze(mean(pp_sum,1)),err,'linestyle','none','linewidth',2,'color','k')
e.CapSize = 10;

set(gca,'xticklabel',{'baseline','streaming','idealized'})
set(gca,'fontname','arial')
set(gca,'fontsize',15);
set(gcf,'unit','centimeters','position',[1.50 1.5 15 10])
box off


addpath 'E:\echo_gjx\others\code\tools'
Bv = 10000;
[O,A]=bootstrap_for_vector(Bv,0.05,@compare_mean,pp_sum(:,[1 2])); 
p=[sum(A>=0)+1]/[Bv+1];
p1(1)=2*[min(sum(A<=0),sum(A>=0))+1]/[Bv+1];

[O,A]=bootstrap_for_vector(Bv,0.05,@compare_mean,pp_sum(:,[3 2])); 
p=[sum(A>=0)+1]/[Bv+1];
p1(2)=2*[min(sum(A<=0),sum(A>=0))+1]/[Bv+1];
[h,q,p_new] = fdr_bh(p1)

%%

fs = 100;
t = -20:1:79;t=t/fs;

figure
sub = [1:12];
% subplot(1,2,2)
h1 = squeeze(rms(h_sum2(:,1,tidu,sub),3))-mean(squeeze(rms(h_sum2(1:21,1,tidu,sub),3)));
err_1 = std(h1,0,2)./sqrt(size(h1,2));
h2 = squeeze(rms(h_sum2(:,2,tidu,sub),3))-mean(squeeze(rms(h_sum2(1:21,2,tidu,sub),3)));
err_2 = std(h1,0,2)./sqrt(size(h2,2));
f1 = fill([t fliplr(t)],[rms(h1,2)+err_1;flipud(rms(h1,2)-err_1)],[195 195 195]./255);
hold on
f1.EdgeColor = 'none';
f1.FaceColor = [250 164 25]./255;
f1.FaceAlpha = 0.3;
f2 = fill([t fliplr(t)],[rms(h2,2)+err_2;flipud(rms(h2,2)-err_2)],[195 195 195]./255);
hold on
f2.EdgeColor = 'none';
f2.FaceColor = [46 75 160]./255;
f2.FaceAlpha = 0.3;

p1 = plot(t,squeeze(rms(h1,2)),'linewidth',2,'color',[250 164 25]./255)
hold on
p2 = plot(t,squeeze(rms(h2,2)),'linewidth',2,'color',[46 75 160]./255)
legend([p1 p2],'direct sound','echo');

set(gca,'fontname','arial')
set(gca,'fontsize',13);
set(gca,'ytick',[0:0.001:0.003])
set(gcf,'unit','centimeters','position',[1.50 1.5 10.5 8])
xlabel('Time (s)')
xlim([-0.2 0.8])
ylim([0 2.1e-3])
title('Exp2 (TRF_{Streaming)}')

%% topography
% echo
addpath '.\tools'
clear p1 p2
pp1_tidu = sqrt((pp_sum1(2:3:306,:).^2+pp_sum1(3:3:306,:).^2)/2);
% first+second
pp2_tidu = sqrt((pp_sum2(2:3:306,:).^2+pp_sum2(3:3:306,:).^2)/2);
% bootstrap
Bv = 10000;

for chan = 1:102
    con = [squeeze(pp1_tidu(chan,:));squeeze(pp2_tidu(chan,:))];
    [O,A]=bootstrap_for_vector(Bv,0.05,@compare_mean,con'); 
    p1(chan)=2*[min(sum(A<=0),sum(A>=0))+1]/[Bv+1];
end

[h,q,p_new] = fdr_bh(p1);
P_value = 0.05;

figure
PKUtopoGradOutline(mean(pp2_tidu(:,:)-pp1_tidu(:,:),2),3);
colormap(flipud(othercolor('RdBu4')))
set(gcf,'unit','centimeters','position',[1.5 1.5 10 8])
colorbar
% title(condition{i})
hold on
PKUtopoChN(find(p_new<P_value),10);
set(gca,'fontsize',13);
set(gca,'fontname','Arial');
title(['Exp2 (topography)'])
clim([-0.02 0.02])
