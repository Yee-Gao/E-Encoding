%% Fig4

%% Experiment 1
clear
clc
% close all
fs = 100;

load(['.\data_availability\TRF\Exp1\idealized_model\wid_trf.mat']);
pp_sum0 = squeeze(pp_new2_2_sum);
h_sum0 = h_new2_2_sum;
load(['.\data_availability\TRF\Exp1\mixture_model\wid_trf.mat']);
pp_sum1 = squeeze(pp_new2_3_sum);
h_sum1 = h_new2_3_sum;
load(['.\data_availability\TRF\Exp1\streaming_model\wid_trf.mat']);
pp_sum2 = squeeze(pp_new2_4_sum);
h_sum2 = h_new2_4_sum;
%% predictive power
figure
tidu = setdiff(1:306,1:3:306);
pp_sum = squeeze(cat(3,squeeze(mean(pp_sum1(tidu,:),1)),squeeze(mean(pp_sum2(tidu,:),1)),squeeze(mean(pp_sum0(tidu,:),1))));
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
ylim([0 0.1])
plot(pp_sum','color',[195 195 195]./255,'linewidth',1)
s1 = scatter(ones(15,1),pp_sum(:,1),'MarkerFaceColor','#228B22','MarkerEdgeColor','none','SizeData',40)
s2 = scatter(ones(15,1)*2,pp_sum(:,2),'MarkerFaceColor',[46 75 160]./255,'MarkerEdgeColor','none','SizeData',40)
s3 = scatter(ones(15,1)*3,pp_sum(:,3),'MarkerFaceColor',[250 164 25]./255,'MarkerEdgeColor','none','SizeData',40)
errorbar(xx,squeeze(mean(pp_sum,1)),err,'linestyle','none','linewidth',2,'color','k','CapSize',10)
box off       
set(gcf,'Units','centimeters','Position',[1.5 1.5 15 10]);
set(gca,'xticklabel',{'baseline','streaming','idealized'})

pp_nodep = pp_sum;

%% significant test

randddd = dec2bin(0:(2^15-1), 15);
rd = double(randddd)-48;

for i = 1:size(rd,1)
    label = rd(i,:);
    pos = find(label == 1);
    pp_chance = pp_sum(:,[1 2]);
    pp_chance(pos,:) = pp_sum(pos,[2 1]);
    pp_chance_sb(:,:,i) = pp_chance;
    clear pp_chance
end

for i = 1:size(rd,1)
%     tic
    label = rd(i,:);
    pos = find(label == 1);
    pp_chance = pp_sum(:,[2 3]);
    pp_chance(pos,:) = pp_sum(pos,[3 2]);
    pp_chance_si(:,:,i) = pp_chance;
    clear pp_chance
end
chance_sb = squeeze(mean(pp_chance_sb(:,2,:)-pp_chance_sb(:,1,:),1));
p1_nodep(1) = (length(find(chance_sb>mean(pp_sum(:,2)-pp_sum(:,1),1)))+1)/(length(chance_sb)+1);

chance_sb = squeeze(mean(pp_chance_si(:,2,:)-pp_chance_si(:,1,:),1));
p1_nodep(2) = (length(find(chance_sb>mean(pp_sum(:,2)-pp_sum(:,3),1)))+1)/(length(chance_sb)+1);


[~,~,p_new] = fdr_bh(p1_nodep)


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
clear p1 p2
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


%% --------------------------------------------------
% Experiment 2
%%---------------------------------------------------

clear
clc
% close all
fs = 100;

load(['.\data_availability\TRF\Exp2\idealized_model\wid_trf.mat']);
pp_sum0 = pp_new2_2_sum(:,1:12);
h_sum0 = h_new2_2_sum(:,:,:,1:12);
load(['.\data_availability\TRF\Exp2\mixture_model\wid_trf.mat']);
pp_sum1 = pp_new2_3_sum(:,1:12);
h_sum1 = h_new2_3_sum(:,:,:,1:12);
load(['.\data_availability\TRF\Exp2\streaming_model\wid_trf.mat']);
pp_sum2 = pp_new2_4_sum(:,1:12);
h_sum2 = h_new2_4_sum(:,:,:,1:12);

figure
tidu = setdiff(1:306,1:3:306);
pp_sum = squeeze(cat(3,squeeze(mean(pp_sum1(tidu,:),1)),squeeze(mean(pp_sum2(tidu,:),1)),squeeze(mean(pp_sum0(tidu,:),1))));
% pp_sum = squeeze(cat(3,squeeze(mean(pp_sum1(tidu,:),1)),squeeze(mean(pp_sum2(tidu,:),1))));
err = squeeze(std(pp_sum,0,1)/sqrt(size(pp_sum,1)));
b = bar(squeeze(mean(pp_sum,1)));
hold on
for i = 1:size(err,1)
    xx(:,i) = b(1,i).XEndPoints';
end
box off
b(1).FaceColor = 'none';
b(1).LineWidth = 2;
set(gca,'fontname','arial')
set(gca,'fontsize',15);
set(gcf,'unit','centimeters','position',[1.50 1.5 15 10])
ylabel('predictive power')
ylim([0 0.1])
set(gca,'xticklabel',{'baseline','streaming','idealized'})
plot(pp_sum','color',[195 195 195]./255,'linewidth',1)
s1 = scatter(ones(12,1),pp_sum(:,1),'MarkerFaceColor','#228B22','MarkerEdgeColor','none','SizeData',40)
s2 = scatter(ones(12,1)*2,pp_sum(:,2),'MarkerFaceColor',[46 75 160]./255,'MarkerEdgeColor','none','SizeData',40)
s3 = scatter(ones(12,1)*3,pp_sum(:,3),'MarkerFaceColor',[250 164 25]./255,'MarkerEdgeColor','none','SizeData',40)
e = errorbar(xx,squeeze(mean(pp_sum,1)),err,'linestyle','none','linewidth',2,'color','k','CapSize',10)
e.CapSize = 10;

%% significant test

randddd = dec2bin(0:(2^12-1), 12);
rd = double(randddd)-48;

for i = 1:size(rd,1)
%     tic
    label = rd(i,:);
    pos = find(label == 1);
    pp_chance = pp_sum(:,[1 2]);
    pp_chance(pos,:) = pp_sum(pos,[2 1]);
    pp_chance_sb(:,:,i) = pp_chance;
    clear pp_chance
end

for i = 1:size(rd,1)
%     tic
    label = rd(i,:);
    pos = find(label == 1);
    pp_chance = pp_sum(:,[2 3]);
    pp_chance(pos,:) = pp_sum(pos,[3 2]);
    pp_chance_si(:,:,i) = pp_chance;
    clear pp_chance
end
chance_sb = squeeze(mean(pp_chance_sb(:,2,:)-pp_chance_sb(:,1,:),1));
p1_nodep(1) = (length(find(chance_sb>mean(pp_sum(:,2)-pp_sum(:,1),1)))+1)/(length(chance_sb)+1);

chance_sb = squeeze(mean(pp_chance_si(:,2,:)-pp_chance_si(:,1,:),1));
p1_nodep(2) = (length(find(chance_sb>mean(pp_sum(:,2)-pp_sum(:,3),1)))+1)/(length(chance_sb)+1);


[~,~,p_new] = fdr_bh(p1_nodep);



%% TRF

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

