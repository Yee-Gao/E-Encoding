% adaptation&TRF

%% exp1
clear
clc
close all
fs = 100;
d = fs*1;
lag = 'all';
shift = 0.2;
i = 1;

load(['.\data_availability\TRF\Ada+TRF\Exp1\ada_idealized_model\wid_trf.mat']);
pp_sum0 = squeeze(pp_new2_2_sum);
load(['.\data_availability\TRF\Ada+TRF\Exp1\ada_mixture_model\wid_trf.mat']);
pp_sum1 = squeeze(pp_new2_3_sum);
load(['.\data_availability\TRF\Ada+TRF\Exp1\ada_streaming_model\wid_trf.mat']);
pp_sum2 = squeeze(pp_new2_4_sum);

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


%% significant test (exp 1)

randddd = dec2bin(0:(2^15-1), 15);
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
p1_dep(1) = (length(find(chance_sb>mean(pp_sum(:,2)-pp_sum(:,1),1)))+1)/(length(chance_sb)+1);

chance_sb = squeeze(mean(pp_chance_si(:,2,:)-pp_chance_si(:,1,:),1));
p1_dep(2) = (length(find(chance_sb>mean(pp_sum(:,2)-pp_sum(:,3),1)))+1)/(length(chance_sb)+1);

[~,~,p_new] = fdr_bh(p1_dep)


%% exp2

clear
clc
fs = 100;
d = fs*1;
lag = 'all';
shift = 0.2;


load(['.\data_availability\TRF\Ada+TRF\Exp2\ada_idealized_model\wid_trf.mat']);
pp_sum0 = pp_new2_2_sum(:,1:12);
load(['.\data_availability\TRF\Ada+TRF\Exp2\ada_mixture_model\wid_trf.mat']);
pp_sum1 = pp_new2_3_sum(:,1:12);
load(['.\data_availability\TRF\Ada+TRF\Exp2\ada_streaming_model\wid_trf.mat']);
pp_sum2 = pp_new2_4_sum(:,1:12);

figure
tidu = setdiff(1:306,1:3:306);
pp_sum = squeeze(cat(3,squeeze(mean(pp_sum1(tidu,:),1)),squeeze(mean(pp_sum2(tidu,:),1)),squeeze(mean(pp_sum0(tidu,:),1))));
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
%% significant test (exp2)

randddd = dec2bin(0:(2^12-1), 12);
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
    label = rd(i,:);
    pos = find(label == 1);
    pp_chance = pp_sum(:,[2 3]);
    pp_chance(pos,:) = pp_sum(pos,[3 2]);
    pp_chance_si(:,:,i) = pp_chance;
    clear pp_chance
end
chance_sb = squeeze(mean(pp_chance_sb(:,2,:)-pp_chance_sb(:,1,:),1));
p1_dep(1) = (length(find(chance_sb>mean(pp_sum(:,2)-pp_sum(:,1),1)))+1)/(length(chance_sb)+1);

chance_sb = squeeze(mean(pp_chance_si(:,2,:)-pp_chance_si(:,1,:),1));
p1_dep(2) = (length(find(chance_sb>mean(pp_sum(:,2)-pp_sum(:,3),1)))+1)/(length(chance_sb)+1);

[~,~,p_new] = fdr_bh(p1_dep)

%% exp3

clear
clc
fs = 100;
tik = 0.1;
d = fs*1;
lag = 'all';
i = 1;

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
set(gca,'ytick',[-0.02:0.02:0.1])
ylim([0 0.1])
box off
end


%% significant test
randddd = dec2bin(0:(2^21-1), 21);
rd = double(randddd)-48;
parfor i = 1:size(rd,1)
    label = rd(i,:);
    pos = find(label == 1);
    pp_chance_nv = squeeze(pp_sum(1,:,[1 2]));
    pp_chance_nv(pos,:) = squeeze(pp_sum(1,pos,[2 1]));
    pp_chance_nv_sb(:,:,i) = pp_chance_nv;
end

parfor i = 1:size(rd,1)
    label = rd(i,:);
    pos = find(label == 1);
    pp_chance_nv = squeeze(pp_sum(1,:,[2 3]));
    pp_chance_nv(pos,:) = squeeze(pp_sum(1,pos,[3 2]));
    pp_chance_nv_si(:,:,i) = pp_chance_nv;
end
chance_sb = squeeze(mean(pp_chance_nv_sb(:,2,:)-pp_chance_nv_sb(:,1,:),1));
p_dep(1,1) = (length(find(chance_sb>mean(pp_sum(1,:,2)-pp_sum(1,:,1),2)))+1)/(length(chance_sb)+1);

chance_sb = squeeze(mean(pp_chance_nv_si(:,2,:)-pp_chance_nv_si(:,1,:),1));
p_dep(1,2) = (length(find(chance_sb>mean(pp_sum(1,:,2)-pp_sum(1,:,3),2)))+1)/(length(chance_sb)+1);
[~,~,p_new] = fdr_bh(p_dep(1,:))

%% unattend
randddd = dec2bin(0:(2^21-1), 21);
rd = double(randddd)-48;
parfor i = 1:size(rd,1)
    label = rd(i,:);
    pos = find(label == 1);
    pp_chance_nv = squeeze(pp_sum(2,:,[1 2]));
    pp_chance_nv(pos,:) = squeeze(pp_sum(2,pos,[2 1]));
    pp_chance_nv_sb(:,:,i) = pp_chance_nv;
end

parfor i = 1:size(rd,1)
    label = rd(i,:);
    pos = find(label == 1);
    pp_chance_nv = squeeze(pp_sum(2,:,[2 3]));
    pp_chance_nv(pos,:) = squeeze(pp_sum(2,pos,[3 2]));
    pp_chance_nv_si(:,:,i) = pp_chance_nv;
end
chance_sb = squeeze(mean(pp_chance_nv_sb(:,2,:)-pp_chance_nv_sb(:,1,:),1));
p_dep(2,1) = (length(find(chance_sb>mean(pp_sum(2,:,2)-pp_sum(2,:,1),2)))+1)/(length(chance_sb)+1);

chance_sb = squeeze(mean(pp_chance_nv_si(:,2,:)-pp_chance_nv_si(:,1,:),1));
p_dep(2,2) = (length(find(chance_sb>mean(pp_sum(2,:,2)-pp_sum(2,:,3),2)))+1)/(length(chance_sb)+1);
[~,~,p_new] = fdr_bh(p_dep(2,:))

%% vocode
parfor i = 1:size(rd,1)
    label = rd(i,:);
    pos = find(label == 1);
    pp_chance_v = squeeze(pp_sum(3,:,[1 2]));
    pp_chance_v(pos,:) = squeeze(pp_sum(3,pos,[2 1]));
    pp_chance_v_sb(:,:,i) = pp_chance_v;
end

parfor i = 1:size(rd,1)
    label = rd(i,:);
    pos = find(label == 1);
    pp_chance_v = squeeze(pp_sum(3,:,[1 3]));
    pp_chance_v(pos,:) = squeeze(pp_sum(3,pos,[3 1]));
    pp_chance_v_bi(:,:,i) = pp_chance_v;
end
chance_sb = squeeze(mean(pp_chance_v_sb(:,2,:)-pp_chance_v_sb(:,1,:),1));
p_dep(3,1) = (length(find(chance_sb>mean(pp_sum(3,:,1)-pp_sum(3,:,2),2)))+1)/(length(chance_sb)+1);

chance_sb = squeeze(mean(pp_chance_v_bi(:,2,:)-pp_chance_v_bi(:,1,:),1));
p_dep(3,2) = (length(find(chance_sb>mean(pp_sum(3,:,1)-pp_sum(3,:,3),2)))+1)/(length(chance_sb)+1);
[~,~,p_new] = fdr_bh(p_dep(3,:))



