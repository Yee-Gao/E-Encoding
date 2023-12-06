%% Fig5
%% exp3 
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
err = squeeze(std(pp_sum,0,2)/sqrt(size(pp_sum,2)));

%% predictive power

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

% attend
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
p_nodep(1,1) = (length(find(chance_sb>mean(pp_sum(1,:,2)-pp_sum(1,:,1),2)))+1)/(length(chance_sb)+1);

chance_sb = squeeze(mean(pp_chance_nv_si(:,2,:)-pp_chance_nv_si(:,1,:),1));
p_nodep(1,2) = (length(find(chance_sb>mean(pp_sum(1,:,2)-pp_sum(1,:,3),2)))+1)/(length(chance_sb)+1);
[~,~,p_new] = fdr_bh(p_nodep(1,:))

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
% save(['first_ProData\pp_chance_nodep_bro\pp_chance_uat.mat'],'pp_chance_nv_sb','pp_chance_nv_si');

chance_sb = squeeze(mean(pp_chance_nv_sb(:,2,:)-pp_chance_nv_sb(:,1,:),1));
p_nodep(2,1) = (length(find(chance_sb>mean(pp_sum(2,:,2)-pp_sum(2,:,1),2)))+1)/(length(chance_sb)+1);

chance_sb = squeeze(mean(pp_chance_nv_si(:,2,:)-pp_chance_nv_si(:,1,:),1));
p_nodep(2,2) = (length(find(chance_sb>mean(pp_sum(2,:,2)-pp_sum(2,:,3),2)))+1)/(length(chance_sb)+1);
[~,~,p_new] = fdr_bh(p_nodep(2,:))

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
p_nodep(3,1) = (length(find(chance_sb>mean(pp_sum(3,:,1)-pp_sum(3,:,2),2)))+1)/(length(chance_sb)+1);

chance_sb = squeeze(mean(pp_chance_v_bi(:,2,:)-pp_chance_v_bi(:,1,:),1));
p_nodep(3,2) = (length(find(chance_sb>mean(pp_sum(3,:,1)-pp_sum(3,:,3),2)))+1)/(length(chance_sb)+1);
[~,~,p_new] = fdr_bh(p_nodep(3,:))


%% TRF
condition = {'Attend to speech','Movie watching','Vocoded speech'};
t = -0.2*fs:1*fs-0.2*fs-1;t=t/fs;
tidu = setdiff(1:306,1:3:306);

for i = 1:3
    if i == 3
    figure
    h1 = squeeze(rms(rms(hh_sum1(:,:,tidu,i,:,:),6),3))-...
        mean(squeeze(rms(rms(hh_sum1(1:20,:,tidu,i,:,:),6),3)));
    err_1 = std(h1,0,2)./sqrt(size(h1,2));
    f = fill([t fliplr(t)],[rms(h1,2)+err_1;flipud(rms(h1,2)-err_1)],[195 195 195]./255);
    hold on
    f.EdgeColor = 'none';
    f.FaceColor = [0 128 0]./255;
    f.FaceAlpha = 0.3;
    p = plot(t,rms(h1,2),'linewidth',2,'color',[0 128 0]./255)

    set(gca,'fontname','arial')
    set(gca,'fontsize',13);
    set(gcf,'unit','centimeters','position',[1.50 1.5 10.5 8])
    xlabel('Time (s)')
    xlim([-0.2 0.8])
    ylim([0 0.003])
    set(gca,'ytick',[0 0.001 0.002 0.003])
    title([condition{i} '_{baseline}'])
    legend([p],'echoic speech')
    else
    figure
    h1 = squeeze(rms(rms(hh_sum2(:,1,tidu,i,:,:),6),3))...
        -mean(squeeze(rms(rms(hh_sum2(1:20,1,tidu,i,:,:),6),3)));
    err_1 = std(h1,0,2)./sqrt(size(h1,2));
    f = fill([t fliplr(t)],[rms(h1,2)+err_1;flipud(rms(h1,2)-err_1)],[195 195 195]./255);
    hold on
    f.EdgeColor = 'none';
    f.FaceColor = [250 164 25]./255;
    f.FaceAlpha = 0.3;
    
    h2 = squeeze(rms(rms(hh_sum2(:,2,tidu,i,:,:),6),3))...
        -mean(squeeze(rms(rms(hh_sum2(1:20,2,tidu,i,:,:),6),3)));
    err_2 = std(h2,0,2)./sqrt(size(h2,2));
    f = fill([t fliplr(t)],[rms(h2,2)+err_2;flipud(rms(h2,2)-err_2)],[195 195 195]./255);
    hold on
    f.EdgeColor = 'none';
    f.FaceColor = [46 75 160]./255;
    f.FaceAlpha = 0.3;

    p1 = plot(t,rms(h1,2),'linewidth',2,'color',[250 164 25]./255)
    hold on
    p2 = plot(t,rms(h2,2),'linewidth',2,'color',[46 75 160]./255)
    legend([p1 p2],'direct sound','echo');
    set(gca,'fontname','arial')
    set(gca,'fontsize',13);
    set(gcf,'unit','centimeters','position',[1.50 1.5 10.5 8])
    xlabel('Time (s)')
    xlim([-0.2 0.8])
    ylim([0 0.002])
    set(gca,'ytick',[0 0.001 0.002 0.003])
    title([condition{i} '_{streaming}'])
    end
end

%% Topography
clear p1 p2
% echo
pp1_tidu = sqrt((pp_sum1(2:3:306,:,:).^2+pp_sum1(3:3:306,:,:).^2)/2);
% first+second
pp2_tidu = sqrt((pp_sum2(2:3:306,:,:).^2+pp_sum2(3:3:306,:,:).^2)/2);
Bv = 10000;
for chan = 1:102
    for i = 1:2
        con = [squeeze(pp1_tidu(chan,i,:)),squeeze(pp2_tidu(chan,i,:))];
        [O,A]=bootstrap_for_vector(Bv,0.05,@compare_mean,con); 
        p1(chan,i)= 2*[min(sum(A<=0),sum(A>=0))+1]/[Bv+1];
    end
    for i = 3
        con = [squeeze(pp1_tidu(chan,i,:)),squeeze(pp2_tidu(chan,i,:))];
        [O,A]=bootstrap_for_vector(Bv,0.05,@compare_mean,con);
        p1(chan,i)= 2*[min(sum(A<=0),sum(A>=0))+1]/[Bv+1];
    end
end

[h,q,p_new] = fdr_bh(p1);
P_value = 0.05;
for i = 1:3
    figure
    PKUtopoGradOutline(mean(pp2_tidu(:,i,:)-pp1_tidu(:,i,:),3),3);
    colormap(flipud(othercolor('RdBu4')))
    set(gcf,'unit','centimeters','position',[1.5 1.5 10 8])
    colorbar
    title(condition{i})
    hold on
    PKUtopoChN(find(p_new(:,i)<P_value),10);
    set(gca,'fontsize',13);
    set(gca,'fontname','Arial');
    title([condition{i}])
    clim([-0.02 0.02])
end
