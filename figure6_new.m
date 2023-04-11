%%

clear
clc
close all
fs = 100;
d = fs*1;
lag = 'all';

tik = 0.24;
i = 1
for sub = [1:9 12 14:24]
%     load(['..\ProData\shift_0.2\trf_new_' lag '\clean\tik' num2str(tik) 'd' num2str(d) '\wid_trf_' num2str(sub) '.mat']);
    load(['.\data_availability\TRF\Exp3\idealized_model\tik' num2str(tik) 'd' num2str(d) '\wid_trf_' num2str(sub) '.mat']);
    pp_sum0(:,1,i) = pp_new2_2(:,1);
    hh_sum0(:,:,:,1,i) = h_new2_2(:,:,:,1);
    i = i+1;
end


tik = 0.26;
i = 1
for sub = [1:9 12 14:24]
%     load(['..\ProData\shift_0.2\trf_new_' lag '\clean\tik' num2str(tik) 'd' num2str(d) '\wid_trf_' num2str(sub) '.mat']);
    load(['.\data_availability\TRF\Exp3\idealized_model\tik' num2str(tik) 'd' num2str(d) '\wid_trf_' num2str(sub) '.mat']);
    pp_sum0(:,2,i) = pp_new2_2(:,2);
    hh_sum0(:,:,:,2,i) = h_new2_2(:,:,:,2);
    i = i+1;
end



tik = inf;
i = 1
for sub = [1:9 12 14:24]
%     load(['..\ProData\shift_0.2\trf_new_' lag '\clean\tik' num2str(tik) 'd' num2str(d) '\wid_trf_' num2str(sub) '.mat']);
    load(['.\data_availability\TRF\Exp3\idealized_model\tik' num2str(tik) 'd' num2str(d) '\wid_trf_' num2str(sub) '.mat']);
    pp_sum0(:,3,i) = pp_new2_2(:,3);
    hh_sum0(:,:,:,3,i) = h_new2_2(:,:,:,3);
    i = i+1;
end


tik = 0.19;
i = 1
for sub = [1:9 12 14:24]
%     load(['..\ProData\shift_0.2\trf_new_' lag '\echo\tik' num2str(tik) 'd' num2str(d) '\wid_trf_' num2str(sub) '.mat']);
    load(['.\data_availability\TRF\Exp3\baseline_model\tik' num2str(tik) 'd' num2str(d) '\wid_trf_' num2str(sub) '.mat']);
    pp_sum1(:,1,i) = pp_new2_3(:,1);
    hh_sum1(:,:,:,1,i) = h_new2_3(:,:,:,1);
    i = i+1;
end

tik = 0.18;
i = 1
for sub = [1:9 12 14:24]
%     load(['..\ProData\shift_0.2\trf_new_' lag '\echo\tik' num2str(tik) 'd' num2str(d) '\wid_trf_' num2str(sub) '.mat']);
    load(['.\data_availability\TRF\Exp3\baseline_model\tik' num2str(tik) 'd' num2str(d) '\wid_trf_' num2str(sub) '.mat']);
    pp_sum1(:,2,i) = pp_new2_3(:,2);
    hh_sum1(:,:,:,2,i) = h_new2_3(:,:,:,2);
    i = i+1;
end

tik = 0.14;
i = 1
for sub = [1:9 12 14:24]
%     load(['..\ProData\shift_0.2\trf_new_' lag '\echo\tik' num2str(tik) 'd' num2str(d) '\wid_trf_' num2str(sub) '.mat']);
    load(['.\data_availability\TRF\Exp3\baseline_model\tik' num2str(tik) 'd' num2str(d) '\wid_trf_' num2str(sub) '.mat']);
    pp_sum1(:,3,i) = pp_new2_3(:,3);
    hh_sum1(:,:,:,3,i) = h_new2_3(:,:,:,3);
    i = i+1;
end

tik = 0.13;
i = 1
for sub = [1:9 12 14:24]

%     load(['..\ProData\shift_0.2\trf_new_' lag '\first+second\tik' num2str(tik) 'd' num2str(d) '\wid_trf_' num2str(sub) '.mat']);
    load(['.\data_availability\TRF\Exp3\streaming_model\tik' num2str(tik) 'd' num2str(d) '\wid_trf_' num2str(sub) '.mat']);
    pp_sum2(:,1,i) = pp_new2_4(:,1);
    hh_sum2(:,:,:,1,i) = h_new2_4(:,:,:,1);
    i = i+1;
end

tik = 0.16;
i = 1
for sub = [1:9 12 14:24]
%     load(['..\ProData\shift_0.2\trf_new_' lag '\first+second\tik' num2str(tik) 'd' num2str(d) '\wid_trf_' num2str(sub) '.mat']);
    load(['.\data_availability\TRF\Exp3\streaming_model\tik' num2str(tik) 'd' num2str(d) '\wid_trf_' num2str(sub) '.mat']);
    pp_sum2(:,2,i) = pp_new2_4(:,2);
    hh_sum2(:,:,:,2,i) = h_new2_4(:,:,:,2);
    i = i+1;
end

tik = 0.23;
i = 1
for sub = [1:9 12 14:24]
%     load(['..\ProData\shift_0.2\trf_new_' lag '\first+second\tik' num2str(tik) 'd' num2str(d) '\wid_trf_' num2str(sub) '.mat']);
    load(['.\data_availability\TRF\Exp3\streaming_model\tik' num2str(tik) 'd' num2str(d) '\wid_trf_' num2str(sub) '.mat']);
    pp_sum2(:,3,i) = pp_new2_4(:,3);
    hh_sum2(:,:,:,3,i) = h_new2_4(:,:,:,3);
    i = i+1;
end

%% baseline streaming and idealized

addpath '.\tools'

tidu = setdiff(1:306,1:3:306);
pp_sum = cat(3,squeeze(mean(pp_sum1(tidu,:,:),1)),squeeze(mean(pp_sum2(tidu,:,:),1)),squeeze(mean(pp_sum0(tidu,:,:),1)));
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
s2 = scatter(2,pp_sum(con,i,2),'MarkerFaceColor',[46 75 160]./255,'MarkerEdgeColor','none','SizeData',40);
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


%% control clean and streaming
addpath '.\tools'
Bv = 10000;
pp_sum = cat(3,squeeze(mean(pp_sum1(tidu,:,:),1)),squeeze(mean(pp_sum2(tidu,:,:),1)),squeeze(mean(pp_sum0(tidu,:,:),1)));

% Attending to speech
[O,A]=bootstrap_for_vector(Bv,0.05,@compare_mean,squeeze(pp_sum(1,:,[1 2]))); 
p1(1)=2*[min(sum(A<=0),sum(A>=0))+1]/[Bv+1];

% A= A-mean(A);
% p1(1)=2*[min(sum(A<=O),sum(A>=O))+1]/[Bv+1];
[O,A]=bootstrap_for_vector(Bv,0.05,@compare_mean,squeeze(pp_sum(1,:,[2 3]))); 
p1(2)=2*[min(sum(A<=0),sum(A>=0))+1]/[Bv+1];

% A= A-mean(A);
% p1(2)=2*[min(sum(A<=O),sum(A>=O))+1]/[Bv+1];
[h,q,p_new] = fdr_bh([p1]);
p_new

% Movie watching
[O,A]=bootstrap_for_vector(Bv,0.05,@compare_mean,squeeze(pp_sum(2,:,[1 2]))); 
% p2=[sum(A<=0)+1]/[Bv+1];
p2(1)=2*[min(sum(A<=0),sum(A>=0))+1]/[Bv+1];
% A= A-mean(A);
% p2(1)=2*[min(sum(A<=O),sum(A>=O))+1]/[Bv+1];
[O,A]=bootstrap_for_vector(Bv,0.05,@compare_mean,squeeze(pp_sum(2,:,[2 3]))); 
% p2=[sum(A<=0)+1]/[Bv+1];
p2(2)=2*[min(sum(A<=0),sum(A>=0))+1]/[Bv+1];
% A= A-mean(A);
% p2(2)=2*[min(sum(A<=O),sum(A>=O))+1]/[Bv+1];
[h,q,p_new] = fdr_bh([p2]);
p_new

% Vocoded speech(Movie watching)
[O,A]=bootstrap_for_vector(Bv,0.05,@compare_mean,squeeze(pp_sum(3,:,[1 2]))); 
p3(1)=2*[min(sum(A<=0),sum(A>=0))+1]/[Bv+1];
% A= A-mean(A);
% p3(1)=2*[min(sum(A<=O),sum(A>=O))+1]/[Bv+1];
[O,A]=bootstrap_for_vector(Bv,0.05,@compare_mean,squeeze(pp_sum(3,:,[1 3]))); 
p3(2)=2*[min(sum(A<=0),sum(A>=0))+1]/[Bv+1];
% A= A-mean(A);
% p3(2)=2*[min(sum(A<=O),sum(A>=O))+1]/[Bv+1];
[h,q,p_new] = fdr_bh([p3]);
p_new


%%
condition = {'Attend to speech','Movie watching','Vocoded speech'};
t = -0.2*fs:1*fs-0.2*fs-1;t=t/fs;
tidu = setdiff(1:306,1:3:306);
addpath '.\tools'

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

%%
clear p1 p2
% echo
pp1_tidu = sqrt((pp_sum1(2:3:306,:,:).^2+pp_sum1(3:3:306,:,:).^2)/2);
% first+second
pp2_tidu = sqrt((pp_sum2(2:3:306,:,:).^2+pp_sum2(3:3:306,:,:).^2)/2);
% bootstrap
Bv = 10000;
addpath '.\tools'

for chan = 1:102
    for i = 1:2
        con = [squeeze(pp1_tidu(chan,i,:)),squeeze(pp2_tidu(chan,i,:))];
        [O,A]=bootstrap_for_vector(Bv,0.05,@compare_mean,con); 
%         A= A-mean(A);
        p1(chan,i)= 2*[min(sum(A<=0),sum(A>=0))+1]/[Bv+1];
%         p1(chan,i)=[sum(A>=0)+1]/[Bv+1];
    end
    for i = 3
        con = [squeeze(pp1_tidu(chan,i,:)),squeeze(pp2_tidu(chan,i,:))];
        [O,A]=bootstrap_for_vector(Bv,0.05,@compare_mean,con);
        p1(chan,i)= 2*[min(sum(A<=0),sum(A>=0))+1]/[Bv+1];
%         p1(chan,i)=[sum(A<=0)+1]/[Bv+1];
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
