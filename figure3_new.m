clear
clc

load('.\data_availability\phase_coherence\Exp1\phy_coh_clean.mat');
phy_coh1 = phy_coh_clean;

load('.\data_availability\phase_coherence\Exp2\phy_coh_clean.mat');
phy_coh2 = phy_coh_clean;

fs = 100;
f = 1:size(phy_coh2,1);f=f-1;f=f/size(phy_coh1,1)*fs;
Bv = 10000;
figure
i = 1;
for frebin = [1 2 4 6]
fp = find(f == frebin);
tidu = setdiff(1:306,1:3:306);

e1_submean(:,1) = squeeze(mean(mean(phy_coh1(fp,tidu,2,:),2),4))./squeeze(mean(mean(phy_coh1(fp,tidu,1,:),2),4));
e1_submean(:,2) = squeeze(mean(mean(phy_coh1(fp,tidu,3,:),2),4))./squeeze(mean(mean(phy_coh1(fp,tidu,1,:),2),4));

e1(:,1) = 10*log10(squeeze(mean(phy_coh1(fp,tidu,2,:),2))./squeeze(mean(phy_coh1(fp,tidu,1,:),2)));
e1(:,2) = 10*log10(squeeze(mean(phy_coh1(fp,tidu,3,:),2))./squeeze(mean(phy_coh1(fp,tidu,1,:),2)));
% e1 = e1_submean;
% err1 = squeeze(mean(phy_coh1(fp,tidu,2,:),2),4))./squeeze(mean(mean(phy_coh1(fp,tidu,1,:),2),4))./sqrt(15);
err1 = std(e1,0,1)./sqrt(15);

e2_submean(:,1) = squeeze(mean(mean(phy_coh2(fp,tidu,2,:),2),4))./squeeze(mean(mean(phy_coh2(fp,tidu,1,:),2),4));
e2_submean(:,2) = squeeze(mean(mean(phy_coh2(fp,tidu,3,:),2),4))./squeeze(mean(mean(phy_coh2(fp,tidu,1,:),2),4));
% e2 = e2_submean;% err2 = exp(sqrt(sum(log(e2(:,1)./geomean(e2(:,1))).^2)/14)/sqrt(15));

e2(:,1) = 10*log10(squeeze(mean(phy_coh2(fp,tidu,2,:),2))./squeeze(mean(phy_coh2(fp,tidu,1,:),2)));
e2(:,2) = 10*log10(squeeze(mean(phy_coh2(fp,tidu,3,:),2))./squeeze(mean(phy_coh2(fp,tidu,1,:),2)));
err2 = std(e2,0,1)./sqrt(12);

s_sum =cat(2,squeeze(mean(e1,1)),squeeze(mean(e2,1)));

err = [err1,err2];

subplot(4,1,i)
b = bar(s_sum);
% b.BaseValue = 1;
hold on
e = errorbar([1:4],s_sum,err,'linestyle','none','LineWidth',2,'CapSize',10);
e.Color = 'k';
e.LineWidth = 2;
b.FaceColor = 'flat';
b.CData(2,:) = [255 153 153]./255;
b.CData(4,:) = [255 153 153]./255;
b.CData(1,:) = [219 50 50]./255;
b.CData(3,:) = [219 50 50]./255;
b.BarWidth = 0.8;
b.LineWidth = 2;
i = i+1;
ylim([-4 4])
ylabel('coherence')
title([num2str(frebin) 'Hz'])
set(gca,'xticklabel',{'125','250','125','250'})
set(gca,'fontsize',13)
set(gca,'ytick',[-4 0 4])

end
set(gcf,'unit','centimeters','position',[0 0 10 25])
sgtitle('MEG response')

%% bootstrap test for MEG response

fs = 100;
f = 1:size(phy_coh2,1);f=f-1;f=f/size(phy_coh1,1)*fs;
addpath 'E:\echo_gjx\others\code\tools'
Bv = 10000;

i = 1;
for frebin = [1 2 4 6]
    fp = find(f == frebin);
    tidu = setdiff(1:306,1:3:306);
    
    e1(:,1) = 10*log10(squeeze(mean(phy_coh1(fp,tidu,2,:),2))./squeeze(mean(phy_coh1(fp,tidu,1,:),2)));
    e1(:,2) = 10*log10(squeeze(mean(phy_coh1(fp,tidu,3,:),2))./squeeze(mean(phy_coh1(fp,tidu,1,:),2)));
    e2(:,1) = 10*log10(squeeze(mean(phy_coh2(fp,tidu,2,:),2))./squeeze(mean(phy_coh2(fp,tidu,1,:),2)));
    e2(:,2) = 10*log10(squeeze(mean(phy_coh2(fp,tidu,3,:),2))./squeeze(mean(phy_coh2(fp,tidu,1,:),2)));
    
    [O,A]=bootstrap_for_vector(Bv,0.05,@compare_mean,e1); 
    % A= A-mean(A);
    % p1(1)=2*[min(sum(A<=O),sum(A>=O))+1]/[Bv+1];
    p_pair(1,i)=2*[min(sum(A<=0),sum(A>=0))+1]/[Bv+1];

    [O,A]=bootstrap_for_vector(Bv,0.05,@compare_mean,e2); 
    % A= A-mean(A);
    % p1(1)=2*[min(sum(A<=O),sum(A>=O))+1]/[Bv+1];
    p_pair(2,i)=2*[min(sum(A<=0),sum(A>=0))+1]/[Bv+1];
    i = i+1;
%     p_pair
end
%% unpaired
i = 1;
for frebin = [1 2 4 6]
    fp = find(f == frebin);
    tidu = setdiff(1:306,1:3:306);
    
    e1(:,1) = 10*log10(squeeze(mean(phy_coh1(fp,tidu,2,:),2))./squeeze(mean(phy_coh1(fp,tidu,1,:),2)));
    e1(:,2) = 10*log10(squeeze(mean(phy_coh1(fp,tidu,3,:),2))./squeeze(mean(phy_coh1(fp,tidu,1,:),2)));
    e2(:,1) = 10*log10(squeeze(mean(phy_coh2(fp,tidu,2,:),2))./squeeze(mean(phy_coh2(fp,tidu,1,:),2)));
    e2(:,2) = 10*log10(squeeze(mean(phy_coh2(fp,tidu,3,:),2))./squeeze(mean(phy_coh2(fp,tidu,1,:),2)));
    
    [O,A_125]=bootstrap_for_vector(Bv,0.05,@compare_mean,[e1(:,1),zeros(15,1)]); 
    [O,A_250]=bootstrap_for_vector(Bv,0.05,@compare_mean,[e1(:,2),zeros(15,1)]); 
    % A= A-mean(A);
    % p1(1)=2*[min(sum(A<=O),sum(A>=O))+1]/[Bv+1];
    p_unpair(1,i)=2*[min(sum(mean(e2(:,1))<=A_125),sum(mean(e2(:,1))>=A_125))+1]/[Bv+1];
    p_unpair(2,i)=2*[min(sum(mean(e2(:,2))<=A_125),sum(mean(e2(:,2))>=A_125))+1]/[Bv+1];

    p_unpair(3,i)=2*[min(sum(mean(e2(:,1))<=A_250),sum(mean(e2(:,1))>=A_250))+1]/[Bv+1];
    p_unpair(4,i)=2*[min(sum(mean(e2(:,2))<=A_250),sum(mean(e2(:,2))>=A_250))+1]/[Bv+1];
%     p_unpair
    i = i+1;
end
%% base 1 sig
i = 1;
for frebin = [1 2 4 6]
    fp = find(f == frebin);
    tidu = setdiff(1:306,1:3:306);
    
    e1(:,1) = 10*log10(squeeze(mean(phy_coh1(fp,tidu,2,:),2))./squeeze(mean(phy_coh1(fp,tidu,1,:),2)));
    e1(:,2) = 10*log10(squeeze(mean(phy_coh1(fp,tidu,3,:),2))./squeeze(mean(phy_coh1(fp,tidu,1,:),2)));
    e2(:,1) = 10*log10(squeeze(mean(phy_coh2(fp,tidu,2,:),2))./squeeze(mean(phy_coh2(fp,tidu,1,:),2)));
    e2(:,2) = 10*log10(squeeze(mean(phy_coh2(fp,tidu,3,:),2))./squeeze(mean(phy_coh2(fp,tidu,1,:),2)));
    
    [O,A]=bootstrap_for_vector(Bv,0.05,@compare_mean,[e1(:,1),zeros(15,1)]); 
    % A= A-mean(A);
    % p1(1)=2*[min(sum(A<=O),sum(A>=O))+1]/[Bv+1];
    p_base(1,i)=2*[min(sum(A<=0),sum(A>=0))+1]/[Bv+1];

    [O,A]=bootstrap_for_vector(Bv,0.05,@compare_mean,[e1(:,2),zeros(15,1)]); 
    % A= A-mean(A);
    % p1(1)=2*[min(sum(A<=O),sum(A>=O))+1]/[Bv+1];
    p_base(2,i)=2*[min(sum(A<=0),sum(A>=0))+1]/[Bv+1];



    [O,A]=bootstrap_for_vector(Bv,0.05,@compare_mean,[e2(:,1),zeros(12,1)]); 
    % A= A-mean(A);
    % p1(1)=2*[min(sum(A<=O),sum(A>=O))+1]/[Bv+1];
    p_base(3,i)=2*[min(sum(A<=0),sum(A>=0))+1]/[Bv+1];

    [O,A]=bootstrap_for_vector(Bv,0.05,@compare_mean,[e2(:,2),zeros(12,1)]); 
    % A= A-mean(A);
    % p1(1)=2*[min(sum(A<=O),sum(A>=O))+1]/[Bv+1];
    p_base(4,i)=2*[min(sum(A<=0),sum(A>=0))+1]/[Bv+1];
    i = i+1;
%     p_base

%     p_pair
end
p_sum = [p_base; p_pair; p_unpair];

[h,q,p1_new_sum] = fdr_bh(p_sum(:,1));
[h,q,p2_new_sum] = fdr_bh(p_sum(:,2));
[h,q,p3_new_sum] = fdr_bh(p_sum(:,3));
[h,q,p4_new_sum] = fdr_bh(p_sum(:,4));

%% adaptation model

clear
clc

fs = 100;

load('.\data_availability\Stimuli\anechoic_env.mat')
sig_env = [];


for i = 1:size(sig_speech_sum,2)
    sig_env = [sig_env;sig_speech_sum{i}];% original anechoic speech
end

for  tao_sd =  [0]
for V_th = [0.15]
for  tao_gn =  [0.3]
load(['.\data_availability\adaptation_model\model1_mes\rsdgn_anechoic.mat'])

r_sdgn1 = r_sdgn;
load(['.\data_availability\adaptation_model\model1_mes\exp2\rsdgn_125.mat'])

r_sdgn2 = r_sdgn;
load(['.\data_availability\adaptation_model\model1_mes\exp2\rsdgn_250.mat'])
r_sdgn3 = r_sdgn;
end
end
end

fs = 100;
sig1 = sig_env;
sig2 = r_sdgn1;

for i = 1:1000
    try x = sig1((i-1)*2*fs+1:i*fs*2);
    catch disp(i);
        break
    end
    x_sum(:,i) = x;

    try y = sig2((i-1)*2*fs+1:i*fs*2);
    catch disp(i);
        break
    end
    y_sum(:,i) = y; 
end
fre_x = fft(x_sum);
fre_y = fft(y_sum);
phy = angle(fre_x)-angle(fre_y);
phy_coh = mean(sin(phy),2).^2+mean(cos(phy),2).^2;
phy_coh_sum2(:,1) = phy_coh;

sig1 = sig_env;
sig2 = r_sdgn2;
clear x_sum y_sum
for i = 1:1000
    try x = sig1((i-1)*2*fs+1:i*fs*2);
    catch disp(i);
        break
    end
    x_sum(:,i) = x;

    try y = sig2((i-1)*2*fs+1:i*fs*2);
    catch disp(i);
        break
    end
    y_sum(:,i) = y; 
end
fre_x = fft(x_sum);
fre_y = fft(y_sum);
phy = angle(fre_x)-angle(fre_y);
phy_coh = mean(sin(phy),2).^2+mean(cos(phy),2).^2;
phy_coh_sum2(:,2) = phy_coh;

sig1 = sig_env;
sig2 = r_sdgn3;
clear x_sum y_sum

for i = 1:1000
    try x = sig1((i-1)*2*fs+1:i*fs*2);
    catch disp(i);
        break
    end
    x_sum(:,i) = x;

    try y = sig2((i-1)*2*fs+1:i*fs*2);
    catch disp(i);
        break
    end
    y_sum(:,i) = y; 
end
fre_x = fft(x_sum);
fre_y = fft(y_sum);
phy = angle(fre_x)-angle(fre_y);
phy_coh = mean(sin(phy),2).^2+mean(cos(phy),2).^2;
phy_coh_sum2(:,3) = phy_coh;


for  tao_sd =  [0]
for V_th = [0.15]
for  tao_gn =  [0.3]
load(['.\data_availability\adaptation_model\model1_mes\rsdgn_anechoic.mat'])

r_sdgn1 = r_sdgn;
load(['.\data_availability\adaptation_model\model1_mes\exp1\rsdgn_125.mat'])
r_sdgn2 = r_sdgn;
load(['.\data_availability\adaptation_model\model1_mes\exp1\rsdgn_250.mat'])
r_sdgn3 = r_sdgn;
end
end
end

fs = 100;
sig1 = sig_env;
sig2 = r_sdgn1;

for i = 1:1000
    try x = sig1((i-1)*2*fs+1:i*fs*2);
    catch disp(i);
        break
    end
    x_sum(:,i) = x;

    try y = sig2((i-1)*2*fs+1:i*fs*2);
    catch disp(i);
        break
    end
    y_sum(:,i) = y; 
end
fre_x = fft(x_sum);
fre_y = fft(y_sum);
phy = angle(fre_x)-angle(fre_y);
phy_coh = mean(sin(phy),2).^2+mean(cos(phy),2).^2;
phy_coh_sum1(:,1) = phy_coh;

sig1 = sig_env;
sig2 = r_sdgn2;
clear x_sum y_sum
for i = 1:1000
    try x = sig1((i-1)*2*fs+1:i*fs*2);
    catch disp(i);
        break
    end
    x_sum(:,i) = x;

    try y = sig2((i-1)*2*fs+1:i*fs*2);
    catch disp(i);
        break
    end
    y_sum(:,i) = y; 
end
fre_x = fft(x_sum);
fre_y = fft(y_sum);
phy = angle(fre_x)-angle(fre_y);
phy_coh = mean(sin(phy),2).^2+mean(cos(phy),2).^2;
phy_coh_sum1(:,2) = phy_coh;

sig1 = sig_env;
sig2 = r_sdgn3;
clear x_sum y_sum

for i = 1:1000
    try x = sig1((i-1)*2*fs+1:i*fs*2);
    catch disp(i);
        break
    end
    x_sum(:,i) = x;

    try y = sig2((i-1)*2*fs+1:i*fs*2);
    catch disp(i);
        break
    end
    y_sum(:,i) = y; 
end
fre_x = fft(x_sum);
fre_y = fft(y_sum);
phy = angle(fre_x)-angle(fre_y);
phy_coh = mean(sin(phy),2).^2+mean(cos(phy),2).^2;
phy_coh_sum1(:,3) = phy_coh;

figure
fs = 100;
f = 1:size(phy_coh_sum2,1);f=f-1;f=f/size(phy_coh_sum1,1)*fs;
i = 1;

for frebin = [1 2 4 6]
fp = find(f == frebin);
tidu = setdiff(1:306,1:3:306);

e1(:,1) = 10*log10(squeeze(phy_coh_sum1(fp,2)./phy_coh_sum1(fp,1)));
e1(:,2) = 10*log10(squeeze(phy_coh_sum1(fp,3)./phy_coh_sum1(fp,1)));
e2(:,1) = 10*log10(squeeze(phy_coh_sum2(fp,2)./phy_coh_sum2(fp,1)));
e2(:,2) = 10*log10(squeeze(phy_coh_sum2(fp,3)./phy_coh_sum2(fp,1)));

s_sum =cat(2,squeeze(mean(e1,1)),squeeze(mean(e2,1)));
subplot(4,1,i)
b = bar(s_sum);
b.FaceColor = 'flat';
b.CData(1,:) = [0 102 204]./255;
b.CData(3,:) = [0 102 204]./255;
b.CData(2,:) = [153 204 255]./255;
b.CData(4,:) = [153 204 255]./255;
b.BarWidth = 0.8;
b.LineWidth = 2;
i = i+1;
ylim([-4 4])
ylabel('coherence')
title([num2str(frebin) 'Hz'])
set(gca,'xticklabel',{'125','250','125','250'})
set(gca,'fontsize',13)
set(gca,'ytick',[-4 0 4])

end
set(gcf,'unit','centimeters','position',[0 0 10 25])
set(gca,'ytick',[-4 0 4])
sgtitle('Neural adaptation model')

%% Stimulus

clear
clc

fs = 100;

load('.\data_availability\Stimuli\anechoic_env.mat');
load('.\data_availability\Stimuli\125_exp2_env.mat');
load('.\data_availability\Stimuli\250_exp2_env.mat');


sig_env_echo = [];
sig_env_echo2 = [];
sig_env = [];

for i = 1:size(sig_speech_echo_sum,2)
    sig_env_echo = [sig_env_echo;sig_speech_echo_sum{i}]; % echoic speech
    sig_env_echo2 = [sig_env_echo2;sig_speech_echo_sum2{i}]; % echoic speech
    sig_env = [sig_env;sig_speech_sum{i}];% original anechoic speech
end



fs = 100;
sig1 = sig_env;
sig2 = sig_env;

for i = 1:1000
    try x = sig1((i-1)*2*fs+1:i*fs*2);
    catch disp(i);
        break
    end
    x_sum(:,i) = x;

    try y = sig2((i-1)*2*fs+1:i*fs*2);
    catch disp(i);
        break
    end
    y_sum(:,i) = y; 
end
fre_x = fft(x_sum);
fre_y = fft(y_sum);
phy = angle(fre_x)-angle(fre_y);
phy_coh = mean(sin(phy),2).^2+mean(cos(phy),2).^2;
phy_coh_sum2(:,1) = phy_coh;

sig1 = sig_env;
sig2 = sig_env_echo;
clear x_sum y_sum
for i = 1:1000
    try x = sig1((i-1)*2*fs+1:i*fs*2);
    catch disp(i);
        break
    end
    x_sum(:,i) = x;

    try y = sig2((i-1)*2*fs+1:i*fs*2);
    catch disp(i);
        break
    end
    y_sum(:,i) = y; 
end
fre_x = fft(x_sum);
fre_y = fft(y_sum);
phy = angle(fre_x)-angle(fre_y);
phy_coh = mean(sin(phy),2).^2+mean(cos(phy),2).^2;
phy_coh_sum2(:,2) = phy_coh;

sig1 = sig_env;
sig2 = sig_env_echo2;
clear x_sum y_sum

for i = 1:1000
    try x = sig1((i-1)*2*fs+1:i*fs*2);
    catch disp(i);
        break
    end
    x_sum(:,i) = x;

    try y = sig2((i-1)*2*fs+1:i*fs*2);
    catch disp(i);
        break
    end
    y_sum(:,i) = y; 
end
fre_x = fft(x_sum);
fre_y = fft(y_sum);
phy = angle(fre_x)-angle(fre_y);
phy_coh = mean(sin(phy),2).^2+mean(cos(phy),2).^2;
phy_coh_sum2(:,3) = phy_coh;



load('.\data_availability\Stimuli\anechoic_env.mat');
load('.\data_availability\Stimuli\125_exp1_env.mat');
load('.\data_availability\Stimuli\250_exp1_env.mat');


sig_env_echo = [];
sig_env_echo2 = [];
sig_env = [];

for i = 1:size(sig_speech_echo_sum,2)
    sig_env_echo = [sig_env_echo;sig_speech_echo_sum{i}]; % echoic speech
    sig_env_echo2 = [sig_env_echo2;sig_speech_echo_sum2{i}]; % echoic speech
    sig_env = [sig_env;sig_speech_sum{i}];% original anechoic speech
end



fs = 100;
sig1 = sig_env;
sig2 = sig_env;

for i = 1:1000
    try x = sig1((i-1)*2*fs+1:i*fs*2);
    catch disp(i);
        break
    end
    x_sum(:,i) = x;

    try y = sig2((i-1)*2*fs+1:i*fs*2);
    catch disp(i);
        break
    end
    y_sum(:,i) = y; 
end
fre_x = fft(x_sum);
fre_y = fft(y_sum);
phy = angle(fre_x)-angle(fre_y);
phy_coh = mean(sin(phy),2).^2+mean(cos(phy),2).^2;
phy_coh_sum1(:,1) = phy_coh;

sig1 = sig_env;
sig2 = sig_env_echo;
clear x_sum y_sum
for i = 1:1000
    try x = sig1((i-1)*2*fs+1:i*fs*2);
    catch disp(i);
        break
    end
    x_sum(:,i) = x;

    try y = sig2((i-1)*2*fs+1:i*fs*2);
    catch disp(i);
        break
    end
    y_sum(:,i) = y; 
end
fre_x = fft(x_sum);
fre_y = fft(y_sum);
phy = angle(fre_x)-angle(fre_y);
phy_coh = mean(sin(phy),2).^2+mean(cos(phy),2).^2;
phy_coh_sum1(:,2) = phy_coh;

sig1 = sig_env;
sig2 = sig_env_echo2;
clear x_sum y_sum

for i = 1:1000
    try x = sig1((i-1)*2*fs+1:i*fs*2);
    catch disp(i);
        break
    end
    x_sum(:,i) = x;

    try y = sig2((i-1)*2*fs+1:i*fs*2);
    catch disp(i);
        break
    end
    y_sum(:,i) = y; 
end
fre_x = fft(x_sum);
fre_y = fft(y_sum);
phy = angle(fre_x)-angle(fre_y);
phy_coh = mean(sin(phy),2).^2+mean(cos(phy),2).^2;
phy_coh_sum1(:,3) = phy_coh;

figure
fs = 100;
f = 1:size(phy_coh_sum2,1);f=f-1;f=f/size(phy_coh_sum1,1)*fs;
i = 1;

for frebin = [1 2 4 6]
fp = find(f == frebin);
tidu = setdiff(1:306,1:3:306);

e1(:,1) = 10*log10(squeeze(phy_coh_sum1(fp,2)./phy_coh_sum1(fp,1)));
e1(:,2) = 10*log10(squeeze(phy_coh_sum1(fp,3)./phy_coh_sum1(fp,1)));
e2(:,1) = 10*log10(squeeze(phy_coh_sum2(fp,2)./phy_coh_sum2(fp,1)));
e2(:,2) = 10*log10(squeeze(phy_coh_sum2(fp,3)./phy_coh_sum2(fp,1)));
s_sum =cat(2,squeeze(mean(e1,1)),squeeze(mean(e2,1)));
subplot(4,1,i)
b = bar(s_sum);
b.FaceColor = 'flat';
b.CData(1,:) = [96 96 96]./255;
b.CData(3,:) = [96 96 96]./255;
b.CData(2,:) = [195 195 195]./255;
b.CData(4,:) = [195 195 195]./255;
b.BarWidth = 0.8;
b.LineWidth = 2;
i = i+1;
ylim([-4 4])
ylabel('coherence')
title([num2str(frebin) 'Hz'])
% legend('Experiment 1','Experiment 2')
set(gca,'xticklabel',{'125','250','125','250'})
set(gca,'fontsize',13)
set(gca,'ytick',[-4 0 4])

end
set(gcf,'unit','centimeters','position',[0 0 10 25])
sgtitle('Stimuli')
