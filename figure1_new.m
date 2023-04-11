%% reverberation cohrence
clear
clc
load('.\data_availability\phase_coherence\coh_reverb.mat');
phy_coh = cohs;
fs = 100;
Color_temp = [[46 75 160]./255;[250 164 25]./255];
f = 1:size(phy_coh,1);f=f-1;f=f/size(phy_coh,1)*fs;
std_rev = std(phy_coh,0,2);

hold on

f1 = fill([f fliplr(f)],[mean(phy_coh,2)+1.96*std_rev;flipud(mean(phy_coh,2)-1.96*std_rev)],[195 195 195]./255);
f1.EdgeColor = 'none';
% f1.FaceColor = 'r';
f1.FaceAlpha = 0.3;
plot(f,mean(phy_coh,2),'color',[195 195 195]./255,'linewidth',2)
xlim([1 14])
ylim([0 1])
xlabel('Frequency(Hz)')
ylabel('Phase coherence')
set(gca,'fontname','arial')
set(gca,'fontsize',13);
set(gcf,'unit','centimeters','position',[1.50 1.5 10.5 8])
title('reverberation')
box on
%% impulse response
clear
clc
load ('.\data_availability\phase_coherence\IR_mat.mat');;
Fs=32000;
cut=63952; %use max_length
IR_mat_abs = abs(IR_mat);
IR_mean = mean(IR_mat_abs, 2);
echo_0125 = zeros(cut, 1);
echo_0125(Fs*0.125+1, 1)=1;
echo_0125(1, 1)=1;

echo_0250 = zeros(cut, 1);
echo_0250(Fs*0.25+1, 1)=1;
echo_0250(1, 1)=1;

IR_mean = IR_mean(1:cut, 1);

x=[];
for i=1:cut
    x = [x, 1/Fs*(i-1)];
end
t = 1:cut;t=t/Fs;

plot(t, IR_mat,'color',[220 220 220]./255);
hold on
plot(t,echo_0125,'LineWidth',2,'color','k');
plot(t,echo_0250,'LineWidth',2,'color',[150 150 150]./255);
ylabel('amplitude')
xlabel('time (s)')

xlim([0 0.5])
set(gca,'ytick',[-1:0.5:1]);
set(gca,'xtick',[0:0.1:0.5]);
set(gca,'fontname','arial')
set(gca,'fontsize',13);
set(gcf,'unit','centimeters','position',[1.50 1.5 10.5 8])
set(gca,'yscale','linear')
box off
title('Impulse response')

%% stimulus phase coherence (0.125s)

clear
clc

fs = 100;

% experiement 1

load(".\data_availability\Stimuli\125_exp1_env.mat");
load(".\data_availability\Stimuli\anechoic_env.mat");

sig_env_echo = [];
sig_env = [];
for i = 1:size(sig_speech_echo_sum,2)
    sig_env_echo = [sig_env_echo;sig_speech_echo_sum{i}]; % echoic speech after depression
    sig_env = [sig_env;sig_speech_sum{i}];% original anechoic speech
end

sig1 = mean(sig_env,2);
sig2 = mean(sig_env_echo,2);

for ind = 1:1000
    try  x = sig1(fs*2*(ind-1)+1:fs*2*ind);
    catch,disp([num2str(ind)]);break;end
    sig_depn(:,ind) = x;
    try  y = sig2(fs*2*(ind-1)+1:fs*2*ind);
    catch,break;end
    sig_envn(:,ind) = y;
end
fre_dep = fft(sig_depn);
fre_env = fft(sig_envn);
phy = angle(fre_dep)-angle(fre_env);
phy_coh = mean(cos(phy),2).^2+mean(sin(phy),2).^2;

f = 1:size(phy_coh,1);f=f-1;f=f/size(phy_coh,1)*fs;
p1 = plot(f,mean(phy_coh,2),'color','k','linewidth',2)
hold on
plot([4 4],[0 1],'color',[176 31 36]./255,'linestyle','--','linewidth',2);

% experiement 2

load(".\data_availability\Stimuli\125_exp2_env.mat");
load(".\data_availability\Stimuli\anechoic_env.mat");

sig_env_echo = [];
sig_env = [];
for i = 1:size(sig_speech_echo_sum,2)
    sig_env_echo = [sig_env_echo;sig_speech_echo_sum{i}]; % echoic speech 
    sig_env = [sig_env;sig_speech_sum{i}];% original anechoic speech
end

sig1 = mean(sig_env,2);
sig2 = mean(sig_env_echo,2);

for ind = 1:1000
    try  x = sig1(fs*2*(ind-1)+1:fs*2*ind);
    catch,disp([num2str(ind)]);break;end
    sig_depn(:,ind) = x;
    try  y = sig2(fs*2*(ind-1)+1:fs*2*ind);
    catch,break;end
    sig_envn(:,ind) = y;
end
fre_dep = fft(sig_depn);
fre_env = fft(sig_envn);
phy = angle(fre_dep)-angle(fre_env);
phy_coh = mean(cos(phy),2).^2+mean(sin(phy),2).^2;

hold on
f = 1:size(phy_coh,1);f=f-1;f=f/size(phy_coh,1)*fs;
p2 = plot(f,mean(phy_coh,2),'color','k','linewidth',2,'linestyle','--');
hold on
plot([4 4],[0 1],'color',[176 31 36]./255,'linestyle','--','linewidth',2);


% reverberation

load('.\data_availability\phase_coherence\coh_reverb.mat');
phy_coh = cohs;
std_rev = std(phy_coh,0,2);
f1 = fill([f fliplr(f)],[mean(phy_coh,2)+1.96*std_rev;flipud(mean(phy_coh,2)-1.96*std_rev)],[195 195 195]./255);
f1.EdgeColor = 'none';

f1.FaceAlpha = 0.3;
p3 = plot(f,mean(phy_coh,2),'color',[195 195 195]./255,'linewidth',2)

xlim([1 14])
ylim([0 1])
xlabel('Frequency(Hz)')
ylabel('Phase coherence')
legend([p1 p2 p3],'Exp1','Exp2','Reverberation');
% set(gca,'xticklabel',{'attend','unattend','vocode'})
set(gca,'fontname','arial')
set(gca,'fontsize',13);
set(gcf,'unit','centimeters','position',[1.50 1.5 13 10])
title('0.125-s echoic speech')


%% stimulus phase coherence (0.25s)

clear
clc

addpath 'E:\echo_gjx\others\code\tools\'
addpath '.\tools'
fs = 100;

% experiement 1

load(".\data_availability\Stimuli\250_exp1_env.mat");
load(".\data_availability\Stimuli\anechoic_env.mat");
sig_env_echo = [];
sig_env = [];
for i = 1:size(sig_speech_echo_sum2,2)
    sig_env_echo = [sig_env_echo;sig_speech_echo_sum2{i}]; % echoic speech after depression
    sig_env = [sig_env;sig_speech_sum{i}];% original anechoic speech
end

sig1 = mean(sig_env,2);
sig2 = mean(sig_env_echo,2);

for ind = 1:1000
    try  x = sig1(fs*2*(ind-1)+1:fs*2*ind);
    catch,disp([num2str(ind)]);break;end
    sig_depn(:,ind) = x;
    try  y = sig2(fs*2*(ind-1)+1:fs*2*ind);
    catch,break;end
    sig_envn(:,ind) = y;
end
fre_dep = fft(sig_depn);
fre_env = fft(sig_envn);
phy = angle(fre_dep)-angle(fre_env);
phy_coh = mean(cos(phy),2).^2+mean(sin(phy),2).^2;

f = 1:size(phy_coh,1);f=f-1;f=f/size(phy_coh,1)*fs;
p1 = plot(f,mean(phy_coh,2),'color','k','linewidth',2)
hold on

% experiement 2

load(".\data_availability\Stimuli\250_exp2_env.mat");
load(".\data_availability\Stimuli\anechoic_env.mat");

sig_env_echo = [];
sig_env = [];
for i = 1:size(sig_speech_echo_sum2,2)
    sig_env_echo = [sig_env_echo;sig_speech_echo_sum2{i}]; % echoic speech
    sig_env = [sig_env;sig_speech_sum{i}];% original anechoic speech
end

sig1 = mean(sig_env,2);
sig2 = mean(sig_env_echo,2);

for ind = 1:1000
    try  x = sig1(fs*2*(ind-1)+1:fs*2*ind);
    catch,disp([num2str(ind)]);break;end
    sig_depn(:,ind) = x;
    try  y = sig2(fs*2*(ind-1)+1:fs*2*ind);
    catch,break;end
    sig_envn(:,ind) = y;
end
fre_dep = fft(sig_depn);
fre_env = fft(sig_envn);
phy = angle(fre_dep)-angle(fre_env);
phy_coh = mean(cos(phy),2).^2+mean(sin(phy),2).^2;


f = 1:size(phy_coh,1);f=f-1;f=f/size(phy_coh,1)*fs;
p2 = plot(f,mean(phy_coh,2),'color','k','linewidth',2,'linestyle','--');
hold on
plot([2 2],[0 1],'color',[176 31 36]./255,'linestyle','--','linewidth',2);
plot([6 6],[0 1],'color',[176 31 36]./255,'linestyle','--','linewidth',2);


% reverberation

load('.\data_availability\phase_coherence\coh_reverb.mat');
phy_coh = cohs;
std_rev = std(phy_coh,0,2);
f1 = fill([f fliplr(f)],[mean(phy_coh,2)+1.96*std_rev;flipud(mean(phy_coh,2)-1.96*std_rev)],[195 195 195]./255);
f1.EdgeColor = 'none';
f1.FaceAlpha = 0.3;
p3 = plot(f,mean(phy_coh,2),'color',[195 195 195]./255,'linewidth',2)

xlim([1 14])
ylim([0 1])
xlabel('Frequency(Hz)')
ylabel('Phase coherence')
legend([p1 p2 p3],'Exp1','Exp2','Reverberation');
% set(gca,'xticklabel',{'attend','unattend','vocode'})
set(gca,'fontname','arial')
set(gca,'fontsize',13);
set(gcf,'unit','centimeters','position',[1.50 1.5 13 10])
title('0.25-s echoic speech')

%% depression and gain control(0.125s)
% 
clear
clc

addpath 'E:\echo_gjx\others\code\tools\'
addpath '.\tools'

tao_sd =  [0];
V_th = [0.15];
tao_gn =  [0.3];
fs = 100;

load(".\data_availability\Stimuli\anechoic_env.mat");
sig_env = [];


for i = 1:size(sig_speech_sum,2)
    sig_env = [sig_env;sig_speech_sum{i}];% original anechoic speech
end
% exp1
load(['.\data_availability\adaptation_model\model1_mes\exp1\rsdgn_125.mat'])
r_sdgn1 = r_sdgn;

% exp2
load(['.\data_availability\adaptation_model\model1_mes\exp2\rsdgn_125.mat'])
r_sdgn2 = r_sdgn;

load(['.\data_availability\adaptation_model\model1_mes\rsdgn_anechoic.mat'])
r_sdgn3 = r_sdgn;

figure
% exp1
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
Color_temp = [[46 75 160]./255;[250 164 25]./255];
f = 1:size(phy_coh,1);f=f-1;f=f/size(phy_coh,1)*fs;
p1 = plot(f,mean(phy_coh,2),'color','k','linewidth',2);
hold on

% exp2
sig1 = sig_env;
sig2 = r_sdgn2;

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
f = 1:size(phy_coh,1);f=f-1;f=f/size(phy_coh,1)*fs;
p2 = plot(f,mean(phy_coh,2),'color','k','linewidth',2,'linestyle','--');
hold on

sig1 = sig_env;
sig2 = r_sdgn3;

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
f = 1:size(phy_coh,1);f=f-1;f=f/size(phy_coh,1)*fs;
plot(f,mean(phy_coh,2),'color',[255 215 0]./255,'linewidth',2)
hold on
plot([4 4],[0 1],'color',[176 31 36]./255,'linestyle','--','linewidth',2);
xlim([1 14])
ylim([0 1])
xlabel('Frequency(Hz)')
ylabel('Coherence Spectrum')
legend([p1 p2],'Exp1','Exp2');
set(gca,'fontname','arial')
set(gca,'fontsize',13);
set(gcf,'unit','centimeters','position',[1.50 1.5 10.5 8])
title('0.125-s echoic speech')


%% depression and gain control(0.25s)
% 
clear
clc

addpath 'E:\echo_gjx\others\code\tools\'
addpath '.\tools'

tao_sd =  [0];
V_th = [0.15];
tao_gn =  [0.3];
fs = 100;

load(".\data_availability\Stimuli\anechoic_env.mat");
sig_env = [];


for i = 1:size(sig_speech_sum,2)
    sig_env = [sig_env;sig_speech_sum{i}];% original anechoic speech
end
% exp1
load(['.\data_availability\adaptation_model\model1_mes\exp1\rsdgn_250.mat'])
r_sdgn1 = r_sdgn;

% exp2
load(['.\data_availability\adaptation_model\model1_mes\exp2\rsdgn_250.mat'])
r_sdgn2 = r_sdgn;

load(['.\data_availability\adaptation_model\model1_mes\rsdgn_anechoic.mat'])
r_sdgn3 = r_sdgn;

figure
% exp1
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
f = 1:size(phy_coh,1);f=f-1;f=f/size(phy_coh,1)*fs;
p1 = plot(f,mean(phy_coh,2),'color','k','linewidth',2);
hold on

% exp2
sig1 = sig_env;
sig2 = r_sdgn2;

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
Color_temp = [[46 75 160]./255;[250 164 25]./255];
f = 1:size(phy_coh,1);f=f-1;f=f/size(phy_coh,1)*fs;
p2 = plot(f,mean(phy_coh,2),'color','k','linewidth',2,'linestyle','--');
hold on

sig1 = sig_env;
sig2 = r_sdgn3;

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
f = 1:size(phy_coh,1);f=f-1;f=f/size(phy_coh,1)*fs;
plot(f,mean(phy_coh,2),'color',[255 215 0]./255,'linewidth',2)
hold on
plot([2 2],[0 1],'color',[176 31 36]./255,'linestyle','--','linewidth',2);
plot([6 6],[0 1],'color',[176 31 36]./255,'linestyle','--','linewidth',2);
xlim([1 14])
ylim([0 1])
xlabel('Frequency(Hz)')
ylabel('Coherence Spectrum')
legend([p1 p2],'Exp1','Exp2');
% set(gca,'xticklabel',{'attend','unattend','vocode'})
set(gca,'fontname','arial')
set(gca,'fontsize',13);
set(gcf,'unit','centimeters','position',[1.50 1.5 10.5 8])
title('0.25-s echoic speech')