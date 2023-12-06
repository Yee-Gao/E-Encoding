%% synaptic depression (David 2009)
clear
clc

addpath '.\tools'
fs = 100;


load('.\data_availability\Stimuli\anechoic_env.mat');
sig_env = [];

for i = 1:size(sig_speech_sum,2)
    sig_env = [sig_env;sig_speech_sum{i}];% original anechoic speech
end
%% 0.125-s echoic speech

tao = 100;
miu = 0.2;

% exp 1
load('.\data_availability\adaptation_model\model2_dav\exp1\125.mat')
sig1 = sig_env;
sig2 = sig_env_dep;


for ind = 1:1000
    try x = sig1(fs*2*(ind-1)+1:fs*2*(ind));
    catch,disp([num2str(ind)]);break;end
    x_sum(:,ind) = x;

    try y = sig2(fs*2*(ind-1)+1:fs*2*(ind));
    catch,disp([num2str(ind)]);break;end
    y_sum(:,ind) = y;

end
fre_x = fft(x_sum);
fre_y = fft(y_sum);
phy = angle(fre_x)-angle(fre_y);
phy_coh_dep= mean(cos(phy),2).^2+mean(sin(phy),2).^2;

figure

f = 1:size(phy_coh_dep,1);f=f-1;f=f/size(phy_coh_dep,1)*fs;
plot(f,mean(phy_coh_dep,2),'color','k','linewidth',2)
xlim([1 14])
ylim([0 1])
set(gca,'fontsize',13)
set(gca,'xtick',0:2:120)
xlabel('Frequency(Hz)')
ylabel('Coherence Spectrum')
set(gca,'fontname','arial')
set(gca,'fontsize',13);
set(gcf,'unit','centimeters','position',[1.50 1.5 10.5 8])
title('depressed anechoic speech')

hold on
% exp 2
load('.\data_availability\adaptation_model\model2_dav\exp2\125.mat');
sig1 = sig_env;
sig2 = sig_env_dep;


for ind = 1:1000
    try x = sig1(fs*2*(ind-1)+1:fs*2*(ind));
    catch,disp([num2str(ind)]);break;end
    x_sum(:,ind) = x;

    try y = sig2(fs*2*(ind-1)+1:fs*2*(ind));
    catch,disp([num2str(ind)]);break;end
    y_sum(:,ind) = y;

end
fre_x = fft(x_sum);
fre_y = fft(y_sum);
phy = angle(fre_x)-angle(fre_y);
phy_coh_dep= mean(cos(phy),2).^2+mean(sin(phy),2).^2;
p1 = plot(f,mean(phy_coh_dep,2),'color','k','linewidth',2,'LineStyle','-.');

load('.\data_availability\adaptation_model\model2_dav\anechoic.mat');
sig1 = sig_env;
sig2 = sig_env_dep;


for ind = 1:1000
    try x = sig1(fs*2*(ind-1)+1:fs*2*(ind));
    catch,disp([num2str(ind)]);break;end
    x_sum(:,ind) = x;

    try y = sig2(fs*2*(ind-1)+1:fs*2*(ind));
    catch,disp([num2str(ind)]);break;end
    y_sum(:,ind) = y;

end
fre_x = fft(x_sum);
fre_y = fft(y_sum);
phy = angle(fre_x)-angle(fre_y);
phy_coh_dep= mean(cos(phy),2).^2+mean(sin(phy),2).^2;
p2 = plot(f,mean(phy_coh_dep,2),'color',[255 215 0]./255,'linewidth',2);

legend([p1 p2],'Exp1','Exp2');
plot([4 4],[0 1],'color',[176 31 36]./255,'linestyle','--','linewidth',2);

xlim([1 14])
ylim([0 1])
set(gca,'fontsize',13)
set(gca,'xtick',0:2:120)
xlabel('Frequency(Hz)')
ylabel('Phase coherence')
set(gca,'fontname','arial')
set(gca,'fontsize',13);
set(gcf,'unit','centimeters','position',[1.5 1.5 13 10])
title('0.125-s echoic speech')
hold on



%% 0.25-s echoic speech

tao = 100;
miu = 0.2;

% exp 1
load('.\data_availability\adaptation_model\model2_dav\exp1\250.mat');
sig1 = sig_env;
sig2 = sig_env_dep;


for ind = 1:1000
    try x = sig1(fs*2*(ind-1)+1:fs*2*(ind));
    catch,disp([num2str(ind)]);break;end
    x_sum(:,ind) = x;

    try y = sig2(fs*2*(ind-1)+1:fs*2*(ind));
    catch,disp([num2str(ind)]);break;end
    y_sum(:,ind) = y;

end
fre_x = fft(x_sum);
fre_y = fft(y_sum);
phy = angle(fre_x)-angle(fre_y);
phy_coh_dep= mean(cos(phy),2).^2+mean(sin(phy),2).^2;

figure

f = 1:size(phy_coh_dep,1);f=f-1;f=f/size(phy_coh_dep,1)*fs;
p1 = plot(f,mean(phy_coh_dep,2),'color','k','linewidth',2);
xlim([1 14])
ylim([0 1])
set(gca,'fontsize',13)
set(gca,'xtick',0:2:120)
xlabel('Frequency(Hz)')
ylabel('Coherence Spectrum')
set(gca,'fontname','arial')
set(gca,'fontsize',13);
set(gcf,'unit','centimeters','position',[1.50 1.5 10.5 8])
title('depressed anechoic speech')

hold on
% exp 2
load('.\data_availability\adaptation_model\model2_dav\exp2\250.mat');
sig1 = sig_env;
sig2 = sig_env_dep;


for ind = 1:1000
    try x = sig1(fs*2*(ind-1)+1:fs*2*(ind));
    catch,disp([num2str(ind)]);break;end
    x_sum(:,ind) = x;

    try y = sig2(fs*2*(ind-1)+1:fs*2*(ind));
    catch,disp([num2str(ind)]);break;end
    y_sum(:,ind) = y;

end
fre_x = fft(x_sum);
fre_y = fft(y_sum);
phy = angle(fre_x)-angle(fre_y);
phy_coh_dep= mean(cos(phy),2).^2+mean(sin(phy),2).^2;
p2=plot(f,mean(phy_coh_dep,2),'color','k','linewidth',2,'LineStyle','-.');

load('.\data_availability\adaptation_model\model2_dav\anechoic.mat');
% load(['..\PNAS_model\Syndep_David\miu' num2str(miu) 'tao' num2str(tao) 'clean.mat'])
sig1 = sig_env;
sig2 = sig_env_dep;


for ind = 1:1000
    try x = sig1(fs*2*(ind-1)+1:fs*2*(ind));
    catch,disp([num2str(ind)]);break;end
    x_sum(:,ind) = x;

    try y = sig2(fs*2*(ind-1)+1:fs*2*(ind));
    catch,disp([num2str(ind)]);break;end
    y_sum(:,ind) = y;

end
fre_x = fft(x_sum);
fre_y = fft(y_sum);
phy = angle(fre_x)-angle(fre_y);
phy_coh_dep= mean(cos(phy),2).^2+mean(sin(phy),2).^2;
plot(f,mean(phy_coh_dep,2),'color',[255 215 0]./255,'linewidth',2)

xlim([1 14])
ylim([0 1])
set(gca,'fontsize',13)
set(gca,'xtick',0:2:120)
xlabel('Frequency(Hz)')
ylabel('Phase coherence')
legend([p1 p2],'Exp1','Exp2');
set(gca,'fontname','arial')
set(gca,'fontsize',13);
set(gcf,'unit','centimeters','position',[1.5 1.5 13 10])
title('0.25-s echoic speech')
hold on
plot([2 2],[0 1],'color',[176 31 36]./255,'linestyle','--','linewidth',2);
plot([6 6],[0 1],'color',[176 31 36]./255,'linestyle','--','linewidth',2);

%% Optimal filter [0.125-s echoic speech]

clear
clc

addpath '.\tools'
fs = 100;

% exp 1

load('.\data_availability\Stimuli\125_exp1_env.mat');
load('.\data_availability\Stimuli\anechoic_env.mat');

sig_env_echo = [];
sig_env = [];


for i = 1:size(sig_speech_echo_sum,2)
    sig_env_echo = [sig_env_echo;sig_speech_echo_sum{i}]; % echoic speech
    sig_env = [sig_env;sig_speech_sum{i}];% original anechoic speech
end

addpath '.\tools\trf'
shift = 0.5;

x = sig_env-mean(sig_env);
y = sig_env_echo-mean(sig_env_echo);
[h,cr] = normRCtik_Z(y',circshift(x,round(shift*fs))',128*2,0.0001);
t = 1:size(h,1);t=t-1;t=t/fs;
t = -round(shift*fs):1:size(h,1)-round(shift*fs)-1;t = t/fs;

x_pred = filter(h,1,y);
x_circ = circshift(x,round(shift*fs));

for ind = 1:1000
    try  x = x_pred(fs*2*(ind-1)+1:fs*2*ind);
    catch,disp([num2str(ind)]);break;end
    sig_pred(:,ind) = x;
    try  y = x_circ(fs*2*(ind-1)+1:fs*2*ind);
    catch,break;end
    sig_circ(:,ind) = y;
end
fre_pred = fft(sig_pred);
fre_circ = fft(sig_circ);
phy = angle(fre_pred)-angle(fre_circ);
phy_coh= mean(cos(phy),2).^2+mean(sin(phy),2).^2;

figure
f = 1:length(phy_coh);f=f-1;f=f/length(phy_coh)*fs;
p1 = plot(f,mean(phy_coh,2),'color','k','linewidth',2);


% exp 2

load('.\data_availability\Stimuli\125_exp2_env.mat');
load('.\data_availability\Stimuli\anechoic_env.mat');

sig_env_echo = [];
sig_env = [];


for i = 1:size(sig_speech_echo_sum,2)
    sig_env_echo = [sig_env_echo;sig_speech_echo_sum{i}]; % echoic speech
    sig_env = [sig_env;sig_speech_sum{i}];% original anechoic speech
end

shift = 0.5;

x = sig_env-mean(sig_env);
y = sig_env_echo-mean(sig_env_echo);
[h,cr] = normRCtik_Z(y',circshift(x,round(shift*fs))',128*2,0.0001);
t = 1:size(h,1);t=t-1;t=t/fs;
t = -round(shift*fs):1:size(h,1)-round(shift*fs)-1;t = t/fs;

x_pred = filter(h,1,y);
x_circ = circshift(x,round(shift*fs));

for ind = 1:1000
    try  x = x_pred(fs*2*(ind-1)+1:fs*2*ind);
    catch,disp([num2str(ind)]);break;end
    sig_pred(:,ind) = x;
    try  y = x_circ(fs*2*(ind-1)+1:fs*2*ind);
    catch,break;end
    sig_circ(:,ind) = y;
end
fre_pred = fft(sig_pred);
fre_circ = fft(sig_circ);
phy = angle(fre_pred)-angle(fre_circ);
phy_coh= mean(cos(phy),2).^2+mean(sin(phy),2).^2;

hold on
f = 1:length(phy_coh);f=f-1;f=f/length(phy_coh)*fs;
p2 = plot(f,mean(phy_coh,2),'color','k','linewidth',2,'LineStyle','-.');
plot([4 4],[0 1],'color',[176 31 36]./255,'linestyle','--','linewidth',2);

xlim([1 14])
ylim([0 1])
hold on
xlabel('Frequency(Hz)')
ylabel('Phase coherence')
legend('Exp1','Exp2');
set(gca,'fontname','arial')
set(gca,'fontsize',13);
set(gcf,'unit','centimeters','position',[1.50 1.5 13 10])
title('0.125-s echoic speech')
legend([p1 p2],'Exp1','Exp2');

%% Optimal filter (0.25-s echoic speech)
clear
clc

% addpath 'E:\echo_gjx\others\code\tools\'
addpath '.\tools'
fs = 100;

% exp 1

load('.\data_availability\Stimuli\250_exp1_env.mat');
load('.\data_availability\Stimuli\anechoic_env.mat');

sig_env_echo = [];
sig_env = [];


for i = 1:size(sig_speech_echo_sum2,2)
    sig_env_echo = [sig_env_echo;sig_speech_echo_sum2{i}]; % echoic speech
    sig_env = [sig_env;sig_speech_sum{i}];% original anechoic speech
end

shift = 0.5;

x = sig_env-mean(sig_env);
y = sig_env_echo-mean(sig_env_echo);
[h,cr] = normRCtik_Z(y',circshift(x,round(shift*fs))',128*2,0.0001);
t = 1:size(h,1);t=t-1;t=t/fs;
t = -round(shift*fs):1:size(h,1)-round(shift*fs)-1;t = t/fs;

x_pred = filter(h,1,y);
x_circ = circshift(x,round(shift*fs));

for ind = 1:1000
    try  x = x_pred(fs*2*(ind-1)+1:fs*2*ind);
    catch,disp([num2str(ind)]);break;end
    sig_pred(:,ind) = x;
    try  y = x_circ(fs*2*(ind-1)+1:fs*2*ind);
    catch,break;end
    sig_circ(:,ind) = y;
end
fre_pred = fft(sig_pred);
fre_circ = fft(sig_circ);
phy = angle(fre_pred)-angle(fre_circ);
phy_coh= mean(cos(phy),2).^2+mean(sin(phy),2).^2;

figure
f = 1:length(phy_coh);f=f-1;f=f/length(phy_coh)*fs;
p1 = plot(f,mean(phy_coh,2),'color','k','linewidth',2);


% exp 2

load('.\data_availability\Stimuli\250_exp2_env.mat');
load('.\data_availability\Stimuli\anechoic_env.mat');

sig_env_echo = [];
sig_env = [];


for i = 1:size(sig_speech_echo_sum2,2)
    sig_env_echo = [sig_env_echo;sig_speech_echo_sum2{i}]; % echoic speech
    sig_env = [sig_env;sig_speech_sum{i}];% original anechoic speech
end

shift = 0.5;

x = sig_env-mean(sig_env);
y = sig_env_echo-mean(sig_env_echo);
[h,cr] = normRCtik_Z(y',circshift(x,round(shift*fs))',128*2,0.0001);
t = 1:size(h,1);t=t-1;t=t/fs;
t = -round(shift*fs):1:size(h,1)-round(shift*fs)-1;t = t/fs;

x_pred = filter(h,1,y);
x_circ = circshift(x,round(shift*fs));

for ind = 1:1000
    try  x = x_pred(fs*2*(ind-1)+1:fs*2*ind);
    catch,disp([num2str(ind)]);break;end
    sig_pred(:,ind) = x;
    try  y = x_circ(fs*2*(ind-1)+1:fs*2*ind);
    catch,break;end
    sig_circ(:,ind) = y;
end
fre_pred = fft(sig_pred);
fre_circ = fft(sig_circ);
phy = angle(fre_pred)-angle(fre_circ);
phy_coh= mean(cos(phy),2).^2+mean(sin(phy),2).^2;

hold on
f = 1:length(phy_coh);f=f-1;f=f/length(phy_coh)*fs;
p2 = plot(f,mean(phy_coh,2),'color','k','linewidth',2,'LineStyle','-.');
plot([2 2],[0 1],'color',[176 31 36]./255,'linestyle','--','linewidth',2);
plot([6 6],[0 1],'color',[176 31 36]./255,'linestyle','--','linewidth',2);

xlim([1 14])
ylim([0 1])
hold on
xlabel('Frequency(Hz)')
ylabel('Phase coherence')
legend('Exp1','Exp2');
set(gca,'fontname','arial')
set(gca,'fontsize',13);
set(gcf,'unit','centimeters','position',[1.50 1.5 13 10])
title('0.25-s echoic speech')
legend([p1 p2],'Exp1','Exp2');
