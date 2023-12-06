%% Fig 1 modulation spectrum
%% exp1
clear
clc
load('.\data_availability\Stimuli\exp1\clean_chunk_as.mat')
for ind = 1:size(v_sum,3)
    v = v_sum(:,:,ind);
    [ms1(:,ind),f]=Modulation_Spectrum(v,200,'log');
end
ms_rms1=sqrt(mean(ms1.^2,2));
nor_label1 = max(ms_rms1([f<32 & f>=.5]));
ms_rms_nor1=ms_rms1/nor_label1;


load('.\data_availability\Stimuli\exp1\echo_125_chunk_as.mat')
for ind = 1:size(v_sum,3)
    v = v_sum(:,:,ind);
    [ms2(:,ind),f]=Modulation_Spectrum(v,200,'log');
end
ms_rms2=sqrt(mean(ms2.^2,2));
nor_label2 = max(ms_rms2([f<32 & f>=.5]));
ms_rms_nor2=ms_rms2/nor_label1;

load('.\data_availability\Stimuli\exp1\echo_250_chunk_as.mat')
for ind = 1:size(v_sum,3)
    v = v_sum(:,:,ind);
    [ms3(:,ind),f]=Modulation_Spectrum(v,200,'log');
end
ms_rms3=sqrt(mean(ms3.^2,2));
nor_label3 = max(ms_rms3([f<32 & f>=.5]));
ms_rms_nor3=ms_rms3/nor_label1;

%%
figure
plot(f(1:size(ms1,1)),10*log10(ms_rms_nor1),'linewidth',2)
hold on
plot(f(1:size(ms2,1)),10*log10(ms_rms2/nor_label1),'linewidth',2)
plot(f(1:size(ms3,1)),10*log10(ms_rms3/nor_label1),'linewidth',2)

xlim([.5 32]);
set(gca,'xscale','log')
set(gca,'xtick',[.5 1 2 4 8 16 32])
xlabel('frequency (Hz)')
ylabel('normalized amplitude')
set(gca,'fontname','calibri');
set(gcf,'unit','centimeters','position',[4 4 15 10])
set(gca,'FontSize',14,'Linewidth',0.8,'GridAlpha',.15)
title('(Exp1) modulation spectrum','FontSize',15,'FontWeight','bold')
ylabel('dB')
%% exp 2 (figure)
clear
clc
load('.\data_availability\Stimuli\exp2\clean_chunk_as.mat')
for ind = 1:size(v_sum,3)
    v = v_sum(:,:,ind);
    [ms1(:,ind),f]=Modulation_Spectrum(v,200,'log');
end
ms_rms1=sqrt(mean(ms1.^2,2));
nor_label1 = max(ms_rms1([f<32 & f>=.5]));
ms_rms_nor1=ms_rms1/nor_label1;


load('.\data_availability\Stimuli\exp2\echo_125_chunk_as.mat')
for ind = 1:size(v_sum,3)
    v = v_sum(:,:,ind);
    [ms2(:,ind),f]=Modulation_Spectrum(v,200,'log');
end
ms_rms2=sqrt(mean(ms2.^2,2));
nor_label2 = max(ms_rms2([f<32 & f>=.5]));
ms_rms_nor2=ms_rms2/nor_label1;

load('.\data_availability\Stimuli\exp2\echo_250_chunk_as.mat')
for ind = 1:size(v_sum,3)
    v = v_sum(:,:,ind);
    [ms3(:,ind),f]=Modulation_Spectrum(v,200,'log');
end
ms_rms3=sqrt(mean(ms3.^2,2));
nor_label3 = max(ms_rms3([f<32 & f>=.5]));
ms_rms_nor3=ms_rms3/nor_label1;

%%
figure
plot(f(1:size(ms1,1)),10*log10(ms_rms_nor1),'linewidth',2)
hold on
plot(f(1:size(ms2,1)),10*log10(ms_rms2/nor_label1),'linewidth',2)
plot(f(1:size(ms3,1)),10*log10(ms_rms3/nor_label1),'linewidth',2)

xlim([.5 32]);
set(gca,'xscale','log')
set(gca,'xtick',[.5 1 2 4 8 16 32])
title('(Exp2) modulation spectrum','FontSize',15,'FontWeight','bold')
xlabel('frequency (Hz)')
ylabel('normalized amplitude')
set(gca,'fontname','calibri');
set(gcf,'unit','centimeters','position',[4 4 15 10])
set(gca,'FontSize',14,'Linewidth',0.8,'GridAlpha',.15)
ylabel('dB')