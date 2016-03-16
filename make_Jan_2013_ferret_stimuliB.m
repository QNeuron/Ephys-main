%
% this script generates stimuli for ferret behavioral pitch discrimination
% experiments
%
% Jan 7 2013

close all

sr =  97656;
stim_dur_ms = 200;

noise_l_co = 100;
noise_h_co = 10000;

f0_1 = 500;
f0_2 = 1000;

rms_scale_factor = 100;


% Stage 2 - tones with high harmonics
%
% 3000-7200 Hz
% tone 1: 300 Hz F0, harmonics 10-24; tone 2: 600 Hz, harmonics 5-12

l_co = 4000;
u_co = 10000;

tone2_hh = [];
for h = l_co/f0_2 : u_co/f0_2
    harm = tone(f0_2*h,stim_dur_ms,0,sr);
    tone2_hh = zadd(tone2_hh,harm);
end
tone2_hh = tone2_hh/rms_scale_factor;
 
tone1_hh = [];
for h = l_co/f0_1 : u_co/f0_1
    harm = tone(f0_1*h,stim_dur_ms,0,sr);
    tone1_hh = zadd(tone1_hh,harm);
end
tone1_hh = tone1_hh/rms(tone1_hh)*rms(tone2_hh); %level difference was ~3dB ==> primaries of tone2 are 3dB higher than those of tone1
%this should be offset by the greater number of harmonics in tone1, leading
%to approximately equal amplitude DPs in the two stimuli
%we assume the DPs to be about 15 dB below the level of the primaries

%determine factor tone_adjustment by which tones must be adjusted to yield 65 dB SPL
%through Kerry's system:

tone_adjustment = 1;

tone1_hh = tone1_hh*tone_adjustment;
tone2_hh = tone2_hh*tone_adjustment;


dp_noise = pnoise(stim_dur_ms,noise_l_co,noise_h_co,-30,0,sr);

% ERBs for ferret: 180 Hz at 300 Hz, 225 Hz at 500 Hz, 250 Hz at 600 Hz, 300 Hz at 1000 Hz
%cf = 300;
%erb_width = 180;
%cf=600; erb_width=250;
%calibrating with level at a higher F0 is more conservative as 
%filterer narrower ==> less noise power
% cf = 600;
% erb_width = 300;
cf = 600;
erb_width = 300;
erb_low = (-erb_width + sqrt(erb_width^2 +4*cf^2))/2;
erb_high = erb_low+erb_width;

dp_noise_fft = fft(dp_noise);

%signal is odd length
nfreqs = (length(dp_noise)-1)/2;
max_freq = sr*(length(dp_noise)-1)/2/length(dp_noise); %max freq is just under nyquist
freqs = [0:max_freq/nfreqs:max_freq];
neg_freqs = fliplr(freqs(2:end));

[temp, low_bin] = min(abs(freqs-erb_low));
[temp, high_bin] = min(abs(freqs-erb_high));

[temp, low_bin_neg] = min(abs(neg_freqs-erb_low));
[temp, high_bin_neg] = min(abs(neg_freqs-erb_high));
low_bin_neg = low_bin_neg+length(freqs);
high_bin_neg = high_bin_neg+length(freqs);

dp_noise_fft_filt = dp_noise_fft;
dp_noise_fft_filt([1:low_bin-1 high_bin+1:high_bin_neg-1 low_bin_neg+1:length(dp_noise)]) = 0;
dp_noise_filt = ifft(dp_noise_fft_filt);

%determine factor by which noise band must be multiplied to be 5 dB below
%primaries, i.e. 10 dB above hypothetical DPs which are supposed to be 15
%dB below primaries (based on Pressnitzer and Patterson data)
%
% P&P data show that DPs grow with number of primaries
% DP amplitude increases ~3dB for each doubling of the number of primaries
% in their experiment primaries were 60 db SPL and with 17 primaries (the
% most they tried), the DP was ~45 db SPL.
% Our high-harmonic tones have 15 or 8 harmonics (low/high f0), so
% accounting for a DP that is 15 dB below the primaries seems conservative.

sample_harmonic = tone(f0_2*10,stim_dur_ms,0,sr);

rms_primary = rms(sample_harmonic)/rms_scale_factor;

noise_adjustment = 1/rms(dp_noise_filt) * rms_primary * 10^(-5/20);

dp_noise_adjusted = dp_noise*noise_adjustment*tone_adjustment;



tone1_high_harm_stim = zadd(tone1_hh, dp_noise_adjusted);
tone2_high_harm_stim = zadd(tone2_hh, dp_noise_adjusted);
% tone1_high_harm_stim = tone1_hh;
% tone2_high_harm_stim = tone2_hh;


% Stage 1 - full bandwidth tones
%
% 600-7200 Hz
% tone 1: 300 Hz F0, harmonics 2-24; tone 2: 600 Hz, harmonics 1-12

l_co = 1000;
u_co = 10000;

tone2_ah = [];
for h = l_co/f0_2 : u_co/f0_2
    harm = tone(f0_2*h,stim_dur_ms,0,sr);
    tone2_ah = zadd(tone2_ah,harm);
end
tone2_ah = tone2_ah/rms(tone2_ah)*rms(tone2_hh);

tone1_ah = [];
for h = l_co/f0_1 : u_co/f0_1
    harm = tone(f0_1*h,stim_dur_ms,0,sr);
    tone1_ah = zadd(tone1_ah,harm);
end
tone1_ah = tone1_ah/rms(tone1_ah)*rms(tone2_ah); %level difference was ~3dB ==> primaries of tone2 are 3dB higher than those of tone1


tone1_all_harm_stim = zadd(tone1_ah, dp_noise_adjusted);
tone2_all_harm_stim = zadd(tone2_ah, dp_noise_adjusted);



% Stage 3 - tones with lower harmonics
%
% 600-3000 Hz
% tone 1: 300 Hz F0, harmonics 2-10; tone 2: 600 Hz, harmonics 1-5

l_co = 1000;
u_co = 4000;

tone2_lh = [];
for h = l_co/f0_2 : u_co/f0_2
    harm = tone(f0_2*h,stim_dur_ms,0,sr);
    tone2_lh = zadd(tone2_lh,harm);
end
tone2_lh = tone2_lh/rms(tone2_lh)*rms(tone2_hh);

tone1_lh = [];
for h = l_co/f0_1 : u_co/f0_1
    harm = tone(f0_1*h,stim_dur_ms,0,sr);
    tone1_lh = zadd(tone1_lh,harm);
end
tone1_lh = tone1_lh/rms(tone1_lh)*rms(tone2_ah); %level difference was ~3dB ==> primaries of tone2 are 3dB higher than those of tone1


tone1_low_harm_stim = zadd(tone1_lh, dp_noise_adjusted);
tone2_low_harm_stim = zadd(tone2_lh, dp_noise_adjusted);


% Stage 4 - ALT phase tones with higher harmonics
%
% 600-3000 Hz
% tone 1: 300 Hz F0, harmonics 2-10; tone 2: 600 Hz, harmonics 1-5
l_co = 4000;
u_co = 10000;


tone2_hh_alt = [];
for h = l_co/f0_2 : u_co/f0_2
    if mod(h,2)==1,
        harm = tone(f0_2*h,stim_dur_ms,0,sr);
    else
        harm = tone(f0_2*h,stim_dur_ms,pi,sr);
    end
    tone2_hh_alt = zadd(tone2_hh_alt,harm);
end
tone2_hh_alt = tone2_hh_alt/rms(tone2_hh_alt)*rms(tone2_hh);

tone1_hh_alt = [];
for h = l_co/f0_1 : u_co/f0_1
    if mod(h,2)==1,
        harm = tone(f0_1*h,stim_dur_ms,0,sr);
    else
        harm = tone(f0_1*h,stim_dur_ms,pi,sr);
    end
    tone1_hh_alt = zadd(tone1_hh_alt,harm);
end
tone1_hh_alt = tone1_hh_alt/rms(tone1_hh_alt)*rms(tone2_hh);


tone1_high_harm_alt_stim = zadd(tone1_hh_alt, dp_noise_adjusted);
tone2_high_harm_alt_stim = zadd(tone2_hh_alt, dp_noise_adjusted);
% tone1_high_harm_alt_stim = tone1_hh_alt;
% tone2_high_harm_alt_stim = tone2_hh_alt;


% Stage 5 - RAND phase tones with higher harmonics
%
% 600-3000 Hz
% tone 1: 300 Hz F0, harmonics 2-10; tone 2: 600 Hz, harmonics 1-5
l_co = 4000;
u_co = 10000;


tone2_hh_rand = [];
for h = l_co/f0_2 : u_co/f0_2
    thisphase=rand*2*pi;
%     thisphase=rand;
    harm = tone(f0_2*h,stim_dur_ms,thisphase,sr);
    tone2_hh_rand = zadd(tone2_hh_rand,harm);
end
tone2_hh_rand = tone2_hh_rand/rms(tone2_hh_rand)*rms(tone2_hh);


tone1_hh_rand = [];
for h = l_co/f0_1 : u_co/f0_1
    thisphase=rand*2*pi;
%     thisphase=rand;
    harm = tone(f0_1*h,stim_dur_ms,thisphase,sr);
    tone1_hh_rand = zadd(tone1_hh_rand,harm);
end
tone1_hh_rand = tone1_hh_rand/rms(tone1_hh_rand)*rms(tone2_hh);

tone1_high_harm_rand_stim = zadd(tone1_hh_rand, dp_noise_adjusted);
tone2_high_harm_rand_stim = zadd(tone2_hh_rand, dp_noise_adjusted);
% tone1_high_harm_rand_stim = tone1_hh_rand;
% tone2_high_harm_rand_stim = tone2_hh_rand;


% Stage 6 - RAND phase tones with all (or variable number of lower) harmonics
l_co = 1000:1000:4000;
u_co = 10000;
tone1_all_harm_rand_stim=[];
tone2_all_harm_rand_stim=[];
for ii=1:length(l_co)
    tone2_ah_rand = [];
    for h = l_co(ii)/f0_2 : u_co/f0_2
        thisphase=rand*2*pi;
    %     thisphase=rand;
        harm = tone(f0_2*h,stim_dur_ms,thisphase,sr);
        tone2_ah_rand = zadd(tone2_ah_rand,harm);
    end
    tone2_ah_rand = tone2_ah_rand/rms(tone2_ah_rand)*rms(tone2_hh);

    tone1_ah_rand = [];
    for h = l_co(ii)/f0_1 : u_co/f0_1
        thisphase=rand*2*pi;
    %     thisphase=rand;
        harm = tone(f0_1*h,stim_dur_ms,thisphase,sr);
        tone1_ah_rand = zadd(tone1_ah_rand,harm);
    end
    tone1_ah_rand = tone1_ah_rand/rms(tone1_ah_rand)*rms(tone2_hh);

    tone1_all_harm_rand_stim(ii,:) = zadd(tone1_ah_rand, dp_noise_adjusted);
    tone2_all_harm_rand_stim(ii,:) = zadd(tone2_ah_rand, dp_noise_adjusted);
end
l_co_rand=l_co;



% envelop sounds
gatelength=5; %5ms onset and offset ramp
tone1_ah=linear_envelope(tone1_ah,gatelength,sr);
tone2_ah=linear_envelope(tone2_ah,gatelength,sr);
tone1_all_harm_stim=linear_envelope(tone1_all_harm_stim,gatelength,sr);
tone2_all_harm_stim=linear_envelope(tone2_all_harm_stim,gatelength,sr);
tone1_high_harm_stim=linear_envelope(tone1_high_harm_stim,gatelength,sr);
tone2_high_harm_stim=linear_envelope(tone2_high_harm_stim,gatelength,sr);
tone1_low_harm_stim=linear_envelope(tone1_low_harm_stim,gatelength,sr);
tone2_low_harm_stim=linear_envelope(tone2_low_harm_stim,gatelength,sr);
tone1_high_harm_alt_stim=linear_envelope(tone1_high_harm_alt_stim,gatelength,sr);
tone2_high_harm_alt_stim=linear_envelope(tone2_high_harm_alt_stim,gatelength,sr);
tone1_high_harm_rand_stim=linear_envelope(tone1_high_harm_rand_stim,gatelength,sr);
tone2_high_harm_rand_stim=linear_envelope(tone2_high_harm_rand_stim,gatelength,sr);
for ii=1:size(tone1_all_harm_rand_stim,1)
    tone1_all_harm_rand_stim(ii,:)=linear_envelope(tone1_all_harm_rand_stim(ii,:),gatelength,sr);
    tone2_all_harm_rand_stim(ii,:)=linear_envelope(tone2_all_harm_rand_stim(ii,:),gatelength,sr);
end

%%
% plotting stimuli
figpos=[1283 -233 1662 420];
xlimits=[0 20];

figure(1);clf %spectra
% set(gcf,'position',figpos)
subplot(2,5,1);
pwelch(tone1_all_harm_stim,[],[],[],sr);
set(gca,'YLim',[-100 -50]);
title('All harmonics')
xlim(xlimits)
subplot(2,5,2);
pwelch(tone1_high_harm_stim,[],[],[],sr);
set(gca,'YLim',[-100 -50]);
title('High harmonics')
xlim(xlimits)
subplot(2,5,3);
pwelch(tone1_low_harm_stim,[],[],[],sr);
set(gca,'YLim',[-100 -50]);
title('Low harmonics')
xlim(xlimits)
% subplot(2,5,4);
% pwelch(tone1_high_harm_alt_stim,[],[],[],sr);
% set(gca,'YLim',[-100 -50]);
% title('ALT phase')
% xlim(xlimits)
subplot(2,5,4);
pwelch(tone1_high_harm_rand_stim,[],[],[],sr);
set(gca,'YLim',[-100 -50]);
title('Rand phase, Hi harm')
xlim(xlimits)
subplot(2,5,5);
pwelch(tone1_all_harm_rand_stim(1,:),[],[],[],sr);
set(gca,'YLim',[-100 -50]);
title('Rand phase, All harm')
xlim(xlimits)
subplot(2,5,6);
pwelch(tone2_all_harm_stim,[],[],[],sr);
set(gca,'YLim',[-100 -50]);
title('')
xlim(xlimits)
subplot(2,5,7);
pwelch(tone2_high_harm_stim,[],[],[],sr);
set(gca,'YLim',[-100 -50]);
title('')
xlim(xlimits)
subplot(2,5,8);
pwelch(tone2_low_harm_stim,[],[],[],sr);
set(gca,'YLim',[-100 -50]);
title('')
xlim(xlimits)
% subplot(2,5,9);
% pwelch(tone2_high_harm_alt_stim,[],[],[],sr);
% set(gca,'YLim',[-100 -50]);
% title('')
% xlim(xlimits)
subplot(2,5,9);
pwelch(tone2_high_harm_rand_stim,[],[],[],sr);
set(gca,'YLim',[-100 -50]);
title('')
xlim(xlimits)
subplot(2,5,10);
pwelch(tone2_all_harm_rand_stim(1,:),[],[],[],sr);
set(gca,'YLim',[-100 -50]);
title('')
xlim(xlimits)


% figure %spectra on a log scale
% set(gcf,'position',figpos)
% subplot(2,5,1);
% pwelch(tone1_all_harm_stim,[],[],[],sr);
% set(gca,'YLim',[-100 -50],'XLim',[.1 8],'XScale','log');
% subplot(2,5,2);
% pwelch(tone1_high_harm_stim,[],[],[],sr);
% set(gca,'YLim',[-100 -50],'XLim',[.1 8],'XScale','log');
% subplot(2,5,3);
% pwelch(tone1_low_harm_stim,[],[],[],sr);
% set(gca,'YLim',[-100 -50],'XLim',[.1 8],'XScale','log');
% subplot(2,5,4);
% pwelch(tone1_high_harm_alt_stim,[],[],[],sr);
% set(gca,'YLim',[-100 -50],'XLim',[.1 8],'XScale','log');
% subplot(2,5,5);
% pwelch(tone1_high_harm_rand_stim,[],[],[],sr);
% set(gca,'YLim',[-100 -50],'XLim',[.1 8],'XScale','log');
% subplot(2,5,6);
% pwelch(tone2_all_harm_stim,[],[],[],sr);
% set(gca,'YLim',[-100 -50],'XLim',[.1 8],'XScale','log');
% subplot(2,5,7);
% pwelch(tone2_high_harm_stim,[],[],[],sr);
% set(gca,'YLim',[-100 -50],'XLim',[.1 8],'XScale','log');
% subplot(2,5,8);
% pwelch(tone2_low_harm_stim,[],[],[],sr);
% set(gca,'YLim',[-100 -50],'XLim',[.1 8],'XScale','log');
% subplot(2,5,9);
% pwelch(tone2_high_harm_alt_stim,[],[],[],sr);
% set(gca,'YLim',[-100 -50],'XLim',[.1 8],'XScale','log');
% subplot(2,5,10);
% pwelch(tone2_high_harm_rand_stim,[],[],[],sr);
% set(gca,'YLim',[-100 -50],'XLim',[.1 8],'XScale','log');



figure(2);clf %waveform
% set(gcf,'position',figpos)
% plotdur=length(tone1_all_harm_stim);%duration to plot
plotdur=.05*sr; %plot just first 'x' ms
subplot(2,5,1);
plot([1:plotdur]./sr,tone1_all_harm_stim(1:plotdur));
set(gca,'YLim',[-.2 .2],'XLim',[0 0.05]);
title('All harmonics')
subplot(2,5,2);
plot([1:plotdur]./sr,tone1_high_harm_stim(1:plotdur));
set(gca,'YLim',[-.2 .2],'XLim',[0 0.05]);
title('High harmonics')
subplot(2,5,3);
plot([1:plotdur]./sr,tone1_low_harm_stim(1:plotdur));
set(gca,'YLim',[-.2 .2],'XLim',[0 0.05]);
title('Low harmonics')
% subplot(2,5,4);
% plot([1:plotdur]./sr,tone1_high_harm_alt_stim(1:plotdur));
% set(gca,'YLim',[-.2 .2],'XLim',[0 0.05]);
% title('ALT phase')
subplot(2,5,4);
plot([1:plotdur]./sr,tone1_high_harm_rand_stim(1:plotdur));
set(gca,'YLim',[-.2 .2],'XLim',[0 0.05]);
title('Rand phase, Hi harm')
subplot(2,5,5);
plot([1:plotdur]./sr,tone1_all_harm_rand_stim(1,1:plotdur));
set(gca,'YLim',[-.2 .2],'XLim',[0 0.05]);
title('Rand phase, All harm')
subplot(2,5,6);
plot([1:plotdur]./sr,tone2_all_harm_stim(1:plotdur));
set(gca,'YLim',[-.2 .2],'XLim',[0 0.05]);
subplot(2,5,7);
plot([1:plotdur]./sr,tone2_high_harm_stim(1:plotdur));
set(gca,'YLim',[-.2 .2],'XLim',[0 0.05]);
subplot(2,5,8);
plot([1:plotdur]./sr,tone2_low_harm_stim(1:plotdur));
set(gca,'YLim',[-.2 .2],'XLim',[0 0.05]);
% subplot(2,5,9);
% plot([1:plotdur]./sr,tone2_high_harm_alt_stim(1:plotdur));
% set(gca,'YLim',[-.2 .2],'XLim',[0 0.05]);
subplot(2,5,9);
plot([1:plotdur]./sr,tone2_high_harm_rand_stim(1:plotdur));
set(gca,'YLim',[-.2 .2],'XLim',[0 0.05]);
subplot(2,5,10);
plot([1:plotdur]./sr,tone2_all_harm_rand_stim(1,1:plotdur));
set(gca,'YLim',[-.2 .2],'XLim',[0 0.05]);


xlabel('Time (s)')


%% listening to phase manipulations

% fullstimHIGH=[tone1_all_harm_stim zeros(1,ceil(sr*.5)) tone1_high_harm_stim zeros(1,ceil(sr*.5))];
% fullstimLOW=[tone1_all_harm_stim zeros(1,ceil(sr*.5)) tone1_low_harm_stim zeros(1,ceil(sr*.5))];
% fullstimALT=[tone1_high_harm_stim zeros(1,ceil(sr*.5)) tone1_high_harm_alt_stim zeros(1,ceil(sr*.5))];
% fullstimRAND=[tone1_high_harm_stim zeros(1,ceil(sr*.5)) tone1_high_harm_rand_stim zeros(1,ceil(sr*.5))];
% 
% soundsc(repmat(fullstimHIGH,1,5),sr)
% 
% soundsc(repmat(fullstimLOW,1,5),sr)
% 
% soundsc(repmat(fullstimALT,1,5),sr)
% 
% soundsc(repmat(fullstimRAND,1,5),sr)

%% Save sounds
% close all
% save JoshSoundsB.mat tone1_ah tone1_all_harm_stim tone1_high_harm_stim tone1_low_harm_stim tone1_high_harm_alt_stim tone1_high_harm_rand_stim...
%     tone2_ah tone2_all_harm_stim tone2_high_harm_stim tone2_low_harm_stim tone2_high_harm_alt_stim tone2_high_harm_rand_stim...
%     sr stim_dur_ms f0_1 f0_2...
%     noise_l_co noise_h_co;

