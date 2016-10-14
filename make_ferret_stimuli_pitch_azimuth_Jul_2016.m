%% Make pitch sounds with varying azimuth
% Stim generation
% Head transfer functions must be in CIPIC format. Use schnup_to_cipic.m to
% convert raw ferret htrf to CIPIC htrf.
%
% Quentin 07 2016
% Adapted from make_ferret_stimuli_physiol_Aug_2013.m
%
% Quentin 10 2016
% Works with Jan HTRF directly
% Resample
% 
% 
% Make sure the parameters of your stimuli fits with those from the non
% azimuth modulated ones !

%% set paths
oldPath = path;
newpath = genpath('D:\Work\Code\KingLab\VAS_test\');
path(oldPath,newpath)


%% Make stim


close all
clear

sr = tdt100k;
stim_dur_ms = 200;

lowcutoff=200;
midcutoff=4000;
highcutoff=10000;

noise_l_co = 200;
noise_h_co = highcutoff;

stepsperoct=2;
numoct=2;
f0 = round(octavesteps(250,stepsperoct,stepsperoct*numoct+1));

rms_scale_factor = 100;
tone_adjustment=1;

% Amplitude adjustment params
wantedDB = 70;
BenwaredBrms1 = 94;
fileFormat = 'wav';

% Head transfer function params
% htfCIPICPath = 'D:\Work\Code\KingLab\VAS_test\ov2_to_cipic\cipic_hrtfs';
% htfName = 'HRTF_9965.mat';
htfCIPICPath = 'D:\Work\Code\KingLab\VAS_test\CIPIC_hrtf_database\standard_hrir_database\subject_003';
htfName = 'hrir_final.mat';
htf = load(fullfile(htfCIPICPath,htfName));
htf.name = htfName;

htf.AzEl = [-90 0; -45 0; 0 0; 45 0; 90 0]; % Required azimuth and elevation angles (degree). [azimuth(1) elevation(1); azimuth(1) elevation(1); ...]


% Stage 2 - tones with high harmonics
%
% 3000-7200 Hz
% tone 1: 300 Hz F0, harmonics 10-24; tone 2: 600 Hz, harmonics 5-12

l_co = midcutoff;
u_co = highcutoff;

tone_hh = [];
for ii=1:length(f0),
    for h = l_co/f0(ii) : u_co/f0(ii)
        harm = tone(f0(ii)*h,stim_dur_ms,0,sr);
        if size(tone_hh,1)<ii
            tone_hh(ii,:) = harm;
        else
            tone_hh(ii,:) = zadd(tone_hh(ii,:),harm);
        end
    end
    if ii==1,
        tone_hh(ii,:) = tone_hh(ii,:)/rms_scale_factor;
    else
        tone_hh(ii,:) = tone_hh(ii,:)/rms(tone_hh(ii,:))*rms(tone_hh(1,:));
    end
end
% snd.tone_hh = tone_hh;
dp_noise = pnoise(stim_dur_ms,noise_l_co,noise_h_co,-30,0,sr);

% ERBs for ferret: 180 Hz at 300 Hz, 225 Hz at 500 Hz, 250 Hz at 600 Hz, 300 Hz at 1000 Hz
%cf = 300;
%erb_width = 180;
%cf=600; erb_width=250;
%calibrating with level at a higher F0 is more conservative as filter there is
%apparently narrower ==> less noise power
cf = 1000;
erb_width = 300;
erb_low = (-erb_width + sqrt(erb_width^2 + 4*cf^2))/2;
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

sample_harmonic = tone(f0(1)*10,stim_dur_ms,0,sr);
rms_primary = rms(sample_harmonic)/rms_scale_factor;
noise_adjustment = 1/rms(dp_noise_filt) * rms_primary * 10^(-5/20);
dp_noise_adjusted = dp_noise*noise_adjustment*tone_adjustment;

tone_high_harm_stim=[];
for ii=1:length(f0),
    tone_high_harm_stim(ii,:) = zadd(tone_hh(ii,:), dp_noise_adjusted);
end
snd.tone_high_harm_stim = tone_high_harm_stim;

% Stage 1 - full bandwidth tones
%
% 600-7200 Hz
% tone 1: 300 Hz F0, harmonics 2-24; tone 2: 600 Hz, harmonics 1-12

l_co = lowcutoff;
u_co = highcutoff;

tone_ah = [];
tone_all_harm_stim=[];
for ii=1:length(f0),
    for h = l_co/f0(ii) : u_co/f0(ii)
        harm = tone(f0(ii)*h,stim_dur_ms,0,sr);
        if size(tone_ah,1)<ii
            tone_ah(ii,:) = harm;
        else
            tone_ah(ii,:) = zadd(tone_ah(ii,:),harm);
        end
    end
    tone_ah(ii,:) = tone_ah(ii,:)/rms(tone_ah(ii,:))*rms(tone_hh(1,:));
    tone_all_harm_stim(ii,:) = zadd(tone_ah(ii,:), dp_noise_adjusted);
end
snd.tone_all_harm_stim = tone_all_harm_stim;
snd.tone_ah = tone_ah;


% Stage 3 - tones with lower harmonics
%
% 600-3000 Hz
% tone 1: 300 Hz F0, harmonics 2-10; tone 2: 600 Hz, harmonics 1-5

l_co = lowcutoff;
u_co = midcutoff;

tone_lh = [];
tone_low_harm_stim=[];
for ii=1:length(f0),
    for h = l_co/f0(ii) : u_co/f0(ii)
        harm = tone(f0(ii)*h,stim_dur_ms,0,sr);
        if size(tone_lh,1)<ii
            tone_lh(ii,:) = harm;
        else
            tone_lh(ii,:) = zadd(tone_lh(ii,:),harm);
        end
    end
    tone_lh(ii,:) = tone_lh(ii,:)/rms(tone_lh(ii,:))*rms(tone_hh(1,:));
    tone_low_harm_stim(ii,:) = zadd(tone_lh(ii,:), dp_noise_adjusted);
end
snd.tone_low_harm_stim = tone_low_harm_stim;

% Stage 4 - ALT phase tones with higher harmonics
%
% 600-3000 Hz
% tone 1: 300 Hz F0, harmonics 2-10; tone 2: 600 Hz, harmonics 1-5
l_co = midcutoff;
u_co = highcutoff;

tone_hh_alt = [];
tone_high_harm_alt_stim=[];
for ii=1:length(f0),
    for h = l_co/f0(ii) : u_co/f0(ii)
        if mod(h,2)==1,
            harm = tone(f0(ii)*h,stim_dur_ms,0,sr);
        else
            harm = tone(f0(ii)*h,stim_dur_ms,pi,sr);
        end
        if size(tone_hh_alt,1)<ii
            tone_hh_alt(ii,:) = harm;
        else
            tone_hh_alt(ii,:) = zadd(tone_hh_alt(ii,:),harm);
        end
    end
    tone_hh_alt(ii,:) = tone_hh_alt(ii,:)/rms(tone_hh_alt(ii,:))*rms(tone_hh(1,:));
    tone_high_harm_alt_stim(ii,:) = zadd(tone_hh_alt(ii,:), dp_noise_adjusted);
end
snd.tone_high_harm_alt_stim = tone_high_harm_alt_stim;

% Stage 5 - RAND phase tones with higher harmonics
%
% 600-3000 Hz
% tone 1: 300 Hz F0, harmonics 2-10; tone 2: 600 Hz, harmonics 1-5
l_co = midcutoff;
u_co = highcutoff;

tone_hh_rand = [];
tone_high_harm_rand_stim=[];

for ii=1:length(f0),
    for h = l_co/f0(ii) : u_co/f0(ii)
        thisphase=rand*2*pi;
    %     thisphase=rand;
        harm = tone(f0(ii)*h,stim_dur_ms,thisphase,sr);
        if size(tone_hh_rand,1)<ii
            tone_hh_rand(ii,:) = harm;
        else
            tone_hh_rand(ii,:) = zadd(tone_hh_rand(ii,:),harm);
        end
    end
    tone_hh_rand(ii,:) = tone_hh_rand(ii,:)/rms(tone_hh_rand(ii,:))*rms(tone_hh(1,:));
    tone_high_harm_rand_stim(ii,:) = zadd(tone_hh_rand(ii,:), dp_noise_adjusted);
end
snd.tone_high_harm_rand_stim = tone_high_harm_rand_stim;

% Just pure tones
tone_pure = [];
tone_pure_stim=[];
for ii=1:length(f0),
    tone_pure(ii,:) = tone(f0(ii),stim_dur_ms,0,sr);
    tone_pure(ii,:) = tone_pure(ii,:)/rms(tone_pure(ii,:))*rms(tone_hh(1,:));
    tone_pure_stim(ii,:) = zadd(tone_pure(ii,:), dp_noise_adjusted);
end
snd.tone_pure_stim = tone_pure_stim;
snd.tone_pure = tone_pure;

% Ramp
gatelength=5; %ramp time in msec
for ii=1:length(f0)
    snd.tone_high_harm_rand_stim(ii,:)=linear_envelope(snd.tone_high_harm_rand_stim(ii,:),gatelength,sr);
    snd.tone_high_harm_alt_stim(ii,:)=linear_envelope(snd.tone_high_harm_alt_stim(ii,:),gatelength,sr);
    snd.tone_low_harm_stim(ii,:)=linear_envelope(snd.tone_low_harm_stim(ii,:),gatelength,sr);
    snd.tone_high_harm_stim(ii,:)=linear_envelope(snd.tone_high_harm_stim(ii,:),gatelength,sr);
    snd.tone_all_harm_stim(ii,:)=linear_envelope(snd.tone_all_harm_stim(ii,:),gatelength,sr);
    snd.tone_pure_stim(ii,:)=linear_envelope(snd.tone_pure_stim(ii,:),gatelength,sr);
    snd.tone_pure(ii,:)=linear_envelope(snd.tone_pure(ii,:),gatelength,sr);
    snd.tone_ah(ii,:)=linear_envelope(snd.tone_ah(ii,:),gatelength,sr);
end

% Azimuth
% snd = structfun(@(x)(permute(repmat(x,[1 1 2]),[1 3 2])),snd,'UniformOutput',false);
ind = [];
for i = 1:size(htf.AzEl,1), % for every azimuth/elevation pairs
%     flt_l = getNearestUCDpulse(htf.AzEl(i,1),htf.AzEl(i,2),htf.hrir_l);
%     [flt_r, htf.RealAzEl(i,1), htf.RealAzEl(i,2)] = getNearestUCDpulse(htf.AzEl(i,1),htf.AzEl(i,2),htf.hrir_r);
    
    [flt_l, flt_r, ind] = getVASfromJan(htf.AzEl(i,1),htf.AzEl(i,2),ind);
    flt_l = resample(flt_l,round(sr),80000);
    flt_r = resample(flt_r,round(sr),80000);
    
    snd_l(i) = structfun(@(x)(convM(x,flt_l)),snd,'UniformOutput',false);
    snd_r(i) = structfun(@(x)(convM(x,flt_r)),snd,'UniformOutput',false);
end


% Adjust the amplitude to match the wantedDB var;
fList = fieldnames(snd_l);
maxVal = [];
for i = 1:length(f0), % F0
        for j = 1:size(htf.AzEl,1) % Azimuths
            for k = 1:length(fList), % Sound types
                snd_l(j).(fList{k})(i,:) = adjustRMSbenware( snd_l(j).(fList{k})(i,:),wantedDB,BenwaredBrms1);
                snd_r(j).(fList{k})(i,:) = adjustRMSbenware( snd_r(j).(fList{k})(i,:),wantedDB,BenwaredBrms1);
                maxVal = max([maxVal max(abs(snd_l(j).(fList{k})(i,:))) max(abs(snd_r(j).(fList{k})(i,:)))]);
            end
        end
end

msg = 'Amplitude above 1 or bellow -1. WAV files will be clipped. Use fileFormat = f32 instead.';
if strcmpi(fileFormat,'wav'); 
    if maxVal>1,
        warning(msg);
    end
end
disp('done')

%% Save sound sequence for Benware V2 04/04/2016
% Juste save 1 stim in one file + silence, Benware takes care of the
% repeats and the randomization. Makes analyzing easier.

% soundtypes=1:6; % 1 tone, 2 all harm, 3 high, 4 low, 5 alt, 6 rand
% soundNames = {'tone', 'all_harm', 'high', 'low', 'alt', 'rand'};

soundtypes=[1 2]; % 1 tone, 2 all harm, 3 high, 4 low, 5 alt, 6 rand
soundNames = {'tone', 'all_harm', 'high', 'low', 'alt', 'rand'};
stimNames = {'tone_pure_stim', 'tone_all_harm_stim', 'tone_high_harm_stim', 'tone_low_harm_stim', 'tone_high_harm_alt_stim', 'tone_high_harm_rand_stim'};
silDur = 0.5;
silVect = zeros(1,round(silDur*sr));
for fz = 1:length(f0),
    for stim = soundtypes,
        for azimuth = 1:size(htf.AzEl,1),
            clear stimulus
            stimulus(:,1) = [snd_l(azimuth).(stimNames{stim})(fz,:) silVect];
            stimulus(:,2) = [snd_r(azimuth).(stimNames{stim})(fz,:) silVect];
            
            % Save
            filename=sprintf('QuentinPitchSounds2016_%s_%dHz_%ddB_%ddeg', soundNames{stim},f0(fz),wantedDB,htf.AzEl(azimuth,1));
            switch fileFormat
                case 'wav'
                    audiowrite([filename '.wav'],stimulus,round(sr));
                case 'f32'
                    fid = fopen([filename '.f32'], 'w');
                    fwrite(fid,stimulus, 'float32');
                    fclose(fid);
                otherwise
                    error('Invalid file format.');
            end
            disp([filename ' saved']);
        end
    end
end

%% restore paths

path(oldPath);

%% Sanity check - listen, plot wav & plot fft of one stim
stimField = 'tone_pure_stim'; % List in snd.()
AzElInd = 1; % List in htf.RealAzEl
F0Ind = 1; % List in f0

fprintf('Checking %s with azimuth %d deg, elevation %d deg, and F0 %d Hz\n',stimField,htf.RealAzEl(AzElInd,1),htf.RealAzEl(AzElInd,2),f0(F0Ind));

% listening to the stim
sig = [snd_l(AzElInd).(stimField)(F0Ind,:);...
    snd_r(AzElInd).(stimField)(F0Ind,:)];
sound(sig,sr)

% plot waveforms
figure;
plot(snd_l(AzElInd).(stimField)(F0Ind,:));
hold on;
plot(snd_r(AzElInd).(stimField)(F0Ind,:));
hold off
legend({'left' 'right'})
title('Waveforms')

% plot spectrums
figure;
plot_fft(snd_l(AzElInd).(stimField)(F0Ind,:),sr);
hold on;
plot_fft(snd_r(AzElInd).(stimField)(F0Ind,:),sr);
hold off
legend({'left' 'right'})


%% Sanity check 2 - listening to stim sequence

stimField = 'tone_pure_stim'; % List in snd.()
AzElInd = [1 2 3 4 5]; % List in htf.RealAzEl
F0Ind = 1; % List in f0

fprintf('Checking %s with F0 %d Hz, [Azimuth | Elevation] ',stimField,f0(F0Ind));
fprintf(' [%d | %d],',htf.RealAzEl(AzElInd,1),htf.RealAzEl(AzElInd,2));
fprintf('\n');

sig = [];
for i = 1:length(AzElInd)
    sig = [sig [snd_l(AzElInd(i)).(stimField)(F0Ind,:) zeros(1,round(0.3*sr));...
        snd_r(AzElInd(i)).(stimField)(F0Ind,:) zeros(1,round(0.3*sr))]];
end
sound(sig,sr)
