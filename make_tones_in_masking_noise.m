

sr = 24414;
stim_dur_ms = 200;

noise_l_co = 100;
noise_h_co = 7500;

f0_1 = 300;
f0_2 = 500;

l_co = 3000;
u_co = 7500;

rms_scale_factor = 100;

tone1 = [];
for h = l_co/f0_1 : u_co/f0_1
    harm = tone(f0_1*h,stim_dur_ms,0,sr);
    tone1 = zadd(tone1,harm);
end
tone1 = tone1/rms_scale_factor;

tone1_full = [];
for h = 1 : floor((sr/2)/f0_1)
    harm = tone(f0_1*h,stim_dur_ms,0,sr);
    tone1_full = zadd(tone1_full,harm);
end
tone1_full = tone1_full/rms(tone1_full)*rms(tone1);

tone2 = [];
for h = l_co/f0_2 : u_co/f0_2
    harm = tone(f0_2*h,stim_dur_ms,0,sr);
    tone2 = zadd(tone2,harm);
end
tone2 = tone2/rms(tone2)*rms(tone1); %level difference was ~2dB ==> primaries of tone2 are 2dB higher than those of tone1
%this should be offset by the greater number of harmonics in tone1, leading
%to approximately equal amplitude DPs in the two stimuli
%we assume the DPs to be about 15 dB below the level of the primaries

%determine factor tone_adjustment by which tones must be adjusted to yield 65 dB SPL
%through Kerry's system

tone_adjustment = 1;

tone1 = tone1*tone_adjustment;
tone2 = tone2*tone_adjustment;


dp_noise = pnoise(stim_dur_ms,noise_l_co,noise_h_co,-30,0,sr);

% ERBs for ferret: 180 Hz at 300 Hz, 210 Hz at 500 Hz
%cf = 300;
%erb_width = 180;
%calibrating with level at 500 Hz is more conservative as filter there is
%apparently narrower ==> less noise power
cf = 500;
erb_width = 210;
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

sample_harmonic = tone(f0_1*10,stim_dur_ms,0,sr);

rms_primary = rms(sample_harmonic)/rms_scale_factor;

noise_adjustment = rms_primary*10^(-5/20) / rms(dp_noise_filt);

dp_noise_adjusted = dp_noise*tone_adjustment*noise_adjustment;




