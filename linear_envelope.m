function stimulus = linear_envelope(stimulus,gatelength,samplefreq)
% Inputs:
% stimulus           = waveform input to ramp on and off
% gatelength         = duration of raised-cosine gates applied to
%                            onset and offset (milliseconds)
% samplefreq         = Hz
gatelength_samples = round(gatelength*samplefreq/1000);  %gate length is in msecs
ramp = ones(size(stimulus));
ramp(1:gatelength_samples) = linspace(0,1,gatelength_samples);
ramp((end-gatelength_samples+1):end) = fliplr(linspace(0,1,gatelength_samples));
stimulus = stimulus.*ramp;
