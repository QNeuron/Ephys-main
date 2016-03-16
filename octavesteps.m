function freqs=octavesteps(InitialFreq,StepsPerOctave,numSteps)

%jkb feb 2006 (from Bronaugh)
%
%freqs=octavesteps(Sf,StepSize,n)
%f1 = starting frequency
%StepsPerOctave = number of steps per octave 
%m is the total number of steps required (negative values for lower freqs)
%returns freqs, a vector of the required frequencies

StepSize=1/StepsPerOctave;
f(1)=InitialFreq;
for ii=1:numSteps-1
    f(ii+1)=f(ii)*2^StepSize;
end
freqs=f;