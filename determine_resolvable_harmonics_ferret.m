function determine_resolvable_harmonics_ferret(f0,fig,l_co,u_co)

% KMMW Aug 2013; updated March 2016
% Plots harmonics of the F0 given as input (f0), indicated which harmonics
% are resolved in the response of the ferret auditory nerve (based on
% Sumner & Palmer 2012).
% 
% Required input:
% f0 - the fundamental frequency of the sound, in Hz
% 
% Optional input: 
% fig - figure number
% l_co - lowest harmonic, in Hz
% h_co - highest harmonic, in Hz


if nargin>3,
    l_co=l_co;
    u_co=u_co;
else
    l_co = 200;
    u_co = 30000;
end

all_cf=[];
all_lowcf=[];
all_highcf=[];
all_erb=[];
all_res_low=[];
all_res_high=[];

ind=1;
if nargin>1
    figure(fig);
else
    figure();
end
clf;
hold on


for h = l_co/f0 : u_co/f0
    if h>=1
        
        %harmonic tone frequency
        cf=f0*h;
        all_cf=[all_cf cf];

        %find bandwidth for this harmonic component, based on ERB in
        %auditory nerve recordings
        erb = 310*((cf/1000).^0.533); % from Sumner & Palmer 2012
        all_erb=[all_erb erb];
        lowcf = (-erb + sqrt(erb^2 + 4*cf^2))/2; %lower edge of rectangular band
        highcf = (erb + sqrt(erb^2 + 4*cf^2))/2; %lower edge of rectangular band
        all_lowcf=[all_lowcf lowcf];
        all_highcf=[all_highcf highcf];


        % plot harmonic center frequency and bandwidth
        x=logspace(log(lowcf)/log(10),log(cf)/log(10),100);
        y=linspace(0,1,100).^.1;
        plot(x,y,'-k');
        x=logspace(log(cf)/log(10),log(highcf)/log(10),100);
        y=linspace(1,0,100).^.1;
        plot(x,y,'-k');
        plot(cf,1,'ob')

        % is this harmonic resolved?
        if ind>1
            all_res_low=[all_res_low lowcf>all_highcf(ind-1)]; % is lower edge of the band resolved?
            all_res_high=[all_res_high all_highcf(ind-1)<lowcf]; % is upper edge of the band resolved?
        else
            all_res_low=[all_res_low 1];
        end
        ind=ind+1;
    end
end
all_res_high=[all_res_high 0]; %setting final upper band to be resolved (will not matter if lower edge is unresolved)

% plot tiddy-up
set(gca,'xscale','log')
xlabel('Frequency (Hz)')
ylabel('Norm energy')
xlim([100 25000])
ylim([0 1.05])
box off
title(['F0 = ' num2str(f0) ' Hz'])

% indicate resolved harmonics with red asterisk
all_res = and(all_res_low,all_res_high);
f=find(all_res);
plot(all_cf(f),ones(1,length(f)).*.5,'*r');
disp('Harmonics (black curves), with resolved harmonics indicated (red asterisks)')


