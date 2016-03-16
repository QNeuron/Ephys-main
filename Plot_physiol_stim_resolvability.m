% Use this to plot figures of the resolvable harmonics for ferrets for a
% number of different F0's
function Plot_physiol_stim_resolvability(f0,l_co,u_co)
% 
% Inputs required:
% f0 is a vector of F0 values (Hz), in increasing order
% 
% Optional input: 
% l_co - lowest harmonic, in Hz
% h_co - highest harmonic, in Hz


if nargin>1,
    lower=l_co;
    upper=u_co;
else
    lower = 200;
    upper = 30000;
end



% stepsperoct=4;
% numoct=2;
% f0 = round(octavesteps(250,stepsperoct,stepsperoct*numoct+1));

% close all hidden

for ii=1:length(f0)
    determine_resolvable_harmonics_ferret(f0(ii),ii,lower,upper);
    set(gcf,'position',[133 41 1402 420])
    pause()
end

figure()
plot(f0,ones(1,length(f0)),'*k')
box off

