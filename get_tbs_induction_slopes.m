%% induction slopes for theta burst experiments 
clear all 
close all
fpath_variables = 'F:\Google Drive\Work\Research Projects\Theta LTP\Matlab Variables\'; % variables

load(strcat(fpath_variables,'slopes.mat')); % slopes

% get average slopes during tetanus for control condition
x_all = zeros(length(slopes{1,1,1,1,1}(1).indSlopes),length(slopes{1,1,1,1,1}));
for i = 1: length(slopes{1,1,1,1,1})
    x_all(:,i) = slopes{1,1,1,1,1}(i).indSlopes/slopes{1,1,1,1,1}(i).indSlopes(1);
end

x = mean(x_all, 2);

burstf = 100;
interf = 5;
pulses=4;
bursts=15;

t_burst = 0:1/burstf:(pulses-1)/burstf;
t = zeros(bursts*pulses,1);
for i=1:bursts
    t((i-1)*pulses+1:i*pulses) = t_burst + (i-1)/interf;
end

figure;plot(t, x)

clearvars -except keppVariables 't' 'x'


