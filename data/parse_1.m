% parse tomic & bays similarity data

clear;

% similarities
% 4 angles with 2 pairs of 2
% first angle of pair (a and c) range from 0:pi/72:(2*pi - 2*pi/72)
% second angle (b and d) are always larger than their a or c
% and range from pi/72:pi/72:(2*pi-pi/72)
% not obvious exactly how they are sampled beyond that

T = readtable('oriSimData_clean.csv');

ds.participant = T{:, 6};
ds.response = 1 - T{:, 1};
ds.a = T{:, 2};
ds.b = T{:, 3};
ds.c = T{:, 4};
ds.d = T{:, 5};
[ds.nTrials, ~] = size(ds.participant);
ds.nParticipants = length(unique(ds.participant));
[ds.stimuli, ~, tmp] = unique([ds.a ds.b ds.c ds.d]);
tmp = reshape(tmp, ds.nTrials, []);
ds.aIdx = tmp(:, 1);
ds.bIdx = tmp(:, 2);
ds.cIdx = tmp(:, 3);
ds.dIdx = tmp(:, 4);
ds.nStimuli = length(ds.stimuli);

save tomicBays ds
