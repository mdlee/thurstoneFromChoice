% parse tomic & bays similarity data

clear;

% similarities
% 4 angles with 2 pairs of 2
% first angle of pair (a and c) range from 0:pi/72:(2*pi - 2*pi/72)
% second angle (b and d) are always larger than their a or c
% and range from pi/72:pi/72:(2*pi-pi/72)
% not obvious exactly how they are sampled beyond that
pi = 3.1415926535;
decimalPoints = 4;

T = readtable('oriSimData_clean.csv');

ds.participant = T{:, 6};
ds.response = 1 - T{:, 1};
% dividing by 2 comes from personal correspondence with authors
ds.a = round(T{:, 2}, decimalPoints)/2;
ds.b = round(T{:, 3}, decimalPoints)/2;
ds.c = round(T{:, 4}, decimalPoints)/2;
ds.d = round(T{:, 5}, decimalPoints)/2;
[ds.nTrials, ~] = size(ds.participant);
ds.nParticipants = length(unique(ds.participant));
[ds.stimuli, ~, tmp] = unique([ds.a ds.b ds.c ds.d]);
tmp = reshape(tmp, ds.nTrials, []);
ds.aIdx = tmp(:, 1);
ds.bIdx = tmp(:, 2);
ds.cIdx = tmp(:, 3);
ds.dIdx = tmp(:, 4);
ds.nStimuli = length(ds.stimuli);

step = diff(ds.stimuli(1:2));
binsE = [(0-step/2):step:(pi+step/2) pi];
binsC = 0:step:pi;

T = readtable('oriPercData_clean.csv');
dp.participant = T{:, 4};
% dp.response = max(0, (round(T{:, 2}, decimalPoints) + pi)/2);
% dp.stimulus = max(0, (round(T{:, 3}, decimalPoints) + pi)/2);
dp.response = mod(round(T{:, 2}, decimalPoints), pi);
dp.stimulus = mod(round(T{:, 3}, decimalPoints), pi);
[dp.nTrials, ~] = size(dp.participant);
dp.nParticipants = length(unique(dp.participant));
dp.sIdx = zeros(dp.nTrials, 1);
for i = 1:dp.nTrials
   dp.sIdx(i) = find(histcounts(dp.stimulus(i), binsE));
   choices = dp.response(i) + [0 -pi pi];
   [~, chooseIdx] = min(abs(choices - dp.stimulus(i)));
   dp.response(i) = choices(chooseIdx);
  % [dp.stimulus(i) dp.sIdx(i) dp.response(i)]
end
dp.sIdx(find(dp.sIdx == 73)) = 1; % half of these are near 0 and half near pi
dp.stimuli = ds.stimuli;
dp.nStimuli = ds.nStimuli;


T = readtable('oriWMData_clean.csv');
dm.response = mod(round(T{:, ExcelColNo('H')}, decimalPoints), pi);
dm.target = mod(round(T{:, ExcelColNo('I')}, decimalPoints), pi);
dm.participant = T{:, ExcelColNo('O')};
[dm.nTrials, ~] = size(dm.participant);
dm.nontarget = mod(round(T{:, ExcelColNo('J'):ExcelColNo('N')}, decimalPoints), pi);
dm.nNontargets = 5 - sum(isnan(dm.nontarget), 2);
dm.setSize = T{:, ExcelColNo('G')}; % same info as nNonTargets
dm.tIdx = zeros(dm.nTrials, 1);
for i = 1:dm.nTrials
   dm.tIdx(i) = find(histcounts(dm.target(i), binsE));
   choices = dm.response(i) + [0 -pi pi];
   [~, chooseIdx] = min(abs(choices - dm.target(i)));
   dm.response(i) = choices(chooseIdx);
end
dm.tIdx(find(dm.tIdx == 73)) = 1; % half of these are near 0 and half near pi
dm.nIdx = nan(dm.nTrials, 5);
for i = 1:dm.nTrials
   for j = 1:5
      if ~isnan(dm.nontarget(i, j))
         dm.nIdx(i, j) = find(histcounts(dm.nontarget(i, j), binsE));
      end
   end
end
dm.nIdx(find(dm.nIdx == 73)) = 1; % half of these are near 0 and half near pi
dm.stimuli = ds.stimuli;
dm.nStimuli = ds.nStimuli;


save tomicBays ds dp dm



