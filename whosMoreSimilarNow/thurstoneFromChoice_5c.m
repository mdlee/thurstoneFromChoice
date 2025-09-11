%% Thurstone from choice
% small number of one-dimensional stimuli
% each trial presents two pairs of stimuli
% and asks "which pair is more similar"

% generate data then look at inference

clear; close all;
preLoad = true;
printFigures = true;

% graphical model script
modelDir = './';
modelName = 'thurstoneFromChoice_5c';


dataList = {...
   %'whoIsSimilar_A'; ...
   'tomicBays'; ...
   };

analysisList = {...
   'recoveryScatter'; ...
   'recoveryCircle'; ...
   };

%% general constants
load pantoneColors pantone;
figureFlag = false;
CIbounds = [2.5 97.5];
pi = 3.1415;

% loop over data
for dataIdx = 1:numel(dataList)
   dataName = dataList{dataIdx};
   switch dataName

      % task properties
      case 'whoIsSimilar_A'
         rng(1)
         muTruth = [0.1 0.3 0.5 0.8 0.85];
         sigmaTruth = 0.1;
         nTrials = 5e2;

         axisScale = 0.5;

         generateData = true;

      case 'tomicBays'
         dataDir = 'data/';
         dataName = 'tomicBays';
         load([dataDir dataName], 'ds');


         % ds.a = mod(ds.a, pi);
         % ds.b = mod(ds.b, pi);
         % ds.c = mod(ds.c, pi);
         % ds.d = mod(ds.d, pi);
         % 
         a = ds.aIdx;
         b = ds.bIdx;
         c = ds.cIdx;
         d = ds.dIdx;
         y = ds.response;
         nTrials = ds.nTrials;
         nStimuli = ds.nStimuli;

         muTruth = mod(ds.stimuli, 2*pi);
         sigmaTruth = 0.01;

         axisScale = 2;

         generateData = false;


   end

   % generate simulated task data
   if generateData
      nStimuli = length(muTruth);
      for t = 1:nTrials
         tmp = randperm(nStimuli);
         a(t) = tmp(1);
         b(t) = tmp(2);
         c(t) = tmp(3);
         d(t) = tmp(4);
         xA(t) = randn*sigmaTruth+muTruth(a(t));
         xB(t) = randn*sigmaTruth+muTruth(b(t));
         xC(t) = randn*sigmaTruth+muTruth(c(t));
         xD(t) = randn*sigmaTruth+muTruth(d(t));
         if (abs(xA(t)-xB(t)) - abs(xC(t)-xD(t))) < 0
            y(t) = 0;
         else
            y(t) = 1;
         end
      end
   end

   %% sampling from graphical model
   % which engine to use
   engine = 'jags';

   % parameters to monitor
   params = {'mu', 'sigma', 'deviance'};

   % MCMC properties
   nChains    = 8;     % number of MCMC chains
   nBurnin    = 1e2;   % number of discarded burn-in samples
   nSamples   = 1e2;   % number of collected samples
   nThin      = 1;     % number of samples between those collected
   doParallel = 1;     % whether MATLAB parallel toolbox parallizes chains

   % assign MATLAB variables to the observed nodes
   data = struct(...
      'a'  , a      , ...
      'b' , b  , ...
      'c'   , c       , ...
      'd'   , d , ...
      'y'   , y   , ...
      'nStimuli', nStimuli, ...
      'nTrials', nTrials);

   % censoring initial values
   for t = 1:nTrials
      if y(t) == 0
         xAinit(t) = 0.6; xBinit(t) = 0.7;
         xCinit(t) = 0.6; xDinit(t) = 0.8;
      else
         xAinit(t) = 0.6; xBinit(t) = 0.8;
         xCinit(t) = 0.6; xDinit(t) = 0.7;
      end
   end

   % generator for initialization
   generator = @()struct('xA', xAinit, 'xB', xBinit, ...
      'xC', xCinit, 'xD', xDinit);

   fileName = sprintf('%s_%s_%s.mat', modelName, dataName, engine);

   if preLoad && isfile(sprintf('storage/%s', fileName))
      fprintf('Loading pre-stored samples for model %s on data %s\n', modelName, dataName);
      load(sprintf('storage/%s', fileName), 'chains', 'stats', 'diagnostics', 'info');
   else
      tic; % start clock
      [stats, chains, diagnostics, info] = callbayes(engine, ...
         'model'           , sprintf('%s/%s_%s.txt', modelDir, modelName, engine)   , ...   , ...
         'data'            , data                                      , ...
         'outputname'      , 'samples'                                 , ...
         'init'            , generator                                 , ...
         'datafilename'    , modelName                                 , ...
         'initfilename'    , modelName                                 , ...
         'scriptfilename'  , modelName                                 , ...
         'logfilename'     , sprintf('/tmp/%s', modelName)              , ...
         'nchains'         , nChains                                   , ...
         'nburnin'         , nBurnin                                   , ...
         'nsamples'        , nSamples                                  , ...
         'monitorparams'   , params                                    , ...
         'thin'            , nThin                                     , ...
         'workingdir'      , sprintf('/tmp/%s', modelName)              , ...
         'verbosity'       , 0                                         , ...
         'saveoutput'      , true                                      , ...
         'modules'         , {'dic'}                                , ...
         'parallel'        , doParallel                                );
      fprintf('%s took %f seconds!\n', upper(engine), toc); % show timing

      

      % convergence of each parameter
      disp('Convergence statistics:')
      grtable(chains, 1.05)

      % basic descriptive statistics
      disp('Descriptive statistics for all chains:')
      codatable(chains);

      fprintf('Saving samples for model %s on data %s\n', modelName, dataName);
      if ~isfolder('storage')
         !mkdir storage
      end
      save(sprintf('storage/%s', fileName), 'chains', 'stats', 'diagnostics', 'info');

   end

   % just convergent enough chains
      [keepChains, rHat] = findKeepChains(chains.mu_3, 2, 1.1);
      fields = fieldnames(chains);
      for i = 1:numel(fields)
         chains.(fields{i}) = chains.(fields{i})(:, keepChains);
      end

   %% common analysis
   mu = codatable(chains, 'mu', @mean);
   sigma = codatable(chains, 'sigma', @mean);
   muBound = nan(nStimuli, 2);
   for idx = 1:nStimuli
      muBounds(idx, :) = prctile(chains.(sprintf('mu_%d', idx))(:), CIbounds);
   end

   % reflect in origin if that fits better
   if sum((mu - muTruth').^2) > sum(((1-mu) - muTruth').^2)
      mu = 1-mu;
      for idx = 1:nStimuli
         muBounds(idx, :) = 1-muBounds(idx, :);
      end
   end

   %% specific analyses
   for analysisIdx = 1:numel(analysisList)
      analysis = analysisList{analysisIdx};

      switch analysis
         case 'recoveryScatter'

            figureFlag = true;

            % constants
            fontSize = 18;

            F = figure; clf; hold on;
            setFigure(F, [0.2 0.2 0.4 0.4], '');

            xLim = [floor(min(mu)) ceil(max(mu))];
            set(gca, ...
               'xlim'       , [0 2*pi]    , ...
               'xtick'      , [0 pi 2*pi]   , ...
               'xticklabelrot', 0, ...
               'ylim'       , [0 2*pi]    , ...
               'ytick'      , [0 pi 2*pi]  , ...
               'box'        , 'off'     , ...
               'tickdir'    , 'out'     , ...
               'layer'      , 'top'     , ...
               'ticklength' , [0.01 0]  , ...
               'layer'      , 'top'     , ...
               'clipping'   , 'off'     , ...
               'fontsize'   , fontSize  );
            axis square;
            xlabel('Inferred', 'fontsize', fontSize+4);
            ylabel('Generating', 'fontsize', fontSize+4);
            xtickformat('%.2f');
            ytickformat('%.2f');
            moveAxis(gca, [1 1 1 1], [0 0.025 0 0]);
            axis square;
            Raxes(gca, 0.01, 0.01);

            for idx = 1:nStimuli
               plot(muBounds(idx, :), muTruth(idx)*ones(1, 2), '-', ...
                  'color', pantone.ClassicBlue);
               plot(mu(idx), muTruth(idx), 'o', ...
                  'markerfacecolor', pantone.ClassicBlue, ...
                  'markeredgecolor', 'w', ...
                  'markersize', 6);
            end

         case 'recoveryCircle'

              figureFlag = true;

            % constants
            fontSize = 18;
            clrNames = fieldnames(pantone);

            F = figure; clf; hold on;
            setFigure(F, [0.2 0.2 0.4 0.4], '');

            %xLim = [floor(min(mu)) ceil(max(mu))];
            xLim = [-1 1];
            set(gca, ...
               'xlim'       , xLim    , ...
               'xtick'      , []   , ...
               'xticklabelrot', 0, ...
               'ylim'       , xLim    , ...
               'ytick'      , []  , ...
               'box'        , 'off'     , ...
               'tickdir'    , 'out'     , ...
               'layer'      , 'top'     , ...
               'ticklength' , [0.01 0]  , ...
               'layer'      , 'top'     , ...
               'clipping'   , 'off'     , ...
               'fontsize'   , fontSize  );
            axis square;
            % xlabel('Inferred', 'fontsize', fontSize+4);
            % ylabel('Generating', 'fontsize', fontSize+4);
            % moveAxis(gca, [1 1 1 1], [0 0.025 0 0]);
            % axis square;
            % Raxes(gca, 0.01, 0.01);
            axis off
r = 1.1;
            for idx = 1:nStimuli
               plot(cos(muTruth(idx)), sin(muTruth(idx)), 'o', ...
                  'markerfacecolor',  pantone.(sprintf('%s', clrNames{idx})), ...
                  'markeredgecolor', 'w');
               plot(r*cos(mu(idx)), r*sin(mu(idx)), '+', ...
                  'color', pantone.(sprintf('%s', clrNames{idx})), ...
                  'markersize', 8, ...
                  'linewidth', 2);
                plot([cos(muTruth(idx)) r*cos(mu(idx))], [sin(muTruth(idx)) r*sin(mu(idx))], '-', ...
                  'color', pantone.GlacierGray);
                pause;
            end

      end

      % print
            if figureFlag && printFigures
               if ~isfolder('figures')
                  !mkdir figures
               end
               warning off;
               print(sprintf('figures/%s_%s_%s.png', analysis, modelName, dataName), '-dpng');
              % print(sprintf('figures/%s_%s_%s.pdf',  analysis, modelName, dataName), '-dpdf');
               print(sprintf('figures/%s_%s_%s.eps',  analysis, modelName, dataName), '-depsc');
               warning on;
               figureFlag = false;
            end
   end
end