%% Thurstone from choice
% small number of one-dimensional stimuli
% each trial presents two pairs of stimuli
% and asks "which pair is more similar"

% generate data then look at inference

clear; close all;
preLoad = false;
printFigures = true;

   % graphical model script
   modelDir = './';
   modelName = 'thurstoneFromChoice_3b';


dataList = {...
  % 'whoIsSimilar_A'; ...
   %  'whoIsSimilar_B'; ...
  %  'whoIsSimilar_C'; ...
   % 'whoIsSimilar_D'; ...
    'whoIsSimilar_E'; ...
  };

analysisList = {...
   'recovery'; ...
   'recoveryThurstone'; ...
   };

%% general constants
load pantoneColors pantone;
figureFlag = false;
CIbounds = [2.5 97.5];

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

      case 'whoIsSimilar_B'
         rng(1)
         muTruth = [0:0.05:0.75 0.85 1];
         sigmaTruth = 0.2;
         nTrials = 1e3;

         axisScale = 2;

       case 'whoIsSimilar_C'
         rng(1)
         muTruth = rand(71, 1);
         sigmaTruth = 0.05;
         nTrials = 1e3;

         axisScale = 2;

          case 'whoIsSimilar_D'
         rng(2)
         muTruth = 0:1/70:1;
         sigmaTruth = 0.05;
         nTrials = 5e2;

         axisScale = 1;

          case 'whoIsSimilar_E'
         rng(2)
         muTruth = 0:1/70:1;
         sigmaTruth = 0.05;
         nTrials = 1e3;

         axisScale = 1;

   end

   % generate simulated task data
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

   %% sampling from graphical model
   % which engine to use
   engine = 'jags';

   % parameters to monitor
 %  params = {'mu', 'sigma'};
   params = {'mu', 'sigma'};

   % MCMC properties
   nChains    = 8;     % number of MCMC chains
   nBurnin    = 1e3;   % number of discarded burn-in samples
   nSamples   = 1e3;   % number of collected samples
   nThin      = 5;     % number of samples between those collected
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
         'modules'         , {'vonmises', 'wfcombopack'}                                , ...
         'parallel'        , doParallel                                );
      fprintf('%s took %f seconds!\n', upper(engine), toc); % show timing

      % just convergent enough chains
      [keepChains, rHat] = findKeepChains(chains.mu_1, 2, 1.1);
      fields = fieldnames(chains);
      for i = 1:numel(fields)
         chains.(fields{i}) = chains.(fields{i})(:, keepChains);
      end

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
         case 'recovery'

            figureFlag = true;

            % constants
            fontSize = 18;

            F = figure; clf; hold on;
            setFigure(F, [0.2 0.2 0.4 0.4], '');

            xLim = [floor(min(mu)) ceil(max(mu))];
            set(gca, ...
               'xlim'       , xLim    , ...
               'xtick'      , 0:0.2:1   , ...
               'xticklabelrot', 0, ...
               'ylim'       , [0 1]    , ...
               'ytick'      , 0:0.2:1  , ...
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

         case 'recoveryThurstone'

            figureFlag = true;

            % constants
            fontSize = 18;
            clrNames = fieldnames(pantone);
            rng(4); clrNames = clrNames(randperm(numel(clrNames)));
            faceAlpha = 0.25;
            xRes = 1e2;
            thresh = 1e-2;
            lineWidth = 1;

            F = figure; clf; hold on;
            setFigure(F, [0.2 0.2 0.6 0.6], '');

            xLim = [floor(min(muTruth)) ceil(max(muTruth))];
            xLim = xLim + [-1 1]*axisScale*range(xLim)/2;
            xx = xLim(1):range(xLim)/xRes:xLim(end);
            subplot(2, 1, 1); cla; hold on;
            set(gca, ...
               'xlim'       , xLim    , ...
               'xtick'      , []   , ...
               'ylim'       , [0 1]    , ...
               'ycolor'      , 'none'  , ...
               'box'        , 'off'     , ...
               'tickdir'    , 'out'     , ...
               'layer'      , 'top'     , ...
               'ticklength' , [0.01 0]  , ...
               'layer'      , 'top'     , ...
               'clipping'   , 'off'     , ...
               'fontsize'   , fontSize  );
            xlabel('Generating', 'fontsize', fontSize+4);
            moveAxis(gca, [1 1 1 1], [0 0 0 0]);

            % anticipating non-constant sigmas
            maxYY = 0;
            for idx = 1:nStimuli
               yy{idx} = exp(-((xx - muTruth(idx)).^2)/(2*3.1415*sigmaTruth^2));
               maxYY = max(maxYY, max(yy{idx}));
            end
            for idx = 1:nStimuli
               clr = pantone.(sprintf('%s', clrNames{idx}));
               keep = find(yy{idx} > thresh);
               plot(xx(keep), yy{idx}(keep)/maxYY, '-', ...
                  'linewidth', lineWidth, ...
                  'color', clr);
               patch(xx(keep), yy{idx}(keep)/maxYY, '-', ...
                  'edgecolor', 'none', ...
                  'facecolor', clr, ...
                  'facealpha', faceAlpha);
            end

            xLim = [floor(min(mu)) ceil(max(mu))];
            xLim = xLim + [-1 1]*axisScale*range(xLim)/2;
            xx = xLim(1):range(xLim)/xRes:xLim(end);
            subplot(2, 1, 2); cla; hold on;
            set(gca, ...
               'xlim'       , xLim    , ...
               'xtick'      , []   , ...
               'ylim'       , [0 1]    , ...
               'ycolor'      , 'none'  , ...
               'box'        , 'off'     , ...
               'tickdir'    , 'out'     , ...
               'layer'      , 'top'     , ...
               'ticklength' , [0.01 0]  , ...
               'layer'      , 'top'     , ...
               'clipping'   , 'off'     , ...
               'fontsize'   , fontSize  );
            xlabel('Inferred', 'fontsize', fontSize+4);
            moveAxis(gca, [1 1 1 1], [0 0 0 0]);

            maxYY = 0;
            for idx = 1:nStimuli
               yy{idx} = exp(-((xx - mu(idx)).^2)/(2*3.1415*sigma^2));
               maxYY = max(maxYY, max(yy{idx}));
            end
            for idx = 1:nStimuli
               clr = pantone.(sprintf('%s', clrNames{idx}));
               keep = find(yy{idx} > thresh);
               plot(xx(keep), yy{idx}(keep)/maxYY, '-', ...
                  'linewidth', lineWidth, ...
                  'color', clr);
               patch(xx(keep), yy{idx}(keep)/maxYY, '-', ...
                  'edgecolor', 'none', ...
                  'facecolor', clr, ...
                  'facealpha',faceAlpha);
            end

            [gca, T] = suplabel(sprintf('%d stimuli, %d trials', nStimuli, nTrials), 't');
            set(T, 'fontsize', fontSize, 'fontweight', 'normal');


            % print
            if figureFlag && printFigures
               if ~isfolder('figures')
                  !mkdir figures
               end
               warning off;
               print(sprintf('figures/%s_%s_bias.png', modelName, dataName), '-dpng');
               print(sprintf('figures/%s_%s_bias.pdf', modelName, dataName), '-dpdf');
               print(sprintf('figures/%s_%s_bias.eps', modelName, dataName), '-depsc');
               warning on;
               figureFlag = false;
            end

      end

   end
end