%% Thurstone from choice
% small number of one-dimensional stimuli, on a circle
% each trial is wrt a target stimulus
% task presents two stimuli, one of which is target
% and asks "which one is the target"

% generate data then look at inference

clear; close all;
preLoad = false;
printFigures = false;

   % graphical model script
   modelDir = './';
   modelName = 'thurstoneFromChoice_2b';

dataList = {...
   'thisOrCircularThat_A'; ...
   %'thisOrCircularThat_B'; ...
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

      % task properties5
      case 'thisOrCircularThat_A'
         rng(4)
         muTruth = [0 1.2 1.5 3.1 -3.1];
         sigmaTruth = 0.25;
         nTrials = 1e3;

         axisScale = 1;

          case 'thisOrCircularThat_B'
         rng(3)
         muTruth = -3.1415:1/(2*3.1415):3.1415;
         sigmaTruth = 0.25;
         nTrials = 1e3;

         axisScale = 1;

   end

   % generate simulated task data
   nStimuli = length(muTruth);
   for t = 1:nTrials
      tmp = randperm(nStimuli);
      s(t) = tmp(1);
      if rand < 0.5
         a(t) = tmp(1); b(t) = tmp(2);
      else
         a(t) = tmp(2); b(t) = tmp(1);
      end
      xA(t) = vmrand(muTruth(a(t)), 1/sigmaTruth^2);
      xB(t) = vmrand(muTruth(b(t)), 1/sigmaTruth^2);
      xS(t) = vmrand(muTruth(s(t)), 1/sigmaTruth^2);
      dAS = min([abs(xA(t)-xS(t)), abs(xA(t)-xS(t) - 2*3.1415), abs(xA(t)-xS(t) + 2*3.1415)]);
      dBS = min([abs(xB(t)-xS(t)), abs(xB(t)-xS(t) - 2*3.1415), abs(xB(t)-xS(t) + 2*3.1415)]);
      if (dAS - dBS) < 0
         y(t) = 0;
      else
         y(t) = 1;
      end
      % [a(t) b(t) s(t)]
      %  muTruth([a(t) b(t) s(t)])
      % [xA(t) xB(t) xS(t)]
   end

   %% sampling from graphical model
   % which engine to use
   engine = 'jags';

   % parameters to monitor
   params = {'mu', 'sigma', 'xA', 'xB', 'xS', 'dAS', 'dBS', 'yP'};
  params = {'mu', 'sigma', 'yP'};

   % MCMC properties
   nChains    = 8;     % number of MCMC chains
   nBurnin    = 1e3;   % number of discarded burn-in samples
   nSamples   = 1e3;   % number of collected samples
   nThin      = 1;     % number of samples between those collected
   doParallel = 1;     % whether MATLAB parallel toolbox parallizes chains

   % assign MATLAB variables to the observed nodes
   data = struct(...
      'a'  , a      , ...
      'b' , b  , ...
      's'   , s       , ...
      'y'   , y   , ...
      'nStimuli', nStimuli, ...
      'nTrials', nTrials);

   % censoring initial values
   for t = 1:nTrials
      xSinit(t) = 0.5;
      if y(t) == 0
         xAinit(t) = 0.6; xBinit(t) = 1.7;
      else
         xBinit(t) = 0.6; xAinit(t) = 1.7;
      end
   end

   % generator for initialization
%   generator = @()struct('xA', xAinit, 'xB', xBinit, 'xS', xSinit);
   generator = @()struct('xAtmp', xAinit, 'xBtmp', xBinit, 'xStmp', xSinit);

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
         'modules'         , {'vonmises'}                                , ...
         'parallel'        , doParallel                                );
      fprintf('%s took %f seconds!\n', upper(engine), toc); % show timing

      % % just convergent enough chains
      % [keepChains, rHat] = findKeepChains(chains.mu_3, 2, 1.1);
      % fields = fieldnames(chains);
      % for i = 1:numel(fields)
      %    chains.(fields{i}) = chains.(fields{i})(:, keepChains);
      % end

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

   % % reflect in origin if that fits better
   % if norm(mu - muTruth, 2) > norm(-mu - muTruth, 2)
   %    mu = -mu;
   %    for idx = 1:nStimuli
   %       muBounds(idx, :) = -muBounds(idx, :);
   %    end
   % end

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

            %xLim = [floor(min(mu)) ceil(max(mu))];
            xLim = [-3.1415 3.1415];
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

            % xLim = [floor(min(muTruth)) ceil(max(muTruth))];
            % xLim = xLim + [-1 1]*axisScale*range(xLim)/2;
            xLim = [-3.1415 3.1415];

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

            % xLim = [floor(min(mu)) ceil(max(mu))];
            % xLim = xLim + [-1 1]*axisScale*range(xLim)/2;
            xLim = [-3.1415 3.1415];

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