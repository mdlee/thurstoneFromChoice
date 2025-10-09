%% Thurstone from choice
% weighted foil memory task

clear; close all;
preLoad = true;
printFigures = true;

% graphical model script
modelDir = './';
modelName = 'thurstoneFromChoice_8c';


dataList = {...
   'tomicBaysMemory'; ...
   % 'tomicBaysLimitedTrials'; ...
   };

analysisList = {...
   % 'recoveryScatter'; ...
   % 'recoveryCircle'; ...
    'bothPanels'; ...
   %'distanceByProbability'; ...
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
  

      case 'tomicBaysMemory'
         dataDir = 'data/';
         dataName = 'tomicBays';
         load([dataDir dataName], 'dm');
         dataName = 'tomicMemory';

         y = dm.response;
         nTrials = dm.nTrials;
         nStimuli = dm.nStimuli;
         target = dm.tIdx;
         foil = dm.nontarget;
         foil(find(isnan(foil))) = 0;
         nFoils = dm.setSize-1;

         muTruth = dm.stimuli;
         sigmaTruth = 0.01;

         axisScale = 2;
   
   end


   %% sampling from graphical model
   % which engine to use
   engine = 'jags';

   % parameters to monitor
   params = {'mu', 'sigma', 'omega2', 'omega5', 'tau', 'yp'};

   % MCMC properties
   nChains    = 8;     % number of MCMC chains
   nBurnin    = 2e3;   % number of discarded burn-in samples
   nSamples   = 2e3;   % number of collected samples
   nThin      = 10;     % number of samples between those collected
   doParallel = 1;     % whether MATLAB parallel toolbox parallizes chains

   % %    % MCMC properties
   nChains    = 8;     % number of MCMC chains
   nBurnin    = 1e3;   % number of discarded burn-in samples
   nSamples   = 1e3;   % number of collected samples
   nThin      = 1;     % number of samples between those collected
   doParallel = 1;     % whether MATLAB parallel toolbox parallizes chains

   % assign MATLAB variables to the observed nodes
   data = struct(...
      'target'  , target      , ...
      'foil', foil , ...
      'nFoils', nFoils, ...
      'y'   , y   , ...
      'nStimuli', nStimuli, ...
      'nTrials', nTrials);

   % generator for initialization
   generator = @()struct('xS', rand(nTrials, 1));%, 'mu', dm.stimuli);

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
         'modules'         , {}                                , ...
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
      save(sprintf('storage/%s', fileName), 'chains', 'stats', 'diagnostics', 'info', '-v7.3');

   end
   %
   % % just convergent enough chains
   % [keepChains, rHat] = findKeepChains(chains.sigma, 2, 1.1);
   % fields = fieldnames(chains);
   % for i = 1:numel(fields)
   %    chains.(fields{i}) = chains.(fields{i})(:, keepChains);
   % end

   %% common analysis
   mu = codatable(chains, 'mu', @mean);
   sigma = codatable(chains, 'sigma', @mean)
   omega2 = codatable(chains, 'omega2', @mean)
   omega5 = codatable(chains, 'omega5', @mean)
   tau = codatable(chains, 'tau', @mean)
   % yp = codatable(chains, 'yp', @mean);
   muBound = nan(nStimuli, 2);
   % need to calculate these CIs on the circle!
   % any samples further than pi/2 from the mean get pi added to them
   for idx = 1:nStimuli
      muBounds(idx, :) = prctile(chains.(sprintf('mu_%d', idx))(:), CIbounds);
   end
   % 
   % % % reflect in origin if that fits better
   % % if sum((mu - muTruth').^2) > sum(((1-mu) - muTruth').^2)
   % %    mu = 1-mu;
   % %    for idx = 1:nStimuli
   % %       muBounds(idx, :) = 1-muBounds(idx, :);
   % %    end
   % % end
   % 
   % 
   % % log likelihood
   % LL = sum(y.*log(yp)) + sum((1-y).*log(1-yp));
   % fprintf('Log-likelihood = %1.4f\n', LL);

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
               'xlim'       , [0 pi]    , ...
               'xtick'      , [0 pi/2 pi]   , ...
               'xticklabelrot', 0, ...
               'ylim'       , [0 pi]    , ...
               'ytick'      , [0 pi/ 2 pi]  , ...
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
               %   pause;
            end

         case 'bothPanels'


            figureFlag = true;

            % constants
            fontSize = 24;

            F = figure; clf; hold on;
            setFigure(F, [0.2 0.2 0.7 0.5], '');

            xLim = [floor(min(mu)) ceil(max(mu))];
            subplot(1, 2, 1); cla; hold on;
            set(gca, ...
               'xlim'       , [0 pi]    , ...
               'xtick'      , [0 pi/4 pi/2 3*pi/4 pi]   , ...
               'xticklabelrot', 0, ...
               'xticklabel' , {'$0$', '$\frac{\pi}{4}$', '$\frac{\pi}{2}$', '$\frac{3\pi}{4}$', '$\pi$'}, ...
               'ylim'       , [0 pi]    , ...
               'ytick'      , [0 pi/4 pi/2 3*pi/4 pi]   , ...
               'yticklabel' , {'$0$', '$\frac{\pi}{4}$', '$\frac{\pi}{2}$', '$\frac{3\pi}{4}$', '$\pi$'}, ...
               'ticklabelinterpreter', 'latex', ...
               'box'        , 'off'     , ...
               'tickdir'    , 'out'     , ...
               'layer'      , 'top'     , ...
               'ticklength' , [0.02 0]  , ...
               'layer'      , 'top'     , ...
               'clipping'   , 'off'     , ...
               'fontsize'   , fontSize  );
            axis square;
            % grid on;
            xlabel('Psychological', 'fontsize', fontSize);
            ylabel('Physical', 'fontsize', fontSize);
            text(-1, 3.5, 'A', 'fontsize', fontSize, 'fontweight', 'bold');
            % xtickformat('%.2f');
            % ytickformat('%.2f');
            moveAxis(gca, [1 1 0.7 0.7], [0 0.2 0 0]);
            Raxes(gca, 0.02, 0.01);

            for i = pi/4:pi/4:3*pi/4
               plot([i i], [0 pi], '-', ...
                  'color', pantone.GlacierGray);
               plot([0 pi], [i i], '-', ...
                  'color', pantone.GlacierGray);
            end


            for idx = 1:nStimuli
               plot(muBounds(idx, :), muTruth(idx)*ones(1, 2), '-', ...
                  'color', pantone.ClassicBlue, ...
                  'linewidth', 1);
               plot(mu(idx), muTruth(idx), 'o', ...
                  'markerfacecolor', pantone.ClassicBlue, ...
                  'markeredgecolor', 'w', ...
                  'linewidth', 0.5, ...
                  'markersize', 4);

            end
            plot([0 pi], [0 pi], '-', ...
               'color', pantone.AuroraRed, 'linewidth', 0.5);

            xLim = [-1 1];
            subplot(1, 2, 2); cla; hold on;
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
            moveAxis(gca, [1 1 1.3 1.3], [-0.1 -0.3 0 0]);
            axis square;
            axis off
            text(-1.15, 1.185, 'B', 'fontsize', fontSize, 'fontweight', 'bold');

            %            r = [0.85 0.95];
            %  for idx = 1:nStimuli
            %         plot([cos(muTruth(idx)) r(2)*cos(mu(idx))], [sin(muTruth(idx)) r(2)*sin(mu(idx))], '-', ...
            %        'color', pantone.GlacierGray);
            % plot(cos(mu(idx)), sin(mu(idx)), 'o', ...
            %        'markerfacecolor',  'w', ...
            %        'markeredgecolor', 'k');
            %     plot(r*cos(muTruth(idx)), r*sin(muTruth(idx)), 'k-', ...
            %        'linewidth', 2);
            %  end

            r = [1.05 1.15];
            for idx = 1:nStimuli
               plot([cos(muTruth(idx)) r(1)*cos(mu(idx))], [sin(muTruth(idx)) r(1)*sin(mu(idx))], '-', ...
                  'color', pantone.GlacierGray);
               plot(cos(muTruth(idx)), sin(muTruth(idx)), 'o', ...
                  'markerfacecolor',  'w', ...
                  'markeredgecolor', 'k');
               plot(r*cos(mu(idx)), r*sin(mu(idx)), 'k-', ...
                  'linewidth', 2);
            end


         case 'distanceByProbability'

            figureFlag = true;

            fontSize = 20;
            dataClr = pantone.DuskBlue;
            modelClr = pantone.AuroraRed;
            markerSize = 6;

            lo = -pi/2; hi = pi/2; eps = 0.1;
            xTick = [lo 0 hi];
            xTickLabel = {'$-\frac{\pi}{2}$', '$0$', '$\frac{\pi}{2}$'};

            binsC = lo+eps/2:eps:hi-eps/2;
            binsE = lo:eps:hi;

            F = figure; clf; hold on;
            setFigure(F, [0.2 0.2 0.4 0.4], '');
            set(gca, ...
               'xlim'     , [lo hi]              , ...
               'xtick'        , xTick , ...
               'xticklabel', xTickLabel, ...
               'ticklabelinterpreter', 'latex', ...
               'ylim'       , [0 1]         , ...
               'ytick'      , 0:0.1:1   , ...
               'box'        , 'off'                 , ...
               'tickdir'    , 'out'                 , ...
               'ticklength' , [0.01 0]              , ...
               'fontsize'   , fontSize              );
            moveAxis(gca, [1 1 1 0.8], [0 0.075 0 0]);
            xlabel('Difference in Distance', 'fontsize', fontSize+2);
            Raxes(gca, 0.01, 0.01);

            countNum = zeros(size(binsC));
            countDen = zeros(size(binsC));
            for i = 1:length(binsC)
               mp{i} = [];
            end
            for t = 1:ds.nTrials
               dab = min(abs(ds.a(t)-ds.b(t)), pi - abs(ds.a(t)-ds.b(t)));
               dcd = min(abs(ds.c(t)-ds.d(t)), pi - abs(ds.c(t)-ds.d(t)));
               dist = dab - dcd;
               W  = histcounts(dist, binsE);
               % [ds.a(t) ds.b(t) ds.c(t) ds.d(t)]
               % d
               % find(W)
               % pause;
               countNum(find(W)) = countNum(find(W)) + ds.response(t);
               countDen(find(W)) = countDen(find(W)) + 1;
               mp{find(W)} = [mp{find(W)} yp(t)];

            end
            keep = find(countDen > 5);
            plot(binsC(keep), countNum(keep)./countDen(keep), 'o-', ...
               'markerfacecolor', dataClr, ...
               'markeredgecolor', 'w', ...
               'color', dataClr, ...
               'markersize', markerSize);
            plot(binsC, cellfun(@mean, mp), 'o-', ...
               'markerfacecolor', modelClr, ...
               'markeredgecolor', 'w', ...
               'color', modelClr, ...
               'markersize', markerSize);
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