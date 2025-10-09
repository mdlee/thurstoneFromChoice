%% Similarity compatison model

clear; close all;
preLoad = true;

% graphical model script
modelDir = './';
modelName = 'similarityComparison';
engine = 'jags';

% data sets
dataList = {...
   'tomicBaysSimilarity'; ...
   };

%% general constants
pi = 3.1415;

% loop over data
for dataIdx = 1:numel(dataList)
   dataName = dataList{dataIdx};
   switch dataName

      case 'tomicBaysSimilarity'
         dataDir = 'data/';
         dataName = 'tomicBays';
         load([dataDir dataName], 'ds');

         a = ds.aIdx;
         b = ds.bIdx;
         c = ds.cIdx;
         d = ds.dIdx;
         y = ds.response;
         nTrials = ds.nTrials;
         nStimuli = ds.nStimuli;
   end

   %% sampling from graphical model
   % parameters to monitor
   params = {'mu', 'sigma'};

   % MCMC properties
   nChains    = 8;     % number of MCMC chains
   nBurnin    = 2e3;   % number of discarded burn-in samples
   nSamples   = 2e3;   % number of collected samples
   nThin      = 10;    % number of samples between those collected
   doParallel = 1;     % whether MATLAB parallel toolbox parallizes chains

   % assign MATLAB variables to the observed nodes
   data = struct(...
      'a'        , a        , ...
      'b'        , b        , ...
      'c'        , c        , ...
      'd'        , d        , ...
      'y'        , y        , ...
      'nStimuli' , nStimuli , ...
      'nTrials'  , nTrials  );

   % censoring initial values so data have likelihood on first sample
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
   % (note intialization of mu, which encourages the
   % better log-likelihood representation)
   generator = @()struct(...
      'xA', xAinit, ...
      'xB', xBinit, ...
      'xC', xCinit, ...
      'xD', xDinit, ...
      'mu', ds.stimuli);

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

      % just convergent enough chains
      [keepChains, rHat] = findKeepChains(chains.sigma, 2, 1.1);
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
      save(sprintf('storage/%s', fileName), 'chains', 'stats', 'diagnostics', 'info', '-v7.3');

   end


   % log likelihood
   if contains(params, 'yp')
      LL = sum(y.*log(yp)) + sum((1-y).*log(1-yp));
      fprintf('Log-likelihood = %1.4f\n', LL);
   end

end