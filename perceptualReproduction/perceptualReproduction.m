%% Perceptual reproduction model

clear; close all;
preLoad = true;

% graphical model script
modelDir = './';
modelName = 'perceptualReproduction';
engine = 'jags';
% engine = 'stan';

% data sets
dataList = {...
   'tomicBaysPerception'; ...
   };

%% general constants
pi = 3.1415;

% loop over data
for dataIdx = 1:numel(dataList)
   dataName = dataList{dataIdx};
   switch dataName

      case 'tomicBaysPerception'
         dataDir = 'data/';
         dataName = 'tomicBays';
         load([dataDir dataName], 'dp');

         y = dp.response;
         s = dp.sIdx;
         nTrials = dp.nTrials;
         nStimuli = dp.nStimuli;
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
      's'  , s      , ...
      'y'   , y   , ...
      'nStimuli', nStimuli, ...
      'nTrials', nTrials);

   % generator for initialization
   generator = @()struct('sigma', rand*pi);

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
         'allowunderscores', 1                                         , ...
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
end
