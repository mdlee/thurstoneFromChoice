# Memory reproduction model

On each trial 3 or 6 oriented line stimuli are presented then removed. One is then probed as the target, and the participant uses an interface to reproduce it.

`memoryReproduction_jags.txt` implements the full memory reproduction comparison model as a graphical model in JAGS.

`memoryReproductionNoSwap_jags.txt` implements a reduced version of the model that removes the swap error process.

`memoryReproduction.m` and `memoryReproductionNoSwap.m` are MATLAB scripts that appy the models to data. The `tomicBays` empirical data set from [here](https://psycnet.apa.org/record/2023-21056-001) is implemented, but others could be added.
