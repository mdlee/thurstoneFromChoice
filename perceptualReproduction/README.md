# Perceptual reproduction model

On each trial an oriented line stimulus is presented and the participant uses an interface to reproduce it.

`perceptualReproduction_jags.txt` and `perceptualReproduction_stan.txt` implement the perceptual reproduction comparison model as a graphical model in JAGS and Stan.

`perceptualReproduction.m` is a MATLAB script that applies the model to data. The `tomicBays` empirical data set from [here](https://psycnet.apa.org/record/2023-21056-001) is implemented, but others could be added.
