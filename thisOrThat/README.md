JAGS graphical models for a "this or that" task for stimuli varying in orientation.

On each trial, a target is specified, and two stimuli are presented: one is the target and the other is a foil. The response is a binary decision about which of the two stimuli is the target.

`thurstoneFromChoice_1.m` allows data sets (`thisOrThat_A`, `thisOrThat_B`, ...) to be generated with various task properties
- number of stimuli
- true orientations
- number of trials in the experiment

Inferences about the generated data are then made using `thurstoneFromChoice_1x_jags.txt` graphical models, where x allows for model variants.
