JAGS graphical models for a "who's more similar now" task for stimuli varying in orientation.

On each trial, two pairs of stimuli are presented. The response is a binary decision about which of the two pairs has more similar stimuli.

# Version 3

Version 3 assumes all the orientation stimuli lie on a bounded number line.

`thurstoneFromChoice_3.m` allows data sets (`whoIsSimilar_A`, `whoIsSimilar_B`, ...) to be generated with various task properties
-number of stimuli
-true orientations
-number of trials in the experiment

Inferences about the generated data are then made using `thurstoneFromChoice_3x_jags.txt` graphical models, where x allows for model variants.
