JAGS graphical models for a "this or that" task for stimuli varying in orientation.

On each trial, a target is specified, and two stimuli are presented: one is the target and the other is a foil. The response is a binary decision about which of the two stimuli is the target.

# Version 1

Version 1 assumes all the orientation stimuli lie on a bounded number line.

`thurstoneFromChoice_1.m` allows data sets (`thisOrThat_A`, `thisOrThat_B`, ...) to be generated with various task properties
- number of stimuli
- true orientations
- number of trials in the experiment

Inferences about the generated data are then made using `thurstoneFromChoice_1x_jags.txt` graphical models, where x allows for model variants.

# Version 2

Version 2 assumes all the orientation stimuli lie on a circle.

`thurstoneFromChoice_2.m` allows data sets (`thisOrCircularThat_A`, `thisOrCircularThat_B`, ...) to be generated with various task properties
- number of stimuli
- true orientations
- number of trials in the experiment

Inferences about the generated data are then made using `thurstoneFromChoice_2x_jags.txt` graphical models, where x allows for model variants.
