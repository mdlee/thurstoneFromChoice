JAGS graphical models for a "who's more similar now" task for stimuli varying in orientation.

On each trial, two pairs of stimuli are presented. The response is a binary decision about which of the two pairs has more similar stimuli.

# Version 3

Version 3 assumes all the orientation stimuli lie on a bounded number line.

`thurstoneFromChoice_3.m` allows data sets (`whoIsSimilar_A`, `whoIsSimilar_B`, ...) to be generated with various task properties
-number of stimuli
-true orientations
-number of trials in the experiment

Inferences about the generated data are then made using `thurstoneFromChoice_3x_jags.txt` graphical models, where x allows for model variants.

# Version 4

Version 4 assumes all the orientation stimuli lie on a circle, and that the circle wraps around.

`thurstoneFromChoice_4.m` allows data sets (`whoIsSimilar_A`, `whoIsSimilar_B`, ...) to be generated with various task properties
-number of stimuli
-true orientations
-number of trials in the experiment

Inferences about the generated data are then made using `thurstoneFromChoice_4x_jags.txt` graphical models, where x allows for model variants
- whether the orientations are represented modulo pi (so angle x and x+pi are regarded as same orientation)
- whether the distances are taken modulo pi (to insure the smaller angle between two orientations is the distance)

# Version 5

Version 5 assumes all the orientation stimuli lie on a circle, but the circle does not wrap around

`thurstoneFromChoice_5.m` allows data sets (`whoIsSimilar_A`, `whoIsSimilar_B`, ...) to be generated with various task properties
-number of stimuli
-true orientations
-number of trials in the experiment

`thurstoneFromChoice_4.m` also allows infererences about the `tomicBays` empirical data set from [here](https://psycnet.apa.org/record/2023-21056-001)

Inferences about the generated data are then made using `thurstoneFromChoice_5x_jags.txt` graphical models, where x allows for model variants
- `5b` introduces identifiability constraints
