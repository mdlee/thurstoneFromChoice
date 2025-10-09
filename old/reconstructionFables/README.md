JAGS graphical models for a perceptual and memory reconstructions tasks for stimuli varying in orientation.

In the perceptual task, on each trial, a stimulus is presented, and the participant uses an interface to reconstruct it.

In the memory task, 3 or 6 stimuli are presented then removed, and then one of them is revealed to be the targer. The participant uses an interface to reconstruct the target.

# Version 7

Version 7a is for the perceptual reconstruction task.

`thurstoneFromChoice_7a_jags.txt` is the JAGS graphical model, and `thurstoneFromChoice_7a_stan.txt` is a Stan implementation

`thurstoneFromChoice_7a.m` allows infererences about the `tomicBays` empirical data set from [here](https://psycnet.apa.org/record/2023-21056-001). Via the trinity `callBayes` function it can use either the JAGS or Stan implementation.

`thurstoneFromChoice_7a.m` and `thurstoneFromChoice_7b_jags.txt` are a non-pursued attempt at a joint (common cause) model of perceptual and memory reconstruction.

# Version 8

Version is for the memory reconstruction task.

Version 8a to 8d are non-pursued attempts at weighted sum models that assume the reconstructed stimulus is a convex combination of the target and foils.

Version 8e implements the swap error model, for just one setsize (2 foils or 5).

`thurstoneFromChoice_8f_jags.txt` is the JAGS graphical model for the swap error model for both set sizes simultaneously. `thurstoneFromChoice_8f_jags.txt`, with the order constraint that reconstruction noise is greater for the 5 foil trials. `thurstoneFromChoice_8f2_jags.txt` omits that order constraint.

`thurstoneFromChoice_8f.m` and `thurstoneFromChoice_8f2.m` allows infererences about the `tomicBays` empirical data set from [here](https://psycnet.apa.org/record/2023-21056-001).

`thurstoneFromChoice_8g_jags.txt` is `thurstoneFromChoice_8f_jags.txt` without the swap error process. It is applied to data by `thurstoneFromChoice_8g.m` 