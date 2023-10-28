# CICASS
initial conditions generator that self-consistently includes streaming velocity of the baryons 
relative to the dark matter, for cosmological simulations

This code is from https://arxiv.org/abs/1204.1344 and https://arxiv.org/abs/1204.1345 and has been used by many
subsequent studies of streaming velocity.  Its chief use is to initialize streaming velocity simulations
well after z=1000 by correctly evolving the linear equations and then initializing.  An alternative 
is to start at z=1000 as in https://academic.oup.com/mnras/article/474/2/2173/4622957, and then because
the baryons are very uniform you can boost the velocity of the baryons without the need for CICASS.
The code was written a long time ago and is likely too inefficient to generate ICs for simulations 
larger than 1024^3.

Simulations that were run with the CICASS ICs are shown here: https://sites.google.com/uw.edu/mattmcquinn/movies


The subdirectory makeCosICs/ contains the primary README file and the
IC generating code.  

The subdirectory vbc_transfer/ contains the code
that outputs the transfer function in the format that is required by
the makeCosICs code.  In the case with a velocity difference between
the baryons and dark matter, it generates a transfer function that has
a net relative velocity.

The subdirectory glass/ stores a glass file that has been evolved for
10^4 in scale factor as suggested in Crocce et al (2006) to minimize transients.
