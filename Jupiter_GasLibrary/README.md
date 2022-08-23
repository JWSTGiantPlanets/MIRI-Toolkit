## Jupiter Spectral Forward Models

This directory contains spectral forward models calculated by the [NEMESIS software](https://github.com/nemesiscode/radtrancode).  These are calculated for the twelve subbands of MIRI/MRS, using the spectral resolution defined by Labiano et al. (2020).

Each gas has two PNG files, SHORT (channels 1A to 2C) and LONG (channels 3A to 4C), although the LONG channels are expected to fully saturate on Jupiter.

Each file contains multiple panels - the forward model for 0.1x, 0.5x, 1.0x, 2.0x and 5.0x the nominal gaseous abundance; then the residual between the 0.1x abundance case and everything else.  These charts provide a quick-look to understand which gases contribute to each channel.
