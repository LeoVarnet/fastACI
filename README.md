# fastACI toolbox
This is the repository of the fast Auditory Classification Images (fastACI) project. The toolbox is controlled using the command line of MATLAB. It does not have (yet) a graphical interface.

With this toolbox you can run listening experiments as used in the studies by Varnet et al. (2013, 2015), Varnet & Lorenzi (2021), and Osses & Varnet (2021). You can also reproduce some of the figures contained in the mentioned references.

| Citation key   | fastACI experiment name     | Type of background noise     |
| :------------- | :----------: | :-----------: |
| **varnet2013** | `speechACI_varnet2013`   | white  |
| **varnet2015** | `speechACI_varnet2015`   | white  |
| **osses2021c** | `speechACI_varnet2013`   | speech shaped noise (SSN) |
| **varnet2021** | `modulationACI`          | white  |

# How to cite this repository
This repository can be cited as follows: The fastACI toolbox was used (Osses & Varnet, 2021).

**If a specific commit is cited (in this example: commit b18d919):**

A. Osses Vecchi & L. Varnet (2021). "fastACI toolbox: Exploring phoneme representations and their applicability using Auditory Classification Images," Github commit b18d919.

# Running a listening experiment
Next we present the command line required to run each of the ACI experiments that are available in our toolbox. These examples assume that the listener will be named `S01` (standing for Subject 01), however any character-based name can be used instead.

### modulationACI: Experiment as in varnet2021 
    fastACI_experiment('modulationACI','S01');
    
### modulationACI_seeds: Experiment as in varnet2021 (alternative without storing files on disk)
    fastACI_experiment('modulationACI_seeds','S01');
    
### speechACI_varnet2013
    fastACI_experiment('speechACI_varnet2013','S01','white'); % to run it as in varnet2013
    fastACI_experiment('speechACI_varnet2013','S01','SSN');   % to run it as in osses2021c
    
### speechACI_varnet2015: Experiment as in varnet2015    
Alda/Alga/Arda/Arga discrimination using a male speaker
    
    fastACI_experiment('speechACI_varnet2015','S01','white');
    
### speechACI_Logatome: Experiments using speech samples from the Logatome corpus in French (under construction)
Aba/Ada discrimination using a female speaker (S41F from the Logatome corpus):

    fastACI_experiment('speechACI_Logatome-abda-S41F','S01','white');

# Simulating a listening experiment
A listening experiment can be simulated using an artificial listener or, in other words, an auditory model. So far, we have validated the use of the models `osses2021` (Osses & Kohlrausch, 2021) and `king2019` (King et al., 2019), both available within AMT 1.0.

To run simulations you only have to use the corresponding model as the subject name. To use `osses2021` in the simulation of the experiment `speechACI_varnet2013` using SSN noises, you need to type in MATLAB:

    fastACI_experiment('speechACI_varnet2013','osses2021','SSN');   % to run it as in osses2021c
       
or, to use `king2019`:       

    fastACI_experiment('speechACI_varnet2013','king2019','SSN');   
    
See also the next session, where all the simulations in **osses2021c** can be reproduced using the `osses2021` model with two different decision back ends.

# Demo: Obtaining the figures from osses2021c
To run the simulations from Osses & Varnet (2021, DAGA) you need to run in the MATLAB command line, and follow the instructions that will appear on the screen:

    publ_osses2021c_DAGA_1_sim;
    
To obtain figures 1 to 4 (all the paper figures) you need to run, either of the following commands:

    publ_osses2021c_DAGA_2_figs('fig1a'); % REQUIRED: manual download of experimental data (see below)
    publ_osses2021c_DAGA_2_figs('fig1b'); % REQUIRED: manual download of experimental data (see below)
    publ_osses2021c_DAGA_2_figs('fig2');
    publ_osses2021c_DAGA_2_figs('fig3a');
    publ_osses2021c_DAGA_2_figs('fig3b');
    publ_osses2021c_DAGA_2_figs('fig4');
    
To obtain Fig 1A or Fig 1B, you require to manually download (in advance) the experimental dataset, which is available on Zenodo (see ref. **osses2021c_data**).

    publ_osses2021c_DAGA_0_checkdata;

# Installation
The following are the general instructions to get the fastACI toolbox for MATLAB operative in your computer. The toolbox has been tested on Windows and Linux, using MATLAB R2020b. The toolbox should be compatible with earlier MATLAB versions.

1. Download or clone the fastACI project to your local computer (one way: press the button 'Code'->Choose 'Download ZIP' and unzip somewhere).
2. This toolbox requires the Auditory Modelling Toolbox v.1.0 (AMT 1.0) that can be downloaded from [here](http://amtoolbox.org/download.php). The preferred location to store AMT 1.0 is under the MATLAB userpath (type userpath in MATLAB command's line).
3. Open and run the script **startup_fastACI.m**. This script will add all the paths under the fastACI toolbox to your local MATLAB path and it will run the script **amt_start.m** to initilise the AMT toolbox. If the AMT toolbox is not found you will be able to indicate your alternative location using a pop-up window.

# References for the fastACI toolbox
|    |  |
| :------------- | :---------- | 
| **king2019**   | A. King, L. Varnet, & C. Lorenzi (2019). **Accounting for masking of frequency modulation by amplitude modulation with the modulation filter-bank concept**. J. Acoust. Soc. Am. 145, p. 2277-2293 (Doi: [10.1121/1.5094344](http://dx.doi.org/10.1121/1.5094344), [Download paper](https://hal.archives-ouvertes.fr/hal-02993025))|
| **osses2021c** | A. Osses Vecchi & L. Varnet (2021). **Consonant-in-noise discrimination using an auditory model with different speech-based decision devices**. DAGA conference. Vienna, Austria. ([Download paper](https://github.com/aosses-tue/fastACI/blob/main/Publications/Manuscripts/Osses-Varnet-2021-DAGA-000623.pdf))|
| **osses2021c_data** | A. Osses Vecchi & L. Varnet (2021). **Noise data for the study of consonant-in-noise discrimination using an auditory model with different speech-based decision devices**. Experimental data for **osses2021c** (Doi: [10.5281/zenodo.5483835](https://doi.org/10.5281/zenodo.5483835)) |
| **osses2021a** | A. Osses Vecchi & A. Kohlrausch (2021). **Perceptual similarity between piano notes: Simulations with a template-based perception model**. J. Acoust. Soc. Am. 149, p. 3534-3552 (Doi: [10.1121/10.0004818](https://asa.scitation.org/doi/abs/10.1121/10.0004818))|
| **varnet2021** | L. Varnet & C. Lorenzi (2021). **Probing temporal modulation detection in white noise using intrinsic envelope fluctuations: A reverse correlation study**. Submitted to J. Acoust. Soc. Am. |
| **varnet2015** | L. Varnet, K. Knoblauch, W. Serniclaes, F. Meunier, & M. Hoen (2015). **A psychophysical imaging method evidencing auditory cue extraction during speech perception: A group analysis of auditory classification images**. PLoS one 3, p. 1-23 ([Download paper](https://hal.archives-ouvertes.fr/hal-01132995))|
| **varnet2013** | L. Varnet, K. Knoblauch, F. Meunier, & M. Hoen (2013). **Using auditory classification images for the identification of fine acoustic cues used in speech perception**. Front. Hum. Neurosci. 7, p. 1-12 ([Download paper](https://hal.archives-ouvertes.fr/hal-00931465))|

# Other references
P. Majdak, C. Hollomey, & R. Baumgartner (2021). **AMT 1.0: The toolbox for reproducible research in auditory modeling**, submitted to Acta Acustica.

# Acknowledgements
The development of the fast-ACI toolbox was funded by the ANR grant ["fastACI"](https://anr.fr/Project-ANR-20-CE28-0004) attributed to Léo Varnet (ANR-20-CE28-0004).
