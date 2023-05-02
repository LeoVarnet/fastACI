# fastACI toolbox
This is the repository of the fast Auditory Classification Images (fastACI) project. The toolbox is controlled using the command line of MATLAB. It does not have (yet) a graphical interface.

With this toolbox you can run listening experiments as used in the studies **varnet2013**, **varnet2015**, **varnet2021**, **osses2021c**, **osses2022b**,  **varnet2022a**, **osses2023a**, **osses2022b**, and **carranante2023** (see the full citations in the section "_References_"). You can also reproduce some of the figures contained in the mentioned references.

| Citation key   | fastACI experiment name     | Type of background noise     | Target sounds |
| :------------- | :----------: | :-----------: | :-----------: |
| **varnet2013** | `speechACI_varnet2013`   | white  | /aba/-/ada/, female speaker |
| **varnet2015** | `speechACI_varnet2015`   | white  | /alda/-/alga/-/arda/-/arga/, male speaker |
| **osses2021c** | `speechACI_varnet2013`   | speech shaped noise (SSN) |  /aba/-/ada/, female speaker |
| **varnet2022a**| `modulationACI`          | white  | modulated or unmodulated tones |
| **osses2022b** | `speechACI_Logatome`     | white, bump, MPS | /aba/-/ada/, male speaker (S43M) from the OLLO database |
| **carranante2023** | `speechACI_Logatome` | bump | Pairs of contrasts using /aba/, /ada/, /aga/, /apa/, /ata/ from the same male speaker (S43M) in OLLO |
| **osses2023a** | `segmentation`           | random prosody | Pairs of contrasts: /l'amie/-/la mie/, /l'appel/-/la pelle/, /l'accroche/-/la croche/, /l'alarme/-/la larme/ |
| **osses2023b** | `toneinnoise_ahumada1975` | white | Tone-in-noise experiment with 100-ms 500-Hz sinusoids temporally centred in Gaussian noises of 500 ms|

Make sure that you follow the steps indicated in the section **Installation** (below) the first time you use the toolbox. 

# How to cite this repository
This repository can be cited as follows: The fastACI toolbox was used (Osses & Varnet, 2022).

**If a model version is cited (in this example: release fastACI v1.2):**

A. Osses & L. Varnet (2022). "fastACI toolbox: the MATLAB toolbox for investigating auditory perception using reverse correlation (v1.2)" [![DOI](https://zenodo.org/badge/335310799.svg)](https://zenodo.org/badge/latestdoi/335310799)

**If a specific commit is cited (in this example: commit cc9d9cf):**

A. Osses & L. Varnet (2022). "fastACI toolbox: the MATLAB toolbox for investigating auditory perception using reverse correlation," Github commit cc9d9cf.

<a rel="license" href="http://creativecommons.org/licenses/by-nc/4.0/"><img alt="Creative Commons Licence" style="border-width:0" src="https://i.creativecommons.org/l/by-nc/4.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by-nc/4.0/">Creative Commons Attribution-NonCommercial 4.0 International License</a>.

# Running a listening experiment
Next we present the command line required to run each of the ACI experiments that are available in our toolbox. These examples assume that the listener will be named `S01` (standing for Subject 01), however any character-based name can be used instead.

### modulationACI: Experiment as in varnet2021 
    fastACI_experiment('modulationACI','S01');
    
### speechACI_varnet2013
    fastACI_experiment('speechACI_varnet2013','S01','white'); % to run it as in varnet2013
    fastACI_experiment('speechACI_varnet2013','S01','SSN');   % to run it as in osses2021c
    
### speechACI_varnet2015: Experiment as in varnet2015    
Alda/Alga/Arda/Arga discrimination using a male speaker
    
    fastACI_experiment('speechACI_varnet2015','S01','white');
    
### speechACI_Logatome: Experiments using speech samples from the Logatome corpus in French (under construction)
Aba/Ada discrimination using a female speaker (S41F from the Logatome corpus):

    fastACI_experiment('speechACI_Logatome-abda-S41F','S01','white');

Aba/Ada discrimination using a male speaker (S43M from the Logatome corpus), as in **osses2022b**:

    fastACI_experiment('speechACI_Logatome-abda-S43M','S01','white');
    
Additional VCV contrasts using a male speaker (S43M from the Logatome corpus), as in **carranante2023**, using bump noises:
    
    fastACI_experiment('speechACI_Logatome-abda-S43M','S01','bumpv1p2_10dB');
    fastACI_experiment('speechACI_Logatome-adga-S43M','S01','bumpv1p2_10dB');
    fastACI_experiment('speechACI_Logatome-apta-S43M','S01','bumpv1p2_10dB');    
    fastACI_experiment('speechACI_Logatome-abpa-S43M','S01','bumpv1p2_10dB');
    fastACI_experiment('speechACI_Logatome-adta-S43M','S01','bumpv1p2_10dB');    
    
# Simulating a listening experiment
A listening experiment can be simulated using an artificial listener or, in other words, an auditory model. So far, we have validated the use of the models `osses2021` (Osses & Kohlrausch, 2021), `osses2022a` (to be published), and `king2019` (King et al., 2019). The models `osses2021` and `king2019` are both available within AMT 1.0 (or more recent), `osses2022a` is exclusively available in our toolbox.

To run simulations you only have to use the corresponding model as the subject name. To use `osses2021` in the simulation of the experiment `speechACI_varnet2013` using SSN noises, you need to type in MATLAB:

    fastACI_experiment('speechACI_varnet2013','osses2021','SSN');   % to run it as in osses2021c
       
or, to use `king2019`:       

    fastACI_experiment('speechACI_varnet2013','king2019','SSN');   
    
More elaborate simulations can be automatically run using the scripts `pres_osses2022_02_AABBA_1_sim.m` and `publ_osses2021c_DAGA_1_sim.m`, among other scripts. In the next section, all the simulations in **osses2021c** can be reproduced using the `osses2021` model with two different decision back ends. This is related to the script `publ_osses2021c_DAGA_1_sim.m`.

# Demo: Obtaining the figures from osses2021c
To run the simulations from Osses & Varnet (2021, DAGA) you need to run in the MATLAB command line, and follow the instructions that will appear on the screen:

    publ_osses2021c_DAGA_1_sim;
    
To obtain figures 1 to 4 (all the paper figures) you need to run, either of the following commands:

    publ_osses2021c_DAGA_2_figs('fig1a'); % REQUIRED: manual download of experimental data (see below)
    publ_osses2021c_DAGA_2_figs('fig1b'); % REQUIRED: manual download of experimental data (see below)
    publ_osses2021c_DAGA_2_figs('fig2');
    publ_osses2021c_DAGA_2_figs('fig3a');
    publ_osses2021c_DAGA_2_figs('fig3b');
    publ_osses2021c_DAGA_2_figs('fig4'); % REQUIRED: manual download of experimental data (see below)
    
To obtain Fig 1A or Fig 1B, you require to manually download (in advance) the experimental dataset, which is available on Zenodo (see ref. **osses2021c_data**).

    publ_osses2021c_DAGA_0_checkdata;

# Installation
The following are the general instructions to get the fastACI toolbox for MATLAB operative in your computer. The toolbox has been tested on Windows and Linux, using MATLAB (versions R2012b-R2020b).

1. Download or clone the fastACI project to your local computer (one way: press the button 'Code'->Choose 'Download ZIP' and unzip somewhere).
2. This toolbox requires the Auditory Modelling Toolbox v.1.0 (AMT 1.0 or higher) that can be downloaded from [here](http://amtoolbox.org/download.php). After the download you are not expected to do anything else, as the AMT toolbox will automatically be initialised in our next step:
3. Open and run the script **startup_fastACI.m**. This script will add all the paths under the fastACI toolbox to your local MATLAB path and it will run the script **amt_start.m** to initilise the AMT toolbox. If the AMT toolbox is not found you will be able to indicate your alternative location using a pop-up window.

# References for the fastACI toolbox
The references are sorted alphabetically (first author's last name) and then more recent first.
|    |  |
| :------------- | :---------- | 
| **carranante2023** | G. Carranante, M. Giavazzi, & L. Varnet (2023). **Auditory reverse correlation applied to the study of place and voicing: Four new phoneme-discrimination tasks**. To be presented at Forum Acusticum 2023. |
| **king2019**   | A. King, L. Varnet, & C. Lorenzi (2019). **Accounting for masking of frequency modulation by amplitude modulation with the modulation filter-bank concept**. J. Acoust. Soc. Am. 145, p. 2277-2293 (doi: [10.1121/1.5094344](http://dx.doi.org/10.1121/1.5094344), [Download paper](https://hal.archives-ouvertes.fr/hal-02993025))|
| **osses2023a** | A. Osses, E. Spinelli, F. Meunier, E. Gaudrain, & L. Varnet (2023). **Prosodic cues to word boundaries in a segmentation task using reverse correlation**. To be submitted. |
| **osses2023a_data** | A. Osses, E. Spinelli, F. Meunier, E. Gaudrain, & L. Varnet (2023). **Raw and post-processed data for the study of prosodic cues to word boundaries in a segmentation task using reverse correlation** (doi: [10.5281/zenodo.7865424](https://doi.org/10.5281/zenodo.7865424))| 
| **osses2023b** | A. Osses, & L. Varnet (2023). **Using auditory models to mimic human listeners in reverse correlation experiments from the fastACI toolbox**. To be presented at Forum Acusticum 2023. |
| **osses2023b_data** | A. Osses, & L. Varnet (2023). **Raw and post-processing data for using auditory models to mimic human listeners in reverse correlation experiments from the fastACI toolbox** (doi: [10.5281/zenodo.7886232](https://doi.org/10.5281/zenodo.7886232))| 
| **osses2022d** | A. Osses, C. Lorenzi, & L. Varnet (2022). **Assessment of individual listening strategies in amplitude-modulation detection and phoneme categorisation tasks**. International Congress on Acoustics, 24-28 October, Gyeongju, Korea ([Download presentation](https://github.com/aosses-tue/fastACI/blob/main/Publications/Manuscripts/Osses2022d-Lorenzi-Varnet-ICA_presentation.pdf), [Download proceedings](https://ica2022korea.org/data/Proceedings_A11.pdf)) |
| **osses2021c** | A. Osses & L. Varnet (2021). **Consonant-in-noise discrimination using an auditory model with different speech-based decision devices**. DAGA conference. Vienna, Austria. ([Download paper](https://github.com/aosses-tue/fastACI/blob/main/Publications/Manuscripts/Osses-Varnet-2021-DAGA-000623.pdf))|
| **osses2021c_data** | A. Osses & L. Varnet (2021). **Noise data for the study of consonant-in-noise discrimination using an auditory model with different speech-based decision devices**. Experimental data for **osses2021c** (doi: [10.5281/zenodo.5483835](https://doi.org/10.5281/zenodo.5483835)) |
| **varnet2022a** | L. Varnet & C. Lorenzi (2022). **Probing temporal modulation detection in white noise using intrinsic envelope fluctuations: A reverse correlation study**. J. Acoust. Soc. Am. 151, p. 1356-1366 (doi: [10.1121/10.0009629](https://doi.org/10.1121/10.0009629))|
| **varnet2022a_data** | L. Varnet (2021). **AM revcorr data**. Experimental data for **varnet2022a** (doi: [10.5281/zenodo.5571719](https://doi.org/10.5281/zenodo.5571719)) |
| **varnet2022b**| L. Varnet, C. Lorenzi, & A. Osses (2022). **Probing amplitude-modulation detection and phoneme categorization with auditory reverse correlation**. Congrès Français d'Acoustique, 11-15 April, Marseille, France ([Download presentation](https://github.com/aosses-tue/fastACI/blob/main/Publications/Manuscripts/Varnet2022b-CFA.pdf))|
| **varnet2015** | L. Varnet, K. Knoblauch, W. Serniclaes, F. Meunier, & M. Hoen (2015). **A psychophysical imaging method evidencing auditory cue extraction during speech perception: A group analysis of auditory classification images**. PLoS one 3, p. 1-23 ([Download paper](https://hal.archives-ouvertes.fr/hal-01132995))|
| **varnet2013** | L. Varnet, K. Knoblauch, F. Meunier, & M. Hoen (2013). **Using auditory classification images for the identification of fine acoustic cues used in speech perception**. Front. Hum. Neurosci. 7, p. 1-12 ([Download paper](https://hal.archives-ouvertes.fr/hal-00931465))|

# Other references
P. Majdak, C. Hollomey, & R. Baumgartner (2022). **AMT 1.x: A toolbox for reproducible research in auditory modeling**, Acta Acustica, 6, 19. (doi: [10.1051/aacus/2022011](https://doi.org/10.1051/aacus/2022011)).

A. Osses, L. Varnet, L. Carney, T. Dau, I. Bruce, S. Verhulst, & P. Majdak (2022). **A comparative study of eight human auditory models of monaural processing**, Acta Acustica, 6, 17 (doi: [10.1051/aacus/2022008](https://doi.org/10.1051/aacus/2022008)).

A. Osses & A. Kohlrausch (2021). **Perceptual similarity between piano notes: Simulations with a template-based perception model**. J. Acoust. Soc. Am. 149, p. 3534-3552 (doi: [10.1121/10.0004818](https://asa.scitation.org/doi/abs/10.1121/10.0004818)).

# Acknowledgements
The development of the fastACI toolbox was funded by the ANR grant ["fastACI"](https://anr.fr/Project-ANR-20-CE28-0004) attributed to Léo Varnet (ANR-20-CE28-0004) and was further supported by the ["FrontCog"](https://anr.fr/ProjetIA-17-EURE-0017) grant (ANR-17-EURE-0017).

<a rel="license" href="http://creativecommons.org/licenses/by-nc/4.0/"><img alt="Creative Commons Licence" style="border-width:0" src="https://i.creativecommons.org/l/by-nc/4.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by-nc/4.0/">Creative Commons Attribution-NonCommercial 4.0 International License</a>.
