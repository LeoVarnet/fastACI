# fastACI
Repository of the fastACI project

# Demos
To run the simulations from Osses & Varnet (2021, DAGA) you need to run in the MATLAB command line, and follow the instructions that will appear on the screen:

    publ_osses2021c_DAGA_1_sim;
    
To obtain figures 1 to 4 (all the paper figures) you need to run, either of the following commands:

    publ_osses2021c_DAGA_2_figs('fig1a');
    publ_osses2021c_DAGA_2_figs('fig1b');
    publ_osses2021c_DAGA_2_figs('fig2');
    publ_osses2021c_DAGA_2_figs('fig3a');
    publ_osses2021c_DAGA_2_figs('fig3b');
    publ_osses2021c_DAGA_2_figs('fig4');

# Installation
The following are the general instructions to get the fastACI toolbox for MATLAB operative in your computer. The toolbox has been tested on Windows and Linux, using MATLAB R2020b. The toolbox should be compatible with earlier MATLAB versions.

1. Download or clone the fastACI project to your local computer (one way: press the button 'Code'->Choose 'Download ZIP' and unzip somewhere).
2. This toolbox requires the Auditory Modelling Toolbox v.1.0 (AMT 1.0) that can be downloaded from [here](http://amtoolbox.org/download.php). The preferred location to store AMT 1.0 is under the MATLAB userpath (type userpath in MATLAB command's line).
3. Open and run the script **startup_fastACI.m**. This script will add all the paths under the fastACI toolbox to your local MATLAB path and it will run the script **amt_start.m** to initilise the AMT toolbox. If the AMT toolbox is not found you will be able to indicate your alternative location using a pop-up window.

# References
A. Osses Vecchi & L. Varnet (2021). **Consonant-in-noise discrimination using an auditory model with different speech-based decision devices**. DAGA conference. Vienna, Austria. ([Download paper](https://github.com/aosses-tue/fastACI/blob/main/Publications/Manuscripts/Osses-Varnet-2021-DAGA-000623.pdf))

P. Majdak, C. Hollomey, & R. Baumgartner (2021). **AMT 1.0: The toolbox for reproducible research in auditory modeling**, submitted to Acta Acustica.
