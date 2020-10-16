```
              _                 _      
  /\/\   ___ | | ___  ___ _   _| | ___ 
 /    \ / _ \| |/ _ \/ __| | | | |/ _ \
/ /\/\ \ (_) | |  __/ (__| |_| | |  __/
\/    \/\___/|_|\___|\___|\__,_|_|\___|
   __           _                                      _   
  /__\ ____   _(_)_ __ ___  _ __  _ __ ___   ___ _ __ | |_ 
 /_\| '_ \ \ / / | '__/ _ \| '_ \| '_ ` _ \ / _ \ '_ \| __|
//__| | | \ V /| | | | (_) | | | | | | | | |  __/ | | | |_ 
\__/|_| |_|\_/ |_|_|  \___/|_| |_|_| |_| |_|\___|_| |_|\__|

```                                                        

Reinforcement learning environment for inverse drug design.

## Set-up

Install dependencies:
```
conda env create -f environment.yml
```
Activate environment:
```
conda activate mol-env
```
Other Libraries That needs to be installed
conda install kivy -c conda-forge
easy_install networkx


The reward function inside the environment has a default reward of 10 for demonstration purposes. This can be modified as you wish. 

An agent class with skeleton code has been provided for you to modify to your needs.

To run the program please this in terminal

Python gui.py