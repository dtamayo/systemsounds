System Sounds - Bringing Planetary Systems to Life
==================================================

Requirements
------------

This code uses the REBOUND N-body package (https://github.com/hannorein/rebound) to generate the animation and sound files, and will only work on Mac and Unix/Linux operating systems.

If you are unfamiliar with virtual environments, you can follow these steps to install all the dependencies without risking messing up your current python installation:

Download anaconda, python 3 version: https://www.continuum.io/downloads
Open terminal and navigate to your home directory::
    $ Cd ~

Add the following line to your `.bashrc` file (if it doesn't exist, make one with that name)
    Export PATH=path/to/anaconda/bin$PATH

Back in terminal::

    source .bashrc

conda create -n venv pip
source activate venv
conda install numpy matplotlib scipy jupyter
pip install pillow rebound midiutil

brew install timidity
brew install ffmpeg --with-fdk-aac --with-ffplay --with-freetype --with-libass --with-libquvi --with-libvorbis --with-libvpx --with-opus --with-x265

git clone https://github.com/dtamayo/systemsounds.git
Cd systemsounds
jupyter notebook

Git pull to get updates
