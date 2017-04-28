System Sounds - Bringing Planetary Systems to Life
==================================================

Requirements
------------

This code uses the REBOUND N-body package (https://github.com/hannorein/rebound) to generate the animation and sound files, and will only work on Mac and Unix/Linux operating systems.

If you are unfamiliar with virtual environments, you can follow these steps to install all the dependencies without risking messing up your current python installation:

Download anaconda, python 3 version: https://www.continuum.io/downloads
Open terminal and navigate to your home directory::

    $ cd ~ 

Add the following line to your `.bashrc` file (if it doesn't exist, make one with that name): `export PATH=path/to/anaconda/bin$PATH`, where you should replace `/path/to/anaconda` to the path where you installed anaconda (don't forget to add bin afterward).

Back in terminal::

    $source .bashrc

Now::

    $ conda create -n venv pip
    $ source activate venv
    $ conda install numpy matplotlib scipy jupyter
    $ pip install pillow rebound midiutil

This will install most of the libraries you need in the `venv` conda environment (you can replace this with any name), which is completely separate from your base python installation. Any time you want to use it, you have to make sure you first type the source command::

    $ source activate venv

The final dependencies you need can be installed with homebrew on Mac. They have install pages you can check for unix/linux. You have to install ffmpeg with the right options to make sure you have the encoders you need::

    $ brew install timidity
    $brew install ffmpeg --with-fdk-aac --with-ffplay --with-freetype --with-libass --with-libquvi --with-libvorbis --with-libvpx --with-opus --with-x265

Now navigate to the directory where you would like to put a systemsounds directory, and::

    git clone https://github.com/dtamayo/systemsounds.git
    cd systemsounds

And::

    jupyter notebook

to launch the jupyter notebook. Check trappist.ipynb for an example.
