System Sounds - Bringing Planetary Systems to Life
==================================================

This code helps you make movies and MIDI files from N-body simulations using REBOUND. 

We have refactored things in a way that makes it more modular and customizable.
You now have more control over how you make your movies. 
The idea is to build up a library of examples that people can modify, rather than continually complicate a rigid framework that may not do what you want it to when you have a great idea for something new. 
Check out what's in jupyter_examples, and if you make something cool, send us a pull request so we can add it into the examples!

Installation
------------

This code uses the REBOUND N-body package (https://github.com/hannorein/rebound) to generate the animation and sound files, and will only work on Mac and Unix/Linux operating systems.

If you are unfamiliar with virtual environments, you can follow these steps to install all the dependencies without risking messing up your current python installation:

Download anaconda, python 3 version: https://www.continuum.io/downloads. The graphical installer is probably simplest and you should choose to have it append anaconda to your PATH if prompted.

Afterward in the terminal::

    $ conda create -n venv pip
    $ source activate venv
    $ conda install numpy matplotlib scipy jupyter
    $ pip install pillow rebound midiutil

This will install the libraries you need in the `venv` conda environment (you can replace this with any name), which is completely separate from your base python installation. Any time you want to use those libraries, you have to make sure you first type the source command::

    $ source activate venv

Now navigate to the directory where you would like to put a systemsounds directory, and::

    git clone https://github.com/dtamayo/systemsounds.git
    cd systemsounds/jupyter_examples/

And::

    jupyter notebook

to launch the jupyter notebook. Check MakingAMovie.ipynb for a simple example.
