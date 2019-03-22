# Modeling-and-Visualisation-in-Physics

## Checkpoint 1 ##
### Ising model ###
Uses Monte Carlo simulation to update the states at each point if energeticaly or probabalisticaly favourable. Due to run time limitations switched over to learning how to use Cython. Finished simulation in Cython_tests\src . Run using;
```
$ python setup.py build_ext --inplace  # Build c code inplace
$ python ising_exe.py  # Run Cython code.
```
## Checkpoint 2 ##
### Game of Life ###
The game of life is implemented using parallel grid updates which decide the state of the system at the next itteration. When a glider is being used in the game of life the mean x and y positions are printed to the terminal. This can be read to a file using python gol.py > output which can then be plotted (using plotting.py) to find the rate of movment of the glider.

### SIRS Model ###
The system is stocasticaly updated changing the state of each cell from 'S' suceptable to 'I' infected if the 'S' cell is neighbouring an infected cell with probability p1. This cell can then recover 'R' at the next step with probability p2, which will then become suceptable 'S' again with probability p3. Optional argument of 'immunity' which leaves some points in a fixed state of 'R' from the begining where they can not recover to being 'S' or become 'I'.

## Checkpoint 3 ##
### Cahn Hilliard Equation ##
Separates an imulsion by movind to change the order parameter of any given point. This fillows the Cahn Hilliard equation to minimise the free energy of the system by reducing the curvature. 
