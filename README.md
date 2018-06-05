# AntColonySystem

This is a (relatively) simple, single-file implementation of the Ant Colony
System algorithm by Dorigo & Gambardella as described in: *Dorigo, Marco, and
Luca Maria Gambardella. "Ant colony system: a cooperative learning approach to
the traveling salesman problem." IEEE Transactions on evolutionary computation
1.1 (1997): 53-66.*

It is meant for educational purposes and lacks generality and uses by default
so called _candidate lists_ to speed up the search process.

It supports (partially) loading of (A)TSP instances from the well-known [TSPLIB] repository.

It is written in C++ 11.

## Building

Run `make` from the command line (it requires GCC version with C++ 11 support).
The default target is `release`, i.e. optimized version but it can be changed to `debug`.

## Running

By default it tries to load `kroA100.tsp` TSP instance from the current directory
but the instance path can be given as a parameter.

As the algorihtm is being executed, the program outputs information about 
each improved solution found.

Example output:

    Read line: NAME: kroA100
    Read line: TYPE: TSP
    Read line: COMMENT: 100-city problem A (Krolak/Felts/Nelson)
    Read line: DIMENSION: 100
    Read line: EDGE_WEIGHT_TYPE : EUC_2D
    Read line: NODE_COORD_SECTION
    Greedy solution cost: 27807
    initial_pheromone: 3.59622e-07
    New best solution found with cost: 27696 at iteration 0
    New best solution found with cost: 26399 at iteration 0
    New best solution found with cost: 25311 at iteration 0
    New best solution found with cost: 25126 at iteration 1
    New best solution found with cost: 25065 at iteration 1
    New best solution found with cost: 24807 at iteration 2
    New best solution found with cost: 24694 at iteration 6
    ...
    New best solution found with cost: 21368 at iteration 58240
    New best solution found with cost: 21282 at iteration 58241

## Alternative implementation in C#

For those who prefer C#:
https://github.com/mbalchanowski/Ant-Colony-System

[TSPLIB]: <http://comopt.ifi.uni-heidelberg.de/software/TSPLIB95/>
