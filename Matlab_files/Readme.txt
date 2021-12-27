Optimization workflow

Functions 

1) main 
2) update control schedules 
3) trajectory,py
4) NPV
5)outrun
6)forward simulation
7)forward simulation and gradient
8) Final rates

Main file

1) Standard setup

1. Clear all variables, close everything
2. start timing
3. Set up directories


2)  Set reservoir dimensions

1. Number of grids, dimension of reservoir


3) Set up grids, load res properties and assign to rock
1. #of grids, geometry
2. Properties

4) Set up controls
1. # of time steps, time, volume
2. Control limits, rate limits


5) Set up wells
1. Well indices, number of wells, coordinates
2. Set up producers and injectors separately


6,7) Permeability and porosity plots


8) Initial Injector and producer rates
1. Specify initial injection rate and production bhp

9) Well control optimization
1. Define bounds



Big picture

1) initialize variable, set stop criteria

2) initialize objective function, gradient, hessian
    1. function for calculation of objective function and its gradients

3) apply optimization in loop
  1. Update control variable using gradient descent (use function calls