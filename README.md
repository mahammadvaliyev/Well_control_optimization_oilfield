# Well_control_optimization_oilfield
Well control optimization of an oilfield operated under waterflood

## Summary

Significant fraction of oil is left in reservoirs when utilizing only primary recovery methods. This necessicates the application of the secondary recovery methods such as water or gas injection. However, even after application of secondary recovery methods, still a large volume of oil remains in the subsurface reservoirs, leaving a room for application of novel methods to increase the recovery. Among many possible solutions, well control optimization is relatively cheap, widely applicable option to use for further optimization of the recovery. In this paper, well control optimization was applied one layer of SPE 10 model. Numerical model was simulated for 10 years and optimization frequency for controls was set to be once per year. Optimization was performed using interior point algorithm with bound constraints on well controls coupled with MRST reservoir simulator and adjoint formulation for gradient calculation. Using NPV as a performance measure, an increase of 9.4% was observed after application of optimization compared to a constant control strategy.

## Workflow components

![workflow_components](https://user-images.githubusercontent.com/68789630/147499479-7a241419-f66f-475b-b833-d19d34aa60eb.jpg)

## Reservoir model

![model_details](https://user-images.githubusercontent.com/68789630/147499508-757e6916-00e9-4e97-b71c-0d54cd497b63.jpg)

## Results

![results](https://user-images.githubusercontent.com/68789630/147499547-3761f0c4-ad43-46f4-ae18-6094fb7a043a.jpg)

## Files
- project presentation
- project report
- Matlab files (main and auxiliary functions for running project code). Note that you need 1) Matlab and 2) MRST to run the code
