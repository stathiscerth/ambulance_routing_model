# ambulance_routing_model
Comparing the Dynamic Casualty Transportation, optimization model, with immediate first rule, and providing ambulance routing for both. Manual changes into codes first section can simulate different mci scenarios. 
To run the simulations, after downloading the python files and createing a python project:
a) download  glpk solver (https://sourceforge.net/projects/winglpk/) and use the paths at main_optimizer lines 19 and 21. create a folder named glpkzip containing the downladed solver folder. add  this folder into the python projects folder
b) use the requirements file to create and activate a virtual environment, containing pyomo module
c) run the MAIN_optimizer file...edit manually the file to simulate different cases of parallel mci events.
