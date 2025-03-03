# Event Constrained Programming Case Studies
Here we provide the source code to the results presented in 
"Event Constrained Programming" by 
Daniel Ovalle, Stefan Mazzadi, Carl D. Laird, Ignacio E. Grossmann, and Joshua L. Pulsipher. 

![iterations](Pandemic_Control_Files/Figures/paper_indicator_iterations_0.9.png)

## Running the Code
All the results from paper were obtained using `ieee14.py` which uses 
[Pyomo](http://www.pyomo.org/) in combination with the Gurobi solver. 
The plots are generated using plot_results.py. These can be run as 
normal Python scripts. `ieee14.py` generates the data files stored in 
the `data/` folder and `plot_results.py` uses these data files to 
summarize our findings.
