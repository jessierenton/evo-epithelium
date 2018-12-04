#Code for Voronoi Tessellation model

##To generate data on fixation probabilities:
- for decoupled update run `python multirun_constant_N.py av 0 10 X_1 X_2 ... `
	- X\_i are values of b up to 2 decimal places; av is for averaged payoffs; 0/10 are start/end batch numbers, i.e. this will run 10 batches x1000 simulations
	- output saved in directory VTpd\_av\_decoupled
	- data in files fixX (X=b to 2dp) for decoupled update
	- data in form #fixed #lost #incomplete (sum=1000 across rows, each row is one batch)
- for death-birth update run `python multirun_death_birth.py av 0 10 X_1 X_2 ... `
	- same command line args and file structure as above
	- output saved in directory VTpd\_av\_db


##To generate interaction data:
- run `python cluster_stats.py 5000` 
	- command line arg is number of simulations
	- produces interaction\_data/raw\_data folder
	- see data/interaction\_data for processing and how to calculate Lambda_CC and fixation probabilities


##To plot tissue objects or create animation
- need to import libs.plot as vplt
- basic plot for tissue object (with torus geometry) use vplt.torus_plot(tissue) 
- to view animation use vplt.animate_torus(tissue)
- to save animation use vplt.save_mpg_torus(history,filename)
- see file libs/plot.py for optional arguments
