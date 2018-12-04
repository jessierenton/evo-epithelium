#Code for Voronoi Tessellation model

##To generate data on fixation probabilities:
- without migration run `python run_parallel_no_migration.py`
	- set the `nlist_file` variable to equal the file containing the neighbour sets for vertices on the graph, e.g. hex_graph or vt_graph_6 (make sure these files which can be found in the graphs folder are in the folder containing the run file)
    - output saved to EGTpd_av_db/no\_migration/nlist_file 
- with migration run `python run_parallel_no_migration.py`
    - set the `nlist_file` variable as above
    - output saved to EGTpd_av_db/migration/nlist_file/m\_X where X is values of m to 1dp

##To calculate fixation probabilities according to equation in Allen et al 2016 'Evolution of cooperation on any population structure':
- copy files from 'graphs' into main folder and run `python find_coalescence_times.py`
- saves fixprobs to EGTpd_av_db/allen_result/
- prints critical ratios to terminal
