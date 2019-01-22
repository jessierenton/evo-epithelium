##Processing data 
The directories MU-*/raw\_data are obtained from VTdyn by running `cluster\_stats\_vary\_MU.py`. To calculate lambda_CC and fixation probabilities (equation 18 in paper) follow the following steps

- run `process_data.py` 
	- produces MU-*\ints\_CC/ints\_CD/wints\_CC/wints\_CD folders each containing files n\_i for cluster sizes i=1 to 99.
		- ints\_CC/CD files contain number of cooperator-cooperator/defector interactions for each n
		- wints\_CC/CD files contain same for weighted interactions
- run `calc_fixprobs.py` to produce files 
	- MU-*/Lambda\_CC (which saves Lambda_CC values for n=1...100 for each MU value)
	- MU-*\_fixprobs (which saves fixation probs according to equation 18 for each MU value)

To produce Figures 2 and 3 in supplementary information run `plot.py`