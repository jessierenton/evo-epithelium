##Processing data 
The directory raw\_data is obtained from VTdyn. To calculate lambda_CC and fixation probabilities (equation 18 in paper) follow the following steps

- run `process_data.py` 
	- produces ints\_CC/ints\_CD/wints\_CC/wints\_CD folders each containing files n\_i for cluster sizes i=1 to 99.
		- ints\_CC/CD files contain number of cooperator-cooperator/defector interactions for each n
		- wints\_CC/CD files contain same for weighted interactions
- run `calc_fixprobs.py` to produce files 
	- Lambda\_CC (which saves Lambda_CC values for n=1...100)
	- VTpd\_av\_decoupled\_theory (which saves fixation probs according to equation 18)

To produce Figure 4 in text run `plot_distributions.py`