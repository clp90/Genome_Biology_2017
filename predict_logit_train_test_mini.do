/*
--------------------------------------------------------------------------------

Copyright 2017 Colette L Picard

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

03/11/2017
Colette L Picard

This script is meant to run a similar, restricted analysis to predict_logit_train_test_maize.do
on the thaliana data. Analysis requires the same full dataset as predict_logit_train_test.do
(full_dataset.txt) and is run separately on each RIL.

--------------------------------------------------------------------------------
*/

	version 14

	local maxiter=10

* Store results here

	clear
	gen ril = ""
	gen dir = ""
	gen model = ""
	gen numobs = .
	gen p_correct = .
	gen sensitivity = .
	gen specificity = .
	tempfile results
	qui save `results'
	clear

* Load the full dataset, drop everything that won't be used

	insheet using  "full_dataset.txt"	
	keep chr-end ril8_cpg-ril495_chh variability	
	keep if variability == "variable" | variability == "unimodal_lo" | variability == "unimodal_hi" 
	
* Save full version, then loop over all the RILs

	tempfile fulldata
	qui save `fulldata'
	clear
	
	foreach ril in 8 22 84 124 242 258 303 332 363 495 {
			
		// define all models (all possible combinations of the 3 predictors)
		local model1 ril`ril'_cpg
		local model2 ril`ril'_chg
		local model3 ril`ril'_chh
		local model4 ril`ril'_cpg ril`ril'_chg
		local model5 ril`ril'_cpg ril`ril'_chh
		local model6 ril`ril'_chg ril`ril'_chh
		local model7 ril`ril'_cpg ril`ril'_chg ril`ril'_chh

		tempfile training
		tempfile testing
		tempfile success
		tempfile background
	
		// loop through all regions and RILs to test each model
		di "Looking at CGs in strain `strain'"
		forvalues m = 1/7 {
			di "Testing model`m'"
			di "`model`m''"
			use `fulldata'
			keep chr start end variability `model`m''

			// drop any rows with missing values
			di "Dropping all rows with missing values"
			gen mis = 0
			di "`model`m''"
			foreach v of varlist `model`m'' {
				replace mis = mis + 1 if `v' == .		
			}
			tab mis
			drop if mis > 0
			drop mis
			count
			
			// get subset of sites 200bp apart
			clonevar pos = start
			bysort chr (start): replace pos = pos[_n-1] if pos - pos[_n-1] < 200
			bysort chr pos (start): gen pick = _n == 1
			keep if pick
			drop pick pos
		
			tempfile ftemp
			qui save `ftemp', replace
				
			// first, test ability to distinguish variable sites from unimodal lo sites
			keep if variability == "variable"
			gen isVariable = 1
			qui count
			local successcount=`r(N)'
			di "`successcount' successes kept"
			
			qui save `success', replace
			clear
		
			// test ability to distinguish from unimodal hi
			use `ftemp'
			keep if variability == "unimodal_hi"
			qui save `background', replace
			clear

			di "Testing classification using unimodal_hi as background"
			forvalues i = 1/`maxiter' {
				di "Running iter `i'"
				use `background'
				gen u = runiform()
				sort u
				keep if _n <= `successcount'
				drop u
				gen isVariable = 0
		
				// append back to successes -> training set
				append using `success'				
				qui save `training', replace
				count
				local numobs = `r(N)'
				
				// repeat to get test set
				di "Drawing background"
				clear
				use `background'
				gen u = runiform()
				sort u
				keep if _n <= `successcount'
				
				drop u
				gen isVariable = 0
				append using `success'
				qui save `testing', replace
				clear
				
				capture noisily {
					di "Training"
					use `training'
					qui logit isVariable `model`m'', iterate(100) difficult
					
					di "Testing"
					use `testing', clear
					qui estat classification, all

					di "Done"
					clear
					set obs 1
					gen ril = "`ril'"
					gen dir = "vs_unimodal_hi"
					gen model = "model`m'"
					gen numobs = `numobs'
					gen p_correct = `r(P_corr)'
					gen sensitivity = `r(P_p1)'
					gen specificity = `r(P_n0)'
					append using `results'
					qui save `results', replace
				}
				clear
			}

			// test ability to distinguish from unimodal lo
			use `ftemp'
			keep if variability == "unimodal_lo"
			qui save `background', replace
			clear

			di "Testing classification using unimodal_lo as background"
			forvalues i = 1/`maxiter' {
				di "Running iter `i'"
				use `background'
				gen u = runiform()
				sort u
				keep if _n <= `successcount'
				drop u
				gen isVariable = 0
		
				// append back to successes -> training set
				append using `success'
				qui save `training', replace
				count
				local numobs = `r(N)'
				
				// repeat to get test set
				di "Drawing background"
				clear
				use `background'
				gen u = runiform()
				sort u
				keep if _n <= `successcount'
				
				drop u
				gen isVariable = 0
				append using `success'
				qui save `testing', replace
				clear
				
				capture noisily {
					di "Training"
					use `training'
					qui logit isVariable `model`m'', iterate(100)
					
					di "Testing"
					use `testing', clear
					qui estat classification, all

					di "Done"
					clear
					set obs 1
					gen ril = "`ril'"
					gen dir = "vs_unimodal_lo"
					gen model = "model`m'"
					gen numobs = `numobs'
					gen p_correct = `r(P_corr)'
					gen sensitivity = `r(P_p1)'
					gen specificity = `r(P_n0)'
					append using `results'
					qui save `results', replace
					
				}
				clear
			}
			

		}		
		use `results'
		save "regression_results_mini_`ril'.dta", replace
		clear
	
	}
	
	use `results'
	save "regression_results_mini.dta", replace
	
	
	exit 1
	
	
	
