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

This scripts takes a dataset containing seven variables - CG location (chr, start, end),
local CG, CHG and CHH methylation, and classification of CG, from maize methylation
data from Li et al. 2015, and attempts to use the methylation data to predict which
CGs are classified as "variable" using a logistic regression approach.

--------------------------------------------------------------------------------
*/

	version 14

	local maxiter=10

* Store results here

	clear
	gen strain = ""		// which strain was used
	gen dir = ""		// whether variable sites were distinguished from uniform lo or hi CGs)
	gen model = ""		// model number (see below)
	gen numobs = .		// number of observations
	gen p_correct = .	// percent of CGs correctly assigned
	gen sensitivity = .
	gen specificity = .
	tempfile results
	qui save `results'
	clear

* Process data from each strain separately

	// define all models (all possible combinations of the 3 predictors)
	local model1 B73_CpG
	local model2 B73_CHG
	local model3 B73_CHH
	local model4 B73_CpG B73_CHG
	local model5 B73_CpG B73_CHH
	local model6 B73_CHG B73_CHH
	local model7 B73_CpG B73_CHG B73_CHH
	
	// load data from this strain
	insheet using "full_maize_B73_data.txt", case
				
	// get subset of sites 200bp apart
	clonevar pos = start
	bysort chr (start): replace pos = pos[_n-1] if pos - pos[_n-1] < 200
	bysort chr pos (start): gen pick = _n == 1
	keep if pick
	drop pick pos
	
	tempfile ftemp
	qui save `ftemp', replace
	
	tempfile training
	tempfile testing
	tempfile success
	tempfile background

	// loop through all regions and RILs to test each model
	di "Looking at CGs in strain B73"
	forvalues m = 1/7 {
		di "Testing model`m'"
		di "`model`m''"
		use `ftemp'
		
		// first, test ability to distinguish variable sites from unimodal lo sites
		keep if class == "variable"
		gen isVariable = 1
		qui count
		local successcount=`r(N)'
		di "`successcount' successes kept"

		qui save `success', replace
		clear
	
		// test ability to distinguish from unimodal hi
		use `ftemp'
		keep if class == "unimodal_hi"
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
				gen strain = "B73"
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
		keep if class == "unimodal_lo"
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
				gen strain = "B73"
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
	save "regression_results_maize_B73.dta", replace
	
	exit 1
	
	
