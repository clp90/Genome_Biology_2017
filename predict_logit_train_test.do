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

Colette L Picard
03/11/2017

This scripts takes as input a dataset containing data for all RILs on which sites
gained or lost methylation in those RILs, RIL parent of origin, local CpG, CHG and CHH
methylation, [C], [CG], [CHG] and [CHH], FPKM values, and motif occurrences. For each
RIL, it keeps the relevant predictors, drops all data with missing values in predictors,
and gets a subset of successes and failures that are >200bp apart at minimum.

Then runs a logit and estat classification to evaluate models. For each model, random
sampling is performed 1000x and accuracy of predictions stored. Different models (nested)
are tried for best model prediction.

Note: to run, you must cd to the directory containing full_dataset.txt (from Additional file 3)
or add full path to full_dataset.txt in line 44.

Then, to run, simply run:

do predict_logit_train_test.do

assuming predict_logit_train_test.do is in your current directory.

--------------------------------------------------------------------------------
*/

	version 14

	local maxiter=10

* Store results here

	clear
	gen ril = .
	gen parent = ""
	gen dir = ""
	gen model = ""
	gen numobs = .
	gen p_correct = .
	gen sensitivity = .
	gen specificity = .
	tempfile results
	qui save `results'
	clear

* Load the full dataset

	use  "full_dataset.txt"
	gen logfpkm_col = log(fpkm_col)
	gen logfpkm_cvi = log(fpkm_cvi)
	encode variability, gen(varcode)
			
* Save full version, then loop over all the RILs

	tempfile fulldata
	qui save `fulldata'
	clear
	
	tempfile training
	tempfile testing
	tempfile success
	tempfile background
	
	foreach ril in 8 22 84 124 242 258 303 332 363 495 {
		foreach parent in "col" "cvi" {
			if "`parent'" == "col" {
				local nonparent = "cvi"
			}
			if "`parent'" == "cvi" {
				local nonparent = "col"
			}
				
			// single var models
			local model1 ril`ril'_cpg
			local model2 ril`ril'_chg
			local model3 ril`ril'_chh
			local model4 `nonparent'hime
			local model5 gc_content
			local model6 frac_cpg_sites
			local model7 frac_chh_sites
			local model8 logfpkm_`parent'
			local model9 occ_tgcwr
			local model10 occ_rcatw
			local model11 varcode
			local model12 col_srnas
			local model13 cvi_srnas
			
			// compound models
			local model14 ril`ril'_cpg ril`ril'_chg ril`ril'_chh
			local model15 gc_content frac_cpg_sites frac_chh_sites occ_tgcwr occ_rcatw			
			local model16 col_srnas cvi_srnas
			local model17 `nonparent'hime varcode
			local model18 ril`ril'_cpg `nonparent'hime
			local model19 ril`ril'_cpg `nonparent'hime varcode
			local model20 ril`ril'_cpg `nonparent'hime col_srnas cvi_srnas
			local model21 ril`ril'_cpg `nonparent'hime varcode col_srnas cvi_srnas
			local model22 ril`ril'_cpg ril`ril'_chg ril`ril'_chh logfpkm_`parent'
			local model23 ril`ril'_cpg ril`ril'_chg ril`ril'_chh col_srnas cvi_srnas
			local model24 gc_content frac_cpg_sites frac_chh_sites occ_tgcwr occ_rcatw varcode		
			local model25 gc_content frac_cpg_sites frac_chh_sites occ_tgcwr occ_rcatw col_srnas cvi_srnas
			local model26 ril`ril'_cpg ril`ril'_chg ril`ril'_chh `nonparent'hime gc_content frac_cpg_sites frac_chh_sites occ_tgcwr occ_rcatw logfpkm_`parent' varcode col_srnas cvi_srnas

			// recreate all models but with indicator var for varcode if necessary
			local imodel1 ril`ril'_cpg
			local imodel2 ril`ril'_chg
			local imodel3 ril`ril'_chh
			local imodel4 `nonparent'hime
			local imodel5 gc_content
			local imodel6 frac_cpg_sites
			local imodel7 frac_chh_sites
			local imodel8 logfpkm_`parent'
			local imodel9 occ_tgcwr
			local imodel10 occ_rcatw
			local imodel11 i.varcode
			local imodel12 col_srnas
			local imodel13 cvi_srnas
			local imodel14 ril`ril'_cpg ril`ril'_chg ril`ril'_chh
			local imodel15 gc_content frac_cpg_sites frac_chh_sites occ_tgcwr occ_rcatw			
			local imodel16 col_srnas cvi_srnas
			local imodel17 `nonparent'hime i.varcode
			local imodel18 ril`ril'_cpg `nonparent'hime
			local imodel19 ril`ril'_cpg `nonparent'hime i.varcode
			local imodel20 ril`ril'_cpg `nonparent'hime col_srnas cvi_srnas
			local imodel21 ril`ril'_cpg `nonparent'hime i.varcode col_srnas cvi_srnas
			local imodel22 ril`ril'_cpg ril`ril'_chg ril`ril'_chh logfpkm_`parent'
			local imodel23 ril`ril'_cpg ril`ril'_chg ril`ril'_chh col_srnas cvi_srnas
			local imodel24 gc_content frac_cpg_sites frac_chh_sites occ_tgcwr occ_rcatw i.varcode		
			local imodel25 gc_content frac_cpg_sites frac_chh_sites occ_tgcwr occ_rcatw col_srnas cvi_srnas
			local imodel26 ril`ril'_cpg ril`ril'_chg ril`ril'_chh `nonparent'hime gc_content frac_cpg_sites frac_chh_sites occ_tgcwr occ_rcatw logfpkm_`parent' i.varcode col_srnas cvi_srnas

			// loop through all regions and RILs to test each model
			di "Looking at CGs in RIL`ril' with genotype `parent' in `region'"
			forvalues m = 1/26 {
				di "Testing model`m'"
				
				clear
				use `fulldata'
				keep if gene == 1 & te == 0
				
				replace parent_`ril' = lower(parent_`ril')
				keep if parent_`ril' == "`parent'"
				keep chr-end `parent'hime `parent'lome `nonparent'hime `parent'v`ril' ril`ril'_cpg ril`ril'_chg ril`ril'_chh gc_content-frac_chh_sites logfpkm_`parent' fpkm_`parent' occ_tgca-occ_catw varcode col_srnas cvi_srnas
			
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

				// look at gains and losses separately
				drop if `parent'v`ril' == .
				gen RILgain = `parent'v`ril' == 1
				tempfile ftemp
				qui save `ftemp'
			
				// start with gains:
				// keep all successes, then sample such that no two obs are within 200bp
				di "Getting all occurrences where RIL`ril' gained methylation relative to `parent'"
				keep if RILgain
				qui count
				local successcount=`r(N)'
				di "`successcount' successes kept"

				qui save `success', replace
				clear
			
				// sample background obs such that no two obs within 200bp, then get subset matched to prev
				use `ftemp'
				keep if `parent'lome & ! RILgain
				qui save `background', replace
				clear

				forvalues i = 1/`maxiter' {
					di "Running iter `i'"
					use `background'
					gen u = runiform()
					sort u
					keep if _n <= `successcount'
					drop u
			
					// append back to successes -> training set
					append using `success'
					qui save `training', replace
					count
					local numobs = `r(N)'
					
					// repeat to get test set
					clear
					use `background'
					gen u = runiform()
					sort u
					keep if _n <= `successcount'
					
					drop u
					append using `success'
					qui save `testing', replace
					clear
					
					capture noisily {
						use `training'
						qui logit RILgain `imodel`m'', iterate(100)
						
						use `testing', clear
						qui estat classification, all

						clear
						set obs 1
						gen ril = `ril'
						gen parent = "`parent'"
						gen dir = "RILgain"
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
								
				// repeat for RILloss
				use `ftemp'
				drop if `parent'v`ril' == .
				gen RILloss = `parent'v`ril' == -1
				qui save `ftemp', replace
								
				di "Getting all occurrences where RIL`ril' lost methylation relative to `parent'"
				keep if RILloss
				qui count
				local successcount=`r(N)'
				di "`successcount' successes kept"

				qui save `success', replace
				clear
			
				// sample background obs such that no two obs within 200bp, then get subset matched to prev
				use `ftemp'
				keep if `parent'hime & ! RILloss
				qui save `background', replace
				clear

				forvalues i = 1/`maxiter' {
					di "Running iter `i'"
					use `background'
					gen u = runiform()
					sort u
					keep if _n <= `successcount'
					drop u
			
					// append back to successes -> training set
					append using `success'
					qui save `training', replace
					count
					local numobs = `r(N)'
					
					// repeat to get test set
					clear
					use `background'
					gen u = runiform()
					sort u
					keep if _n <= `successcount'
					drop u
					append using `success'
					qui save `testing', replace
					clear
					
					capture noisily {
						use `training'
						qui logit RILloss `imodel`m'', iterate(100)
						
						use `testing', clear
						qui estat classification, all

						clear
						set obs 1
						gen ril = `ril'
						gen parent = "`parent'"
						gen region = "`region'"
						gen dir = "RILloss"
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
		}
	}
	
	use `results'
	save "regression_results.dta", replace
	
	exit 1
	
	
