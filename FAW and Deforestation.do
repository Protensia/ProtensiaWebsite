************************************************************
* Author: Protensia Hadunka
* Project: Charcoal Production and Fall Armyworm (FAW)
* Purpose: Cleaning and data analysis 
************************************************************

use "C:\Users\hadunka2\Box\HICPS Cleaning 07_12_19\PNAS Paper",clear
clear matrix
clear mata

************************************************************
* SECTION 1: Yield and farm size cleaning
************************************************************

* Total maize production = harvested + left in the field

gen Total = qharvested + qleft

* Calculate total farm size across 5 planting plots
egen farmsize_W01 = rowtotal(plot_1 plot_2 plot_3 plot_4 plot_5), miss
label var farmsize "Maize farm size: ha"

* Drop extreme farm sizes (above 8 ha)
replace farmsize_W01 = . if farmsize_W01 > 8

* Log-transform farm size (raw and winsorized)
gen lnfarmsize = log(farmsize)
label var lnfarmsize "Log of maize farm size"

gen lnfarmsize_W01 = log(farmsize_W01)
label var lnfarmsize_W01 "Log of maize farm size (winsorized at 1%)"


************************************************************
* SECTION 2: Clean farmland data and impute outliers
************************************************************

* Create clean copy of farmland
gen farmland_clean = farmland

* Replace large farmland values with camp-level medians
bys camp: egen med2_farmland = median(farmland)
replace farmland_clean = med2_farmland if farmland > 13

* Compare and adjust if cleaned size is less than reported
count if farmland_clean < farmsize_W01
list farmland_clean farmsize_W01 if farmland_clean < farmsize_W01
bys year: count if farmland_clean < farmsize_W01 & farmsize_W01 != .
replace farmland_clean = farmsize_W01 if farmland_clean < farmsize_W01 & farmsize_W01 != .

* Label cleaned variable
label var farmland_clean "Total area of farmland: Ha"

gen farmsize_clean = farmsize_W01


************************************************************
* SECTION 3: Input variables - seed and fertilizer
************************************************************

* Total maize seed used across plots
foreach var of varlist qseed_1 qseed_2 qseed_3 qseed_4 qseed_5 {
    quietly recode `var' (. = 0)
}

gen quaseed = qseed_1 + qseed_2 + qseed_3 + qseed_4 + qseed_5
rename quaseed qseed

* Total top-dressing fertilizer used
foreach var of varlist qtop_1 qtop_2 qtop_3 qtop_4 qtop_5 {
    quietly recode `var' (. = 0)
}
gen qfer = qtop_1 + qtop_2 + qtop_3 + qtop_4 + qtop_5
rename qfer qtop

***************************************************************
* SECTION 4: Create Fertilizer Variables and Basic Aggregates
***************************************************************

* Aggregate quantities of basal fertilizer into a single variable
foreach var of varlist qbasal_1 qbasal_2 qbasal_3 qbasal_4 qbasal_5 {
    quietly recode `var' (. = 0)
}
gen quafert = qbasal_1 + qbasal_2 + qbasal_3 + qbasal_4 + qbasal_5
rename quafert qbasal 

* Combine basal and top dressing fertilizers
gen fert = qtop + qbasal

* Encode the 'camp' variable as numeric
encode camp, gen(camp1)
gen camp2 = camp1  // Duplicate for robustness or alternate coding schemes

***************************************************************
* SECTION 5: Winsorize Key Variables by Year to Handle Outliers
***************************************************************

winsor2 Total, cuts(1 99) by(year)
winsor2 cultivown_land, cuts(1 99) by(year)


* Generate yield variable: output (Total) divided by cultivated land area (cultivown_land)
*drop yield
gen yield = Total_w / cultivown_land_w

label variable yield "Maize yield (Output per hectare)"

* Step 1: Calculate the median excluding zeros
summarize cultivown_land_w if cultivown_land_w > 0, detail
local med = r(p50)

* Step 2: Replace zero values with the calculated median
replace cultivown_land_w = `med' if cultivown_land_w == 0

summarize hh_head_age, detail
scalar median_age = r(p50)
* Step 2: Replace values below 15 with the median
replace hh_head_age = median_age if hh_head_age < 15

***************************************************************
* SECTION 6: Creating the Total Income Variable
* Purpose: Aggregate all income sources to compute total household income
***************************************************************

* Generate total household income by summing across all income sources
gen Total_income = income_piecework + income_salary + income_smallbusiness + ///
                   income_charcoal + income_gardening + income_forestproduct + ///
                   income_livestock + income_remittance + income_other

label variable Total_income "Total household income (sum of all income sources)"


* Purpose: Generate a binary variable indicating household access to credit
***************************************************************

* Initialize access to credit variable (0 = No access, 1 = Has access)
gen Accs_credit = 0
label variable Accs_credit "Household has access to credit (binary indicator)"

* Update Accs_credit to 1 if household borrowed any amount
replace Accs_credit = 1 if borrow500 == 1 | borrow2500 == 1 | borrow10000 == 1

***************************************************************
* SECTION 7: Instrumental Variables (IV) for Intensity of Armyworm Response
* Purpose: Construct IVs at camp level based on household reports of severity
***************************************************************

* Step 1: Generate binary indicators for severity levels

* Low severity (army_aff == 1)
gen low = (army_aff == 1)
label variable low "Armyworm - Low severity (binary)"

* Medium severity (army_aff == 2)
gen med = (army_aff == 2)
label variable med "Armyworm - Medium severity (binary)"

* Severe (army_aff == 3)
gen seve = (army_aff == 3)
label variable seve "Armyworm - Severe severity (binary)"


* Step 2: Sum severity indicators at camp-year level

* Sum of low severity households in each camp-year
egen iv_camp = sum(low), by(camp year)

* Sum of medium severity households in each camp-year
egen iv_camp2 = sum(med), by(camp year)

* Sum of severe severity households in each camp-year
egen iv_camp3 = sum(seve), by(camp year)


* Step 3: Create household-minus-self IV variables

* Exclude the household's own severity from camp sums
gen hh_iv_low  = iv_camp  - low
gen hh_iv_med  = iv_camp2 - med
gen hh_iv_seve = iv_camp3 - seve


* Step 4: Number of households reporting each severity level

egen num_in_camp1 = sum(low), by(camp year)
egen num_in_camp2 = sum(med), by(camp year)
egen num_in_camp3 = sum(seve), by(camp year)

* Subtract one to exclude the household itself
gen camp_1 = num_in_camp1 - 1
gen camp_2 = num_in_camp2 - 1
gen camp_3 = num_in_camp3 - 1

* Total number of households in each camp-year
egen number_in_camp = count(HHID), by(camp year)

* Exclude the household itself
gen num_1 = number_in_camp - 1

* Step 5: Final camp-level IVs for each severity


gen camp_IV1 = camp_1 / num_1
gen camp_IV2 = camp_2 / num_1
gen camp_IV3 = camp_3 / num_1

label variable camp_IV1 "IV: Share of other households - Low severity"
label variable camp_IV2 "IV: Share of other households - Medium severity"
label variable camp_IV3 "IV: Share of other households - Severe severity"


* Step 6: Alternative corrected IVs


gen IV1 = (sum(low) - low) / num_1
gen IV2 = (sum(med) - med) / num_1
gen IV3 = (sum(seve) - seve) / num_1

* Step 7: Binary armyworm incidence (regardless of severity)


egen num_in_camp4 = sum(army), by(camp year)
gen camp_4 = num_in_camp4 - 1
gen camp_IV4 = camp_4 / num_1

* Corrected IV for binary armyworm
gen IV4 = (sum(army) - army) / num_1
label variable camp_IV4 "IV: Share of other households - Armyworm incidence"

* Step 8: Camp-level sums of severity and incidence


egen sum_camp  = sum(army_aff), by(camp year)
egen sum_camp2 = sum(army), by(camp year)

* Household value minus camp mean (severity)
gen hh_minus_camp_mean  = sum_camp  - army_aff
gen hh_minus_camp_mean2 = sum_camp2 - army

* Step 9: IHS Transformation of Differences

* Apply IHS transformation to difference in severity
ihstrans hh_minus_camp_mean

* Rename transformed variable
rename ihs_hh_minus_camp_mean army_diff
label variable army_diff "IHS of deviation from camp mean severity"

* Step 10: Alternative camp IVs using army incidence


egen number_in_camp2 = count(army), by(camp year)
gen camp_one  = number_in_camp - 1
gen camp_one2 = number_in_camp2 - 1

* Camp-level IVs for army incidence
gen camIV  = hh_minus_camp_mean / camp_one
gen camIV2 = hh_minus_camp_mean2 / camp_one2

* Additional IV alternative
gen cam3 = hh_minus_camp_mean2 / num_1

label variable camIV  "Camp IV - deviation in severity"
label variable camIV2 "Camp IV - deviation in incidence"

***************************************************************
* SECTION 8: Construction of Household Asset Indices
* Purpose: Generate asset indices using PCA and ILRI weighting methods
***************************************************************

* Reference for PCA method:
* https://jbhender.github.io/Stats506/F17/Projects/G18.html

***************************************************************
* PART A: Principal Component Analysis (PCA) Approach
***************************************************************

* Step 1: Inspect key asset variables

browse oxen breeding_bull donkey female_cattle_number goat_sheep_number ///
       poultry_number pigs_number asset_phone tv radio bike motorcycle ///
       ox_carts vehicle water_pump plough sprayers

* Step 2: Prepare variables for PCA

* Encode pigs_number if it is not already numeric
encode pigs_number, gen(pigs_count)

* Step 3: Conduct PCA on asset variables

pca oxen breeding_bull donkey female_cattle_number goat_sheep_number ///
    poultry_number pigs_count asset_phone tv radio bike motorcycle ///
    ox_carts vehicle water_pump plough sprayers

* Notes:
* - Component 1 explains ~25.9% of the variation
* - Components 1-13 together explain ~90% of the variation

* Step 4: Generate PCA-based Asset Index

* Predict first principal component as asset index
predict assetind_pca
label variable assetind_pca "Household Asset Index (PCA method)"

* Summarize asset index
summarize assetind_pca
tabulate assetind_pca

* Plot scree plot with confidence intervals
screeplot, ci(asympt level(95)) mean scheme(lean2)

* Plot factor loadings
loadingplot

***************************************************************
* PART B: ILRI Weighted Index Approach
***************************************************************

* Step 1: Prepare asset variables

* Replace missing values with zero for asset ownership variables
foreach v in oxen breeding_bull donkey female_cattle_number goat_sheep_number ///
              poultry_number pigs_count asset_phone tv radio bike motorcycle ///
              ox_carts vehicle water_pump plough sprayers {
    replace `v' = 0 if missing(`v')
}


* Step 2: Construct Asset Sub-Indices following ILRI Method

* Animal Asset Index
gen animalind = 10*(oxen + breeding_bull + donkey + female_cattle_number) ///
                + 3*goat_sheep_number + 1*poultry_number + 2*pigs_count
label variable animalind "Animal Asset Index"

* Domestic Asset Index
gen domesticind = 3*asset_phone + 4*tv + 2*radio
label variable domesticind "Domestic Asset Index"

* Transport Asset Index
* Note: ox carts included in transport assets based on ILRI documentation
gen transportind = 6*bike + 48*motorcycle + 12*ox_carts + 160*vehicle
label variable transportind "Transport Asset Index"

* Productive Asset Index
* Note: Water pumps correspond to treadle pumps; sprayers treated as ploughs
gen productiveind = 6*water_pump + 4*(plough + sprayers)
label variable productiveind "Productive Asset Index"


* Step 3: Construct Total Household Asset Indices


* Full ILRI method asset index
gen assetind_ilri = domesticind + animalind + transportind + productiveind
label variable assetind_ilri "Household Asset Index (ILRI method)"

* Alternative index (excluding domestic assets if missing for some years)
gen assetind = animalind + productiveind
label variable assetind "Alternative Household Asset Index"


************************************************************
* SECTION 9: Transformations for Normalization and Modeling
* Purpose: Normalize skewed variables using IHS and log transformations 
*          and prepare quadratic terms for nonlinear effects
************************************************************

* Rename totalrain to rainfall for consistency
rename totalrain rainfall

************************************************************
* INVERSE HYPERBOLIC SINE (IHS) TRANSFORMATIONS
* Used when log-transforming zero or negative values is not possible
************************************************************

* IHS transform of household size (to reduce right skewness)
ihstrans hh_num

* IHS transform of rainfall (used as a continuous explanatory variable)
ihstrans rainfall

* Generate squared term of raw rainfall for capturing nonlinear rainfall effects
gen sqrainfall = (rainfall)^2

* Generate squared term of IHS-transformed rainfall
gen lsqrainfall = (ihs_rainfall)^2

* IHS transform of temperature to manage skewness in climatic variables
ihstrans temperature

* Generate squared term of raw temperature to model curvature effects
gen sqtemp = (temperature)^2

* Generate squared term of IHS-transformed temperature
gen lsqtemp = (ihs_temperature)^2

* IHS transform of household income, accounting for outliers and zeros
ihstrans Total_income

* IHS transform of total farmland size
ihstrans farmland

* IHS transform of yield (e.g., maize yield per hectare)
ihstrans yield

* IHS transform of cultivated area owned by the household
ihstrans cultivown_land

* IHS transform of fertilizer use (helps manage zero users and skewness)
ihstrans fert

* IHS transform of quantity of seed used
ihstrans qseed

************************************************************
* LOG TRANSFORMATIONS (log(x)) FOR NON-ZERO VARIABLES
************************************************************

* Log of Food Consumption Score (FCS), often right-skewed
gen lnfcs = ln(fcs)
label variable lnfcs "Log of Food Consumption Score (FCS)"

* Generate log and squared log of rainfall to test diminishing or nonlinear returns
gen lnrain = ln(rainfall)
gen ln2_rainfall = (ln(rainfall))^2

* Generate log and squared log of temperature
gen lntemp = ln(temperature)
gen ln2_temperature = (ln(temperature))^2

* Log of squared temperature (alternative nonlinear transformation)
gen lnsqtemp = ln(sqtemp)

* Log of fertilizer use (zero fertilizer users should be treated separately)
gen lnfert = ln(fert)

* Log of seed quantity
gen lnqseed = ln(qseed)

* Log transformation of cultivated land owned (cultivown_land)
gen lncultivation = ln(cultivown_land_w)
label variable lncultivation "Log of Cultivated Land Area Owned"

************************************************************
* LOG TRANSFORMATIONS FOR FOOD SECURITY INDICATORS
************************************************************

* Log of FCS categorical (discrete food security status levels)
gen lnfcs_cat = ln(fcs_cat)
label variable lnfcs_cat "Log of Food Consumption Score (FCS) categorical"

* Log of Household Dietary Diversity Score
gen lhdds = ln(hdds)
label variable lhdds "Log of Household Dietary Diversity Score (HDDS)"

* Log of Coping Strategy Index (proxy for food insecurity stress)
gen lrcsi = ln(rcsi_cat)
label variable lrcsi "Log of Coping Strategy Index (CSI)"


* Food Consumption Score (FCS)
ihstrans fcs
rename ihs_fcs hpsfcs
label variable hpsfcs "Hyperbolic sine transformed Food Consumption Score"

* Food Consumption Score (FCS) - categorical score
ihstrans fcs_cat
rename ihs_fcs_cat hps_fcs_cat
label variable hps_fcs_cat "Hyperbolic sine transformed fcs (categorical score)"

* Household Dietary Diversity Score (HDDS)
ihstrans hdds
rename ihs_hdds hpshdds
label variable hpshdds "Hyperbolic sine transformed Household Dietary Diversity Score"

* Food Consumption Score (FCS) - categorical score
ihstrans hdds_cat
rename ihs_hdds_cat hps_hdds_cat
label variable hps_hdds_cat "Hyperbolic sine transformed hdds (categorical score)"

* Reduced Coping Strategies Index (rCSI) - raw score
ihstrans rcsi1
rename ihs_rcsi1 hpsrcsi
label variable hpsrcsi "Hyperbolic sine transformed rCSI (raw score)"

* Reduced Coping Strategies Index (rCSI) - categorical score
ihstrans rcsi_cat
rename ihs_rcsi_cat hps_rcsi
label variable hps_rcsi "Hyperbolic sine transformed rCSI (categorical score)"


***************************************************************
* SECTION 10: Setting Panel Data Structure
* Purpose: Sort and declare the dataset as panel data
***************************************************************

* Sort the dataset by household ID and year for panel data consistency
sort HHID year

* Declare the panel structure: HHID as the panel identifier, year as the time variable
xtset HHID year

***************************************************************
* SECTION 11: Baseline Descriptive Statistics and Data Preparation
* Purpose: Prepare variables and summarize key baseline characteristics (2016)
***************************************************************

***************************************************************
* PART A: Variable Preparation
***************************************************************

* Log-transform crop yield variable
gen logyield = ln(yield)
label variable logyield "Log of Crop Yield (Output per Hectare)"



* Calculate distance to forest based on time to collection (converted to years)
gen dist_forest = time_coll / 12
label variable dist_forest "Distance to Forest (Years, based on Time to Collection)"

* Replace household head age values below 15 with the median age (excluding values < 15)
* Step 1: Calculate the median age and store it in a scalar
summarize hh_head_age if hh_head_age >= 15, detail
scalar median_age = r(p50)

* Step 2: Replace values below 15 with the median
replace hh_head_age = median_age if hh_head_age < 15

***************************************************************
* SECTION 12: Baseline Summary Statistics (2016)
***************************************************************
*  Summarize Key Household Characteristics for Baseline Year

* Household head's age
summarize hh_head_age if year == 2016
* Household head's gender
summarize hh_head_edu if year == 2016
* Household head's education
summarize hh_head_sex if year == 2016
* Household size
summarize hh_num if year == 2016
* Charcoal production
summarize Charc if year == 2016
* Cultivated land size
summarize cultivown_land_w if year == 2016
* Total farmland owned
summarize farmland if year == 2016
* Maize or crop yield
summarize yield if year == 2016
*Total income_charcoal
summarize Total_income if year == 2016
* Rainfall (mm)
summarize rainfall if year == 2016
* Temperature (degrees)
summarize temperature if year == 2016
* Access to credit (binary indicator)
summarize Accs_credit if year == 2016
* Distance to forest 
summarize dist_forest if year == 2016
* Food consumption score 
summarize fcs_cat if year == 2016
* Household dietary diversity
summarize hdds_cat if year == 2016

***************************************************************
* SECTION 13: Analysis of Maize Yields and FAW (Fall Armyworm)
* Purpose: Estimate the effect of FAW incidence on maize yields
***************************************************************

*Log-Transformed Variables
*--------------------------------------------------------------

* 1. OLS regression: Direct effect of FAW (army_aff) on log yield
reghdfe logyield army_aff qseed fert cultivown_land_w ///
        rainfall sqrainfall temperature sqtemp, ///
        absorb(i.HHID i.year) cl(camp1)
outreg2 using results, excel replace 
estimates store yieldFE

* Notes:
* - ihs_yield is the outcome (Log-Transformed yield)
* - Key independent variable: army_aff (self-reported FAW incidence)
* - Household FE (HHID) and Year FE included via absorb()
* - Clustered standard errors at camp level (camp1)

* 2. OLS regression: Using camp-level IV (camIV)
reghdfe logyield camIV  qseed fert cultivown_land_w ///
        rainfall sqrainfall temperature sqtemp, ///
        absorb(i.HHID i.year) cl(camp1)
outreg2 using results, excel replace
estimates store yieldITT 

* Notes:
* - camIV = camp-level armyworm exposure excluding own household
* - Helps address endogeneity of self-reported FAW incidence

* 3. IV regression: Instrumenting army_aff with camp-level FAW incidence (camIV)
ivreghdfe logyield  qseed fert cultivown_land_w rainfall ///
          sqrainfall temperature sqtemp (army_aff = camIV), ///
          absorb(i.HHID i.year) cl(camp1) first
outreg2 using results, excel append
estimates store yieldIV

* Notes:
* - First stage: army_aff instrumented by camp-level FAW (camIV)
* - Second stage: impact on ihs_yield
* - Controls, fixed effects, and clustering same as above



***************************************************************
* SECTION 14: Analysis of Food Security Outcomes and FAW (Fall Armyworm)
* Objective: Quantify the impact of FAW exposure on household food security
*            using a combination of instrumental variable and control function approaches
* Outcomes analyzed:
*   - hpsfcs  : Household-level Food Consumption Score (FCS)
*   - hpshdds : Household Dietary Diversity Score (HDDS)
*   - hps_rcsi: Household Coping Strategy Index (CSI)
***************************************************************

***************************************************************
* 1. FOOD CONSUMPTION SCORE (FCS) MODELS
***************************************************************

* BASELINE MODEL: Tobit regression with Random Effects
* Key regressor: camIV — binary instrument for FAW exposure
* Rationale: Tobit is used due to censoring at zero; RE handles unobserved heterogeneity across households

xttobit hpsfcs army_aff rainfall sqrainfall temperature sqtemp, ll(0) re
margins, dydx(*) predict(ystar(0,.))
outreg2 using results, excel append
estimates store FCS_FE 


xttobit hpsfcs camIV rainfall sqrainfall temperature sqtemp, ll(0) re
margins, dydx(*) predict(ystar(0,.))
outreg2 using results, excel append
estimates store FCS_ITT 

* FIRST-STAGE MODEL: Estimate determinants of actual FAW exposure

* Fixed Effects control for time-invariant household-specific confounders
xtreg army_aff camIV rainfall sqrainfall temperature sqtemp, fe

* Generate residuals (vhat) from the first-stage regression
* These residuals represent the endogenous variation in army_aff
predict vhat, resid	

* IV MODEL (CONTROL FUNCTION APPROACH): Second-stage Tobit regression
* Includes vhat as a control to address endogeneity of army_aff
* Uses RE Tobit with censoring at 0 for hpsfcs
xttobit hpsfcs army_aff rainfall sqrainfall temperature sqtemp vhat, ll(0) re
margins, dydx(*) predict(ystar(0,.))
outreg2 using results, excel append
estimates store FCS_IV


***************************************************************
* 2. HOUSEHOLD DIETARY DIVERSITY SCORE (HDDS) MODELS
***************************************************************

* BASELINE MODEL: RE Tobit regression for HDDS
* Outcome: hpshdds
* Regressors: camIV + IHS-transformed rainfall and temperature (to reduce skewness)

xttobit hpshdds army_aff ihs_rainfall lsqrainfall ihs_temperature lsqtemp, ll(0)	
margins, dydx(*) predict(ystar(0,.))
outreg2 using results, excel append
estimates store HDDS_FE 

xttobit hpshdds camIV ihs_rainfall lsqrainfall ihs_temperature lsqtemp, ll(0)	
margins, dydx(*) predict(ystar(0,.))
outreg2 using results, excel append
estimates store HDDS_ITT 

* IV MODEL: RE Tobit with endogenous FAW exposure (army_aff)
* Includes vhat from first-stage regression to control for unobservables
xttobit hpshdds army_aff rainfall sqrainfall temperature sqtemp vhat, ll(0) re
margins, dydx(*) predict(ystar(0,.))
outreg2 using results, excel append
estimates store HDDS_IV


***************************************************************
* 3. COPING STRATEGY INDEX (CSI) MODELS
***************************************************************

* BASELINE MODEL: RE Tobit regression for CSI (proxy for food stress)
* Outcome: hps_rcsi
* Regressors: army_aff and IHS-transformed weather controls

xttobit hps_rcsi army_aff rainfall sqrainfall temperature sqtemp , ll(0)		   		  	   
margins, dydx(*) predict(ystar(0,.))
outreg2 using results, excel append

xttobit hps_rcsi camIV rainfall sqrainfall temperature sqtemp , ll(0)		   		  	   
margins, dydx(*) predict(ystar(0,.))
outreg2 using results, excel append
estimates store CSI_ITT


xttobit hps_rcsi army_aff rainfall sqrainfall temperature sqtemp, ll(0) re
margins, dydx(*) predict(ystar(0,.))
outreg2 using results, excel append
estimates store CSI_IV



***************************************************************
* SECTION 16: Cultivated Land Analysis
* Purpose: Estimate the impact of FAW (Fall Armyworm) incidence on cultivated land area
***************************************************************

***************************************************************
* PART A: Generate Log-Transformed Cultivated Land Variable
***************************************************************

* Purpose: Estimate whether FAW exposure affects future land-use decisions
gen forwardcultiv = F.lncultivation

***************************************************************
* PART B: Estimate Intent-to-Treat Effect of FAW on Future Cultivation
***************************************************************
reghdfe forwardcultiv army_aff rainfall sqrainfall temperature sqtemp, ///
       absorb(i.year i.HHID) cl(camp1)
	   
* 1: Reghdfe with Fixed Effects (Intent-to-Treat)
reghdfe forwardcultiv camIV rainfall sqrainfall temperature sqtemp, ///
       absorb(i.year i.HHID) cl(camp1)
outreg2 using results, excel append
estimates store Cult_ITT

***************************************************************
* REGRESSION 2: IV Regression — Addressing Endogeneity of FAW Exposure
***************************************************************
ivreghdfe forwardcultiv rainfall sqrainfall temperature sqtemp ///
          (army_aff = camIV), ///
          absorb(i.year i.HHID) cl(camp1)
outreg2 using results, excel append
estimates store Cult_IV


		
***************************************************************
* SECTION 16: Lagged FAW Incidence and Impact on Charcoal Production
* Objective: Estimate the effect of lagged FAW exposure on household charcoal production (Charc)
*            while accounting for unobserved heterogeneity and potential endogeneity
***************************************************************

***************************************************************
* Step 1: Generate lagged variables of interest
***************************************************************

* Create 1-period lag of FAW severity variable at the household-year level
* This captures the effect of past (rather than current) FAW incidence on current outcomes
gen l_FAW = L.army_aff
label variable l_FAW "Lagged FAW incidence (previous year)"

* Create 1-period lag of instrumental variable (IV) — camIV
gen l_camIV = L.camIV




***************************************************************
* Step 2: Estimate OLS effect of lagged FAW exposure
***************************************************************
quietly regress Charc l_FAW rainfall ///
                 sqrainfall temperature sqtemp

* Retain only estimation sample for downstream operations (to align samples)
keep if e(sample)

***************************************************************
* Step 3: Generate household-level means (CRE strategy)
* These time averages are used to implement the Correlated Random Effects (CRE) estimator
* Purpose: CRE allows us to control for household-level unobserved heterogeneity by conditioning on means
***************************************************************

foreach v in l_FAW  rainfall ///
            sqrainfall temperature sqtemp {
    * Create household-level mean of each time-varying regressor
    * This allows inclusion of within-household deviations and household-specific averages (Mundlak, 1978)
    bysort HHID: egen double mean_`v' = mean(`v')
}

* Compute marginal effects
mfx 

estimates store charcFE
outreg2 using results, excel replace


***************************************************************
* Step 2: Estimate Intent-to-Treat (ITT) effect of lagged FAW exposure
***************************************************************

reghdfe Charc l_camIV ///
        rainfall mean_rainfall ///
        sqrainfall mean_sqrainfall temperature mean_temperature ///
        sqtemp mean_sqtemp
estimates store charcITT



***************************************************************
* Step 4: Control Function / CRE Adjustment for Endogeneity
***************************************************************

* First-stage: regress lagged FAW incidence (l_FAW) on camIV and its mean + weather controls and their means
* This is the first stage of a control function strategy under CRE framework
* Purpose: Extract residual variation in l_FAW not explained by observables and the instrument
reghdfe l_FAW camIV  ///
        rainfall mean_rainfall ///
        sqrainfall mean_sqrainfall temperature mean_temperature ///
        sqtemp mean_sqtemp, resid

* Predict residuals from first-stage regression — these will capture endogeneity
* Inclusion of vhatCRE in second-stage regression adjusts for endogenous selection into FAW exposure
predict vhatCRE, resid	

***************************************************************
* Step 5: Second-stage regression — outcome model with control function
***************************************************************

* Regress Charc on the mean of l_camIV (CRE), weather controls and their means, and the residual (vhatCRE)
* The residual serves as a control function term to remove the endogenous part of l_FAW
* This gives a consistent estimate of the causal effect under CRE assumptions

ivreghdfe Charc  ///
         rainfall mean_rainfall ///
         sqrainfall mean_sqrainfall temperature mean_temperature ///
         sqtemp mean_sqtemp (l_FAW = l_camIV) ,first
estimates store charcIV		 
		 
/*reghdfe Charc l_FAW ///
         rainfall mean_rainfall ///
         sqrainfall mean_sqrainfall temperature mean_temperature ///
         sqtemp mean_sqtemp vhatCRE		 
outreg2 using results, excel append
estimates store charcIV*/




***************************************************************
* SECTION 17: Lagged FAW Incidence and Impact on Quantity of Charcoal Produced
* Objective: Estimate the impact of previous year's FAW (Fall Armyworm) exposure 
*            on household-level charcoal production quantities. 
* Methodology: Use Tobit models due to the presence of left-censoring in charcoal output 
*              and address endogeneity via a control function approach.
***************************************************************

* Step 1: Convert Qcharc to quantity in kilograms for interpretation consistency
gen qchar = Qcharc * 50
label variable qchar "Charcoal production quantity (kg)"

* Step 2: Log-transform Qcharc (+1 is often used, but here direct log is applied)
* This transformation helps normalize the right-skewed distribution
gen log_Qcharc1 = log(Qcharc + 1)

***************************************************************
* Step 3: Tobit Regression 
* Purpose: Estimate the effect of lagged FAW (via camIV) on charcoal production
* Model: Tobit handles left-censoring (many zeros in Qcharc)
***************************************************************
tobit log_Qcharc army_aff rainfall sqrainfall temperature sqtemp, ll(0) 
margins, dydx(*) predict(ystar(0,.)) post
estimates store QcharcFE 


tobit log_Qcharc l_camIV rainfall sqrainfall temperature sqtemp, ll(0) 
margins, dydx(*) predict(ystar(0,.)) post
estimates store QcharcITT
 
***************************************************************
* Step 4: First-stage regression for control function approach

***************************************************************
xtreg army_aff camIV rainfall sqrainfall temperature sqtemp, fe

* Step 5: Predict residuals (vhat) from first stage
* These residuals capture the endogenous variation in army_aff not explained by instrument
predict vhatcharc, resid

***************************************************************
* Step 6: Tobit Regression with Control Function Residual (Second Stage)
* Purpose: Include residual (vhatcharc) in Tobit model to control for endogeneity


tobit log_Qcharc l_FAW rainfall sqrainfall temperature sqtemp vhatcharc, ll(0) 
margins, dydx(*) predict(ystar(0,.)) post
estimates store QcharcIV





***Coefficient plot***

coefplot (yieldITT, asequation(Yield (ITT)) \ yieldIV, asequation(Yield (LATE)) \FCS_ITT, asequation(FCS (ITT)) \ FCS_IV, asequation(FCS (LATE))\ HDDS_ITT, asequation(HDDS (ITT)) \ HDDS_IV, asequation(HDDS (LATE))\ CSI_ITT, asequation(rCSI (ITT)) \ CSI_IV, asequation(rCSI (LATE))\charcITT, asequation(Charcoal 1 = Yes (ITT))\ charcIV, asequation(Charcoal 1 = Yes (LATE))\ QcharcITT, asequation(Qcharcoal (ITT))\ QcharcIV, asequation(Qcharcoal (LATE))\Cult_ITT, asequation(Cultivated land (ITT))\ Cult_IV, asequation(Cultivated land (LATE))\) ///
, drop(_cons) keep(army_aff camIV l_FAW l_camIV mean_l_FAW) ///
mlabel(cond(@pval<.01, "***", cond(@pval<.05, "**", cond(@pval<.1, "*", "")))) mlabc(red) mlabsize(medium) mlabgap(6pt) mlabp(9) /// 
eqlabel(, labsize(small)) msymbol(circle) mcolor(ebblue) ciopts(recast(rcap) lcolor(ebblue)) xline(0,lcolor(red)) ylabel(, labsize(vsmall)) omitted baselevels /// 
bylabel( "Yield") ///
byopts(compact rows(1) note("p-values shown alongside markers" "*** p<.01, ** p<.05, * p<.1" )) ///
name (debt_fl, replace) ysize (30) xsize(50) 




************************************************************
* Test of Leads for Armyworm (FAW) Exposure
* Objective: Evaluate the validity of instruments and the causal timing
* Purpose:
*   - Use a one-period lead of the treatment (FAW exposure) to test for 
*     anticipatory effects or reverse causality.
*   - A significant lead coefficient would suggest that treatment affects
*     outcomes *before* it occurs—implying a violation of causal assumptions.
************************************************************

************************************************************
* STEP 1: Generate One-Period Leads of Key Endogenous Variables
************************************************************

* Generate a one-period lead of actual FAW exposure (army_aff)
* This is used as a placebo variable to test whether future FAW predicts current outcomes
gen lead = F1.army_aff
label variable lead "Lead of army_aff (used as placebo/instrument)"

* Generate a one-period lead of the instrument (camIV)
* This tests for exogeneity of the instrument—future instrument values 
* should not predict current outcomes if instrument is valid
gen lead_camIV = F1.camIV
label variable lead_camIV "Lead of camIV (used as placebo/instrument)"

************************************************************
* STEP 2: Regress Outcome (Yield) on Lead and Current Instruments
************************************************************
reghdfe logyield army_aff lead_camIV lead qseed fert cultivown_land ///
        rainfall sqrainfall temperature sqtemp, ///
        absorb(i.HHID i.year) cl(camp1)
		
* Baseline Fixed Effects Regression (Placebo Test)

reghdfe logyield camIV lead_camIV lead qseed fert cultivown_land ///
        rainfall sqrainfall temperature sqtemp, ///
        absorb(i.HHID i.year) cl(camp1)

outreg2 using results, excel replace 
estimates store yieldITTlead

************************************************************
* STEP 3: Instrumental Variable Regression with Leads
************************************************************
ivreghdfe logyield  qseed fert cultivown_land rainfall ///
          sqrainfall temperature sqtemp ///
          (army_aff lead = camIV lead_camIV), ///
          absorb(i.HHID i.year) cl(camp1) first

* Save and append results to the same Excel output
outreg2 using results, excel append
estimates store yieldIVlead


***** Test of Leads for Armyworm (FAW) Exposure****

xttobit hpsfcs army_aff lead rainfall sqrainfall temperature sqtemp, ll(0) re
margins, dydx(*) predict(ystar(0,.))
outreg2 using results, excel replace
estimates store FCS_FE_lead  

xttobit hpsfcs camIV lead_camIV rainfall sqrainfall temperature sqtemp, ll(0) re
margins, dydx(*) predict(ystar(0,.))
outreg2 using results, excel replace
estimates store FCS_ITT_lead 

* FIRST-STAGE MODEL: Estimate determinants of actual FAW exposure

xtreg army_aff lead camIV lead_camIV rainfall sqrainfall temperature sqtemp, fe

* Generate residuals (vhat) from the first-stage regression
* These residuals represent the endogenous variation in army_aff
predict vhatlead, resid	

* IV MODEL (CONTROL FUNCTION APPROACH): Second-stage Tobit regression

xttobit hpsfcs army_aff lead rainfall sqrainfall temperature sqtemp vhatlead, ll(0) re
margins, dydx(*) predict(ystar(0,.))
outreg2 using results, excel append
estimates store FCS_IV_lead

**HHDS**
xttobit hpshdds army_aff lead_camIV rainfall sqrainfall temperature sqtemp, ll(0) re
margins, dydx(*) predict(ystar(0,.))
outreg2 using results, excel replace
estimates store HDDS_FE_lead  

xttobit hpshdds camIV lead_camIV rainfall sqrainfall temperature sqtemp, ll(0) re
margins, dydx(*) predict(ystar(0,.))
outreg2 using results, excel replace
estimates store HDDS_ITT_lead 

* IV MODEL (CONTROL FUNCTION APPROACH): Second-stage Tobit regression

xttobit hpshdds army_aff lead rainfall sqrainfall temperature sqtemp vhatlead, ll(0) re
margins, dydx(*) predict(ystar(0,.))
outreg2 using results, excel append
estimates store FCS_IV_lead

**CSI**
xttobit hps_rcsi army_aff lead_camIV rainfall sqrainfall temperature sqtemp, ll(0) re
margins, dydx(*) predict(ystar(0,.))
outreg2 using results, excel replace
estimates store CSI_FE_lead  

xttobit hps_rcsi camIV lead_camIV rainfall sqrainfall temperature sqtemp, ll(0) re
margins, dydx(*) predict(ystar(0,.))
outreg2 using results, excel replace
estimates store CSI_ITT_lead 

* IV MODEL (CONTROL FUNCTION APPROACH): Second-stage Tobit regression

xttobit hps_rcsi army_aff lead rainfall sqrainfall temperature sqtemp vhatlead, ll(0) re
margins, dydx(*) predict(ystar(0,.))
outreg2 using results, excel append
estimates store FCS_IV_lead


***Test of Leads – Impact of Future FAW Exposure****

***************************************************************
* 1. Reghdfe with Fixed Effects (Intent-to-Treat with Lead Terms)
***************************************************************
reghdfe forwardcultiv army_aff lead rainfall sqrainfall temperature sqtemp, ///
       absorb(i.year i.HHID) cl(camp1)
outreg2 using results, excel replace
estimates store Cult_FE

reghdfe forwardcultiv camIV lead_camIV rainfall sqrainfall temperature sqtemp, ///
       absorb(i.year i.HHID) cl(camp1)
outreg2 using results, excel replace
estimates store Cult_ITT_lead

***************************************************************
* 2. Instrumental Variables Regression (Two-Stage with Lead Test)
***************************************************************
ivreghdfe forwardcultiv rainfall sqrainfall temperature sqtemp ///
          (army_aff lead = camIV lead_camIV), ///
          absorb(i.year i.HHID) cl(camp1) first
outreg2 using results, excel append
estimates store Cult_IV_lead


***Test of leads***

quietly regress Charc l_FAW lead  rainfall ///
                 sqrainfall temperature sqtemp

* Retain only estimation sample for downstream operations (to align samples)
keep if e(sample)

foreach v in l_FAW lead rainfall ///
            sqrainfall temperature sqtemp {
    bysort HHID: egen double mean_`v' = mean(`v')
}

mfx 

estimates store charcFE_lead
outreg2 using results, excel replace


quietly regress Charc l_camIV lead_camIV  rainfall ///
                 sqrainfall temperature sqtemp

* Retain only estimation sample for downstream operations (to align samples)
keep if e(sample)

foreach v in l_camIV lead_camIV rainfall ///
            sqrainfall temperature sqtemp {
    bysort HHID: egen double mean_`v' = mean(`v')
}

mfx 

estimates store charcFE_lead
outreg2 using results, excel replace


reghdfe Charc l_FAW lead ///
        rainfall mean_rainfall ///
        sqrainfall mean_sqrainfall temperature mean_temperature ///
        sqtemp mean_sqtemp,resid

* Predict residuals from first-stage regression — these will capture endogeneity
* Inclusion of vhatCRE in second-stage regression adjusts for endogenous selection into FAW exposure	
predict vhatCRElead, resid

reghdfe Charc l_FAW lead ///
         rainfall mean_rainfall ///
         sqrainfall mean_sqrainfall temperature mean_temperature ///
         sqtemp mean_sqtemp vhatCRElead		 
outreg2 using results, excel append
estimates store charcIV_lead 



***Test of leads***
tobit log_Qcharc army_aff lead  rainfall sqrainfall temperature sqtemp, ll(0) 
* Step 3: Post-estimation marginal effects
margins, dydx(*) predict(ystar(0,.))
outreg2 using results, excel append
estimates store QcharcFE_leads

tobit log_Qcharc l_camIV lead_camIV  rainfall sqrainfall temperature sqtemp, ll(0) 
* Step 3: Post-estimation marginal effects
margins, dydx(*) predict(ystar(0,.))
outreg2 using results, excel append
estimates store QcharcITT_leads

xtreg army_aff lead camIV lead_camIV rainfall sqrainfall temperature sqtemp, fe

* Step 5: Predict residuals (vhat) from first stage
* These residuals capture the endogenous variation in army_aff not explained by instrument
predict vhatcharclead, resid 

tobit log_Qcharc l_FAW lead rainfall sqrainfall temperature sqtemp vhatcharclead, ll(0) 
* Step 6: Marginal effects for control function Tobit model
* This shows effect of each variable on latent outcome y*
margins, dydx(*) predict(ystar(0,.))
outreg2 using results, excel append
estimates store QcharcIV_lead




************************************************************
* SECTION 18: Heterogeneous Effects of FAW 
* Subtopics: Distance to Tarmac and Asset-Based Wealth
* Purpose: Explore whether the effects of FAW on cultivated land 
*          and food security outcomes differ by household wealth 
*          (asset-based) and biomass availability (proxy for fuel/environmental endowment).
************************************************************ 


************************************************************
* I. WEALTH HETEROGENEITY: Asset-Based Interactions
************************************************************

**************************************************************
* STEP 1: Construct Asset Quantiles Based on 2016 Values
**************************************************************

* Retain asset index for 2016 to ensure consistency in wealth comparison across years
gen assetind_ilri16 = assetind_ilri if year == 2016

* Forward-fill 2016 asset values to subsequent years (panel format)
* This allows stratification by baseline wealth while estimating FAW impacts post-2016
bysort HHID (year): replace assetind_ilri16 = assetind_ilri16[_n-1] if missing(assetind_ilri16)

* Keep filled-in values for years of interest (2016–2019)
gen assetind_ilri16_filled = assetind_ilri16 if inrange(year, 2016, 2019)

* Create asset index quartiles to stratify households into wealth groups
xtile assetq_ = assetind_ilri16_filled, nq(4)
label variable assetq_ "Quantile of Asset Index (2016)"

* Create dummy variables for each quartile
tabulate assetq_, generate(qassetq_)

***Log assets***
gen log_assets  = ln(assetind_ilri16)

bysort HHID (year): replace log_assets = log_assets[_n-1] if missing(log_assets)

gen log_assets_filled = log_assets if inrange(year, 2016, 2019)

xtile logassetq_ = log_assets_filled, nq(4)
label variable logassetq_ "Quantile of Asset Index (2016)"

tabulate logassetq_, generate(logassetq)

gen faw_logassetq1 = logassetq1 * army_aff
gen faw_logassetq2 = logassetq2 * army_aff
gen faw_logassetq3 = logassetq3 * army_aff
gen faw_logassetq4 = logassetq4 * army_aff

gen IVfaw_logassetq1 = logassetq1 * camIV
gen IVfaw_logassetq2 = logassetq2 * camIV
gen IVfaw_logassetq3 = logassetq3 * camIV
gen IVfaw_logassetq4 = logassetq4 * camIV

reghdfe Charc IVfaw_logassetq1 IVfaw_logassetq2 IVfaw_logassetq3 IVfaw_logassetq4 ///
                 rainfall sqrainfall temperature sqtemp, ///
                 absorb(i.year) cluster(camp1)



**************************************************************
* STEP 2: Interact FAW Exposure with Asset Quartiles
**************************************************************

* Interaction terms to estimate conditional effect of FAW by wealth quantile
gen faw_assetq1 = qassetq_1 * army_aff
gen faw_assetq2 = qassetq_2 * army_aff
gen faw_assetq3 = qassetq_3 * army_aff
gen faw_assetq4 = qassetq_4 * army_aff

gen IVfaw_assetq1 = qassetq_1 * camIV
gen IVfaw_assetq2 = qassetq_2 * camIV
gen IVfaw_assetq3 = qassetq_3 * camIV
gen IVfaw_assetq4 = qassetq_4 * camIV

**************************************************************
* STEP 3: Estimate FAW Effect on Cultivated Land by Wealth Group
**************************************************************

 ivreghdfe forwardcultiv ///
    (faw_assetq1 faw_assetq2 faw_assetq3 faw_assetq4 = ///
     IVfaw_logassetq1 IVfaw_logassetq2 IVfaw_logassetq3 IVfaw_logassetq4) ///
    rainfall sqrainfall temperature sqtemp, ///
    absorb(i.year) cluster(camp1)
estimates store hete_cult2

reghdfe forwardcultiv IVfaw_logassetq1 IVfaw_logassetq2 IVfaw_logassetq3 IVfaw_logassetq4 ///
        rainfall sqrainfall temperature sqtemp, ///
        absorb(i.year) cluster(camp1)

ivreghdfe forwardcultiv (faw_logassetq1 faw_logassetq2 faw_logassetq3 faw_logassetq4 = ///
                 IVfaw_logassetq1 IVfaw_logassetq2 IVfaw_logassetq3 IVfaw_logassetq4) ///
                 rainfall sqrainfall temperature sqtemp, ///
                 absorb(i.year) cluster(camp1)
estimates store hete_cult2

***Assets results***
reghdfe forwardcultiv faw_assetq1 faw_assetq2 faw_assetq3 faw_assetq4 ///
                 rainfall sqrainfall temperature sqtemp, ///
                 absorb(i.HHID i.year) cluster(camp1)
estimates store hete_cult3_FE

reghdfe forwardcultiv IVfaw_assetq1 IVfaw_assetq2 IVfaw_assetq3 IVfaw_assetq4 ///
                 rainfall sqrainfall temperature sqtemp, ///
                 absorb(i.HHID i.year) cluster(camp1)
estimates store hete_cult3_ITT


ivreghdfe	forwardcultiv (faw_assetq1 faw_assetq2 faw_assetq3 faw_assetq4 = ///
		IVfaw_assetq1 IVfaw_assetq2 IVfaw_assetq3 IVfaw_assetq4) ///
		rainfall sqrainfall temperature sqtemp, ///
		absorb(i.HHID i.year) cluster(camp1)

estimates store hete_cult3
 
* Generate time means of the endogenous variables and controls
foreach var in faw_assetq1 faw_assetq2 faw_assetq3 faw_assetq4 ///
               rainfall sqrainfall temperature sqtemp {
    egen mean_`var' = mean(`var'), by(HHID)
}

ivreghdfe Charc ///
    faw_assetq1 faw_assetq2 faw_assetq3 faw_assetq4 ///
	mean_faw_assetq1 mean_faw_assetq2 mean_faw_assetq3 mean_faw_assetq4 ///
    mean_rainfall mean_sqrainfall mean_temperature mean_sqtemp ///
    rainfall sqrainfall temperature sqtemp ///
    , absorb(i.year) cluster(camp1)
	

reghdfe Charc ///
     IVfaw_assetq1 IVfaw_assetq2 IVfaw_assetq3 IVfaw_assetq4 ///
    mean_rainfall mean_sqrainfall mean_temperature mean_sqtemp ///
    rainfall sqrainfall temperature sqtemp ///
    , absorb(i.year) cluster(camp1)
	
	
 
ivreghdfe Charc ///
    (faw_assetq1 faw_assetq2 faw_assetq3 faw_assetq4 = ///
     IVfaw_assetq1 IVfaw_assetq2 IVfaw_assetq3 IVfaw_assetq4) ///
    mean_rainfall mean_sqrainfall mean_temperature mean_sqtemp ///
    rainfall sqrainfall temperature sqtemp ///
    , absorb(i.year) cluster(camp1)
	
estimates store hete_charc3


reghdfe forwardcultiv faw_assetq1 faw_assetq2 faw_assetq3 faw_assetq4 ///
        rainfall sqrainfall temperature sqtemp, ///
        absorb(i.year) cluster(camp1)
outreg2 using results, excel replace
estimates store hete_cult

**************************************************************
* STEP 4: Biomass-Based Heterogeneity (Environmental Constraint)
**************************************************************
gen base_treecov = treecov if year == 2016

bysort HHID (year): replace base_treecov = base_treecov[_n-1] if missing(base_treecov)

gen base_treecov_filled = base_treecov if inrange(year, 2016, 2019)	

* Create biomass quartiles
xtile b_treecov = base_treecov_filled, nq(4)
label variable b_treecov "Quantile of Biomass (2016)"

* Generate dummies for biomass quantiles
tabulate b_treecov, generate(b_treecovq)

* Interact biomass quartile dummies with FAW exposure
gen faw_b_treecov1 = b_treecovq1 * army_aff
gen faw_b_treecov2 = b_treecovq2 * army_aff
gen faw_b_treecov3 = b_treecovq3 * army_aff
gen faw_b_treecov4 = b_treecovq4 * army_aff	

gen biomassIV1 = b_treecovq1 * camIV
gen biomassIV2 = b_treecovq2 * camIV
gen biomassIV3 = b_treecovq3 * camIV
gen biomassIV4 = b_treecovq4 * camIV

***Biomass results**
 ivreghdfe forwardcultiv ///
    faw_b_treecov1 faw_b_treecov2 faw_b_treecov3 faw_b_treecov4 ///
    rainfall sqrainfall temperature sqtemp, ///
    absorb(i.HHID i.year) cluster(camp1)	
estimates store bio_cultFE


 ivreghdfe forwardcultiv ///
    biomassIV1 biomassIV2 biomassIV3 biomassIV4 ///
    rainfall sqrainfall temperature sqtemp, ///
    absorb(i.HHID i.year) cluster(camp1)	
estimates store bio_cultITT

 ivreghdfe forwardcultiv ///
    (faw_b_treecov1 faw_b_treecov2 faw_b_treecov3 faw_b_treecov4 = ///
    biomassIV1 biomassIV2 biomassIV3 biomassIV4) ///
    rainfall sqrainfall temperature sqtemp, ///
    absorb(i.HHID i.year) cluster(camp1)	
estimates store bio_cult3


 ivreghdfe lncultivation ///
    (faw_b_treecov1 faw_b_treecov2 faw_b_treecov3 faw_b_treecov4 = ///
    biomassIV1 biomassIV2 biomassIV3 biomassIV4) ///
    rainfall sqrainfall temperature sqtemp, ///
    absorb(i.year) cluster(camp1)	
estimates store bio_cult4


foreach var in faw_b_treecov1 faw_b_treecov2 faw_b_treecov3 faw_b_treecov4 ///
               {
    egen mean_`var' = mean(`var'), by(HHID)
}

 ivreghdfe Charc ///
    faw_b_treecov1 faw_b_treecov2 faw_b_treecov3 faw_b_treecov4 ///
	    mean_rainfall mean_sqrainfall mean_temperature mean_sqtemp ///
    rainfall sqrainfall temperature sqtemp, ///
    absorb(i.year) cluster(camp1)


 ivreghdfe Charc ///
    biomassIV1 biomassIV2 biomassIV3 biomassIV4 ///
	    mean_rainfall mean_sqrainfall mean_temperature mean_sqtemp ///
    rainfall sqrainfall temperature sqtemp, ///
    absorb(i.year) cluster(camp1)

 ivreghdfe Charc ///
    (faw_b_treecov1 faw_b_treecov2 faw_b_treecov3 faw_b_treecov4 = ///
     biomassIV1 biomassIV2 biomassIV3 biomassIV4) ///
    mean_rainfall mean_sqrainfall mean_temperature mean_sqtemp ///
    rainfall sqrainfall temperature sqtemp ///
    , absorb(i.year) cluster(camp1)
estimates store bio_charc3		
	
	
* Construct Potential yield Quantiles Based on 2016 Values
**************************************************************
*egen mean_yield_camp_year = mean(yield_maize), by(camp year)
gen yield_maize16 = yield_maize if year == 2016

* Forward-fill 2016 asset values to subsequent years (panel format)
* This allows stratification by baseline wealth while estimating FAW impacts post-2016
bysort HHID (year): replace yield_maize16 = yield_maize16[_n-1] if missing(yield_maize16)

* Keep filled-in values for years of interest (2016–2019)
gen yield_maize16_filled = yield_maize16 if inrange(year, 2016, 2019)

* Create asset index quartiles to stratify households into wealth groups
xtile yieldq = yield_maize16_filled, nq(4)
label variable yieldq "Quantile of Yield (2016)"

* Create dummy variables for each quartile
tabulate yieldq, generate(yieldq_)

gen faw_yieldq1 = yieldq_1 * army_aff
gen faw_yieldq2 = yieldq_2 * army_aff
gen faw_yieldq3 = yieldq_3 * army_aff
gen faw_yieldq4 = yieldq_4 * army_aff

gen IVfaw_yieldq1 = faw_yieldq1 * camIV
gen IVfaw_yieldq2 = faw_yieldq2 * camIV
gen IVfaw_yieldq3 = faw_yieldq3 * camIV
gen IVfaw_yieldq4 = faw_yieldq4 * camIV


***Distance results***

 ivreghdfe forwardcultiv ///
    faw_yieldq1 faw_yieldq2 faw_yieldq3 faw_yieldq4 
    rainfall sqrainfall temperature sqtemp, ///
    absorb( i.year) cluster(camp1)
estimates store yield_cultFE

 ivreghdfe forwardcultiv ///
    (faw_yieldq1 faw_yieldq2 faw_yieldq3 faw_yieldq4 = ///
     IVfaw_yieldq1 IVfaw_yieldq2 IVfaw_yieldq3 IVfaw_yieldq4) ///
    rainfall sqrainfall temperature sqtemp, ///
    absorb( i.year) cluster(camp1)
estimates store yield_cult2



 ivreghdfe lncultivation ///
    (faw_yieldq1 faw_yieldq2 faw_yieldq3 faw_yieldq4 = ///
     IVfaw_yieldq1 IVfaw_yieldq2 IVfaw_yieldq3 IVfaw_yieldq4) ///
    rainfall sqrainfall temperature sqtemp, ///
    absorb(i.year) cluster(camp1)
estimates store yield_cult3



foreach var in faw_yieldq1 faw_yieldq2 faw_yieldq3 faw_yieldq4  ///
                {
    egen mean_`var' = mean(`var'), by(HHID)
}

 ivreghdfe Charc ///
    (faw_yieldq1 faw_yieldq2 faw_yieldq3 faw_yieldq4 = ///
     IVfaw_yieldq1 IVfaw_yieldq2 IVfaw_yieldq3 IVfaw_yieldq4) ///
    rainfall sqrainfall temperature sqtemp ///
    mean_faw_yieldq1 mean_faw_yieldq2 mean_faw_yieldq1 mean_faw_yieldq1 ///
    mean_rainfall mean_sqrainfall mean_temperature mean_sqtemp, ///
    absorb(i.year) cluster(camp1)
estimates store yield_charc2

	

*****************************************************
* CLEANING AND PREPARING DISTANCE TO TARMAC VARIABLE
* Purpose: Create clean, normalized, and interaction-ready
*          version of the distance to tarmac variable
*****************************************************

* Step 0: Sort dataset in descending year order per household to support backfilling from latest year
gsort HHID -year

* Step 1: Convert distance-to-tarmac variable to numeric format
* 'force' is used here to ensure non-numeric strings are treated as missing
destring tarmac_dist, generate(dist_tarmac) force

* Step 2: Check for duplicated household-year observations to ensure panel integrity
duplicates report HHID year

* Step 3: Replace common non-system missing value codes (e.g., 99, -999) with system missing
gen dist_tarmac_filled = dist_tarmac
replace dist_tarmac_filled = . if inlist(dist_tarmac_filled, 99, 999, -99, -999)

* Step 4: Fill missing 2019 values using the camp-level mean (camp as geographic cluster)
egen camp_mean_2019 = mean(dist_tarmac_filled) if year == 2019, by(camp)
replace dist_tarmac_filled = camp_mean_2019 if missing(dist_tarmac_filled) & year == 2019

*************************************************************
* BACKFILL 2019 DISTANCE VALUES TO EARLIER YEARS PER HHID
* Logic: Distance is assumed time-invariant or updated less frequently
*************************************************************

* Step 5: Generate variable for 2019 values only
gen tarmac_2019 = dist_tarmac_filled if year == 2019

* Step 6: Backfill 2019 value across earlier years for same household
* _n+2 used under the assumption of three-year panel (e.g., 2017–2019)

bysort HHID (year): replace tarmac_2019 = tarmac_2019[_n+1] if missing(tarmac_2019)
bysort HHID (year): replace tarmac_2019 = tarmac_2019[_n+2] if missing(tarmac_2019)



* Step 7: Log-transform distance to reduce skewness and normalize
gen logdist = ln(tarmac_2019)

************************************************************
* CREATE DISTANCE QUANTILES BASED ON LOG-TRANSFORMED VALUES
************************************************************

* Step 8: Create quartile indicators for distance to tarmac
xtile mean_dist_tarmac_q = logdist, nq(4)

* Step 9: Generate dummy variables for each quantile (for interaction terms)
tabulate mean_dist_tarmac_q, generate(qtq_)

************************************************************
* Step 10–11: Generate Interaction Terms for Heterogeneous Effects
************************************************************

* Step 11: Interact FAW variable with each distance quartile
* These will allow estimating heterogeneous effects of FAW by remoteness
gen faw_q1 = qtq_1 * army_aff
gen faw_q2 = qtq_2 * army_aff
gen faw_q3 = qtq_3 * army_aff
gen faw_q4 = qtq_4 * army_aff


gen faw_q1_IV = qtq_1 * camIV
gen faw_q2_IV = qtq_2 * camIV
gen faw_q3_IV = qtq_3 * camIV
gen faw_q4_IV = qtq_4 * camIV


***************************************
* MAIN REGRESSIONS AND MARGINAL EFFECTS
***************************************

* Step 13: Estimate the effect of FAW*distance interaction on cultivated land
* Absorbing household and year fixed effects; clustering at camp level
reghdfe forwardcultiv faw_q1 faw_q2 faw_q3 faw_q4 rainfall ///
        sqrainfall temperature sqtemp, absorb(i.HHID i.year) cl(camp1)
estimates store dist_cult_FE


reghdfe forwardcultiv faw_q1_IV faw_q2_IV faw_q3_IV faw_q4_IV rainfall ///
        sqrainfall temperature sqtemp, absorb(i.HHID i.year) cl(camp1)
estimates store dist_cult_ITT

 ivreghdfe forwardcultiv ///
    (faw_q1 faw_q2 faw_q3 faw_q4 = ///
    faw_q1_IV faw_q2_IV faw_q3_IV faw_q4_IV) ///
    rainfall sqrainfall temperature sqtemp, ///
    absorb(i.HHID i.year) cluster(camp1)
estimates store dist_cult3
	
 ivreghdfe lncultivation ///
    (faw_q1 faw_q2 faw_q3 faw_q4 = ///
    faw_q1_IV faw_q2_IV faw_q3_IV faw_q4_IV) ///
    rainfall sqrainfall temperature sqtemp, ///
    absorb( i.year) cluster(camp1)
estimates store dist_cult4

foreach var in faw_q1 faw_q2 faw_q3 faw_q4 ///
               {
    egen mean_`var' = mean(`var'), by(HHID)
}


 ivreghdfe Charc ///
        faw_q1 faw_q2 faw_q3 faw_q4  ///
    rainfall sqrainfall temperature sqtemp ///
    mean_rainfall mean_sqrainfall mean_temperature mean_sqtemp, ///
    absorb(i.year) cluster(camp1)
	

 ivreghdfe Charc ///
        faw_q1_IV faw_q2_IV faw_q3_IV faw_q4_IV  ///
    rainfall sqrainfall temperature sqtemp ///
    mean_rainfall mean_sqrainfall mean_temperature mean_sqtemp, ///
    absorb(i.year) cluster(camp1)	

 ivreghdfe Charc ///
        (faw_q1 faw_q2 faw_q3 faw_q4 = ///
    faw_q1_IV faw_q2_IV faw_q3_IV faw_q4_IV) ///
    rainfall sqrainfall temperature sqtemp ///
    mean_rainfall mean_sqrainfall mean_temperature mean_sqtemp, ///
    absorb(i.year) cluster(camp1)
estimates store dist_charc3	


 ivreghdfe forwardcultiv ///
    (faw_q1 faw_q2 faw_q3 faw_q4 = ///
    faw_q1_IV faw_q2_IV faw_q3_IV faw_q4_IV) ///
    rainfall sqrainfall temperature sqtemp, ///
    absorb(i.year) cluster(camp1)
estimates store dist_cult3
	


 ivreghdfe Charc ///
    (faw_q1 faw_q2 faw_q3 faw_q4 = ///
     faw_q1_IV faw_q2_IV faw_q3_IV faw_q4_IV) ///
    rainfall sqrainfall temperature sqtemp ///
    ,absorb(i.year) cluster(camp1)
estimates store dist_charc3



label var faw_assetq1 "1"
label var faw_assetq2 "2"
label var faw_assetq3 "3"
label var faw_assetq4 "4"


label var faw_b_treecov1 "1"
label var faw_b_treecov2 "2"
label var faw_b_treecov3 "3"
label var faw_b_treecov4 "4"

label var faw_q1 "1"
label var faw_q2 "2"
label var faw_q3 "3"
label var faw_q4 "4"

label var faw_yieldq1 "1"
label var faw_yieldq2 "2"
label var faw_yieldq3 "3"
label var faw_yieldq4 "4"





coefplot ///
    (hete_cult3, asequation("Assets")) ///
    || (bio_cult3, asequation("Biomass")) ///
    || (dist_cult3, asequation("Distance")) ///
    || (hete_charc3, asequation("Assets")) ///
    || (bio_charc3, asequation("Biomass")) ///
    || (dist_charc3, asequation("Distance")) ///
    , vertical ///
    drop(_cons) ///
    keep(faw_assetq1 faw_assetq2 faw_assetq3 faw_assetq4 ///
         faw_b_treecov1 faw_b_treecov2 faw_b_treecov3 faw_b_treecov4 ///
         faw_q1 faw_q2 faw_q3 faw_q4 ///
         faw_yieldq1 faw_yieldq2 faw_yieldq3 faw_yieldq4) ///
    coeflabels( ///
        faw_assetq1 = "Q1" faw_assetq2 = "Q2" faw_assetq3 = "Q3" faw_assetq4 = "Q4" ///
        faw_b_treecov1 = "Q1" faw_b_treecov2 = "Q2" faw_b_treecov3 = "Q3" faw_b_treecov4 = "Q4" ///
        faw_q1 = "Q1" faw_q2 = "Q2" faw_q3 = "Q3" faw_q4 = "Q4" ///
        faw_yieldq1 = "Q1" faw_yieldq2 = "Q2" faw_yieldq3 = "Q3" faw_yieldq4 = "Q4" ///
    ) ///
    mlabc(red) mlabsize(medium) mlabgap(6pt) mlabp(9) ///
    eqlabel(, labsize(small)) ///
    msymbol(circle) mcolor(ebblue) ///
    ciopts(recast(rarea) color(ebblue%30)) ///
    yline(0, lcolor(gray)) ///
    ylabel(, labsize(vsmall)) ///
    omitted baselevels ///
    bylabels("Assets - Cultivated land" ///
             "Biomass - Cultivated land" ///
             "Distance - Cultivated land" ///
             "Assets - Charcoal" ///
             "Biomass - Charcoal" ///
             "Distance - Charcoal") ///
    byopts(cols(3) graphregion(color(white)) ///
           note("") title("Heterogeneity in Land Use")) ///
    xscale(range(0.5 4.5)) ///
    name(heterogeneity_boxed, replace) ///
    ysize(30) xsize(50)

	
coefplot ///
    (hete_cult3, asequation("Assets")) ///
    || (bio_cult3, asequation("Biomass")) ///
    || (dist_cult3, asequation("Distance")) ///
    || (hete_charc3, asequation("Assets")) ///
    || (bio_charc3, asequation("Biomass")) ///
    || (dist_charc3, asequation("Distance")) ///
    , vertical ///
    drop(_cons) ///
    keep(faw_assetq1 faw_assetq2 faw_assetq3 faw_assetq4 ///
         faw_b_treecov1 faw_b_treecov2 faw_b_treecov3 faw_b_treecov4 ///
         faw_q1 faw_q2 faw_q3 faw_q4 ///
         faw_yieldq1 faw_yieldq2 faw_yieldq3 faw_yieldq4) ///
    coeflabels( ///
        faw_assetq1 = "Q1" faw_assetq2 = "Q2" faw_assetq3 = "Q3" faw_assetq4 = "Q4" ///
        faw_b_treecov1 = "Q1" faw_b_treecov2 = "Q2" faw_b_treecov3 = "Q3" faw_b_treecov4 = "Q4" ///
        faw_q1 = "Q1" faw_q2 = "Q2" faw_q3 = "Q3" faw_q4 = "Q4" ///
        faw_yieldq1 = "Q1" faw_yieldq2 = "Q2" faw_yieldq3 = "Q3" faw_yieldq4 = "Q4" ///
    ) ///
    mlabc(red) mlabsize(medium) mlabgap(6pt) mlabp(9) ///
    eqlabel(, labsize(small)) ///
    msymbol(circle) mcolor(ebblue) ///
    ciopts(recast(rarea) color(ebblue%30)) ///
    yline(0, lcolor(gray)) ///
    ylabel(, labsize(vsmall)) ///
    omitted baselevels ///
    bylabels("Assets - Cultivated land" ///
             "Biomass - Cultivated land" ///
             "Distance - Cultivated land" ///
             "Assets - Charcoal" ///
             "Biomass - Charcoal" ///
             "Distance - Charcoal") ///
    byopts(cols(3) graphregion(color(white)) ///
           xrescale ///
           note("") title("Heterogeneity in Land Use")) ///
    name(heterogeneity_boxed, replace) ///
    ysize(30) xsize(50)
	
	
coefplot ///
    (hete_cult3, asequation("Assets")) ///
    || (bio_cult3, asequation("Biomass")) ///
    || (dist_cult3, asequation("Distance")) ///
    || (hete_charc3, asequation("Assets")) ///
    || (bio_charc3, asequation("Biomass")) ///
    || (dist_charc3, asequation("Distance")) ///
    , vertical ///
    drop(_cons) ///
    keep(faw_assetq1 faw_assetq2 faw_assetq3 faw_assetq4 ///
         faw_b_treecov1 faw_b_treecov2 faw_b_treecov3 faw_b_treecov4 ///
         faw_q1 faw_q2 faw_q3 faw_q4 ///
         faw_yieldq1 faw_yieldq2 faw_yieldq3 faw_yieldq4) ///
    coeflabels( ///
        faw_assetq1 = "Q1" faw_assetq2 = "Q2" faw_assetq3 = "Q3" faw_assetq4 = "Q4" ///
        faw_b_treecov1 = "Q1" faw_b_treecov2 = "Q2" faw_b_treecov3 = "Q3" faw_b_treecov4 = "Q4" ///
        faw_q1 = "Q1" faw_q2 = "Q2" faw_q3 = "Q3" faw_q4 = "Q4" ///
        faw_yieldq1 = "Q1" faw_yieldq2 = "Q2" faw_yieldq3 = "Q3" faw_yieldq4 = "Q4" ///
    ) ///
    mlabc(red) mlabsize(medium) mlabgap(6pt) mlabp(9) ///
    eqlabel(, labsize(small)) ///
    msymbol(circle) mcolor(ebblue) ///
    ciopts(recast(rarea) color(ebblue%30)) ///
    yline(0, lcolor(gray)) ///
    ylabel(, labsize(vsmall)) ///
    omitted baselevels ///
    bylabels("Assets - Cultivated land" ///
             "Biomass - Cultivated land" ///
             "Distance - Cultivated land" ///
             "Assets - Charcoal" ///
             "Biomass - Charcoal" ///
             "Distance - Charcoal") ///
    byopts(cols(3) graphregion(color(white)) ///
           xrescale ///
           note("") title("Heterogeneity in Land Use")) ///
    xlabel(1 "Q1" 2 "Q2" 3 "Q3" 4 "Q4", labsize(vsmall)) ///
    xscale(range(1 4.5) noextend) ///
    name(heterogeneity_boxed, replace) ///
    ysize(30) xsize(50)
	
	
coefplot ///
    || (dist_cult3, asequation("Distance")) ///
    || (dist_charc3, asequation("Distance")) ///
    , vertical ///
    drop(_cons) ///
    keep(faw_assetq1 faw_assetq2 faw_assetq3 faw_assetq4 ///
         faw_q1 faw_q2 faw_q3 faw_q4) ///
    coeflabels( ///
        faw_q1 = "Q1" faw_q2 = "Q2" faw_q3 = "Q3" faw_q4 = "Q4") ///
    mlabc(red) mlabsize(medium) mlabgap(6pt) mlabp(9) ///
    eqlabel(, labsize(small)) ///
    msymbol(circle) mcolor(ebblue) ///
    ciopts(recast(rarea) color(ebblue%30)) ///
    yline(0, lcolor(gray)) ///
    ylabel(, labsize(vsmall)) ///
    omitted baselevels ///
    xlabel(1 "Q1" 2 "Q2" 3 "Q3" 4 "Q4", labsize(vsmall)) ///
    xscale(range(0.5 4.5) noextend) ///
    bylabels("Distance - Cultivated land") ///
    byopts(cols(2) xrescale ///
           graphregion(color(white) margin(t=1 b=1 l=1 r=1)) ///
           note("") title("Heterogeneity in Land Use")) ///
    name(biomass_distance_boxed, replace) ///
    ysize(25) xsize(40)
	
	

	
coefplot ///
    (hete_cult3, asequation("Assets")) ///
    || (bio_cult3, asequation("Biomass")) ///
    || (dist_cult3, asequation("Distance")) ///
    || (hete_charc3, asequation("Assets")) ///
    || (bio_charc3, asequation("Biomass")) ///
    || (dist_charc3, asequation("Distance")) ///
    , vertical ///
    drop(_cons) ///
    keep(faw_assetq1 faw_assetq2 faw_assetq3 faw_assetq4 ///
         faw_b_treecov1 faw_b_treecov2 faw_b_treecov3 faw_b_treecov4 ///
         faw_q1 faw_q2 faw_q3 faw_q4 ///
         faw_yieldq1 faw_yieldq2 faw_yieldq3 faw_yieldq4) ///
    coeflabels( ///
        faw_assetq1 = "Q1" faw_assetq2 = "Q2" faw_assetq3 = "Q3" faw_assetq4 = "Q4" ///
        faw_b_treecov1 = "Q1" faw_b_treecov2 = "Q2" faw_b_treecov3 = "Q3" faw_b_treecov4 = "Q4" ///
        faw_q1 = "Q1" faw_q2 = "Q2" faw_q3 = "Q3" faw_q4 = "Q4" ///
        faw_yieldq1 = "Q1" faw_yieldq2 = "Q2" faw_yieldq3 = "Q3" faw_yieldq4 = "Q4" ///
    ) ///
    mlabc(red) mlabsize(medium) mlabgap(6pt) mlabp(9) ///
    eqlabel(, labsize(small)) ///
    msymbol(circle) mcolor(ebblue) ///
    ciopts(recast(rarea) color(ebblue%30)) ///
    yline(0, lcolor(gray)) ///
    ylabel(, labsize(vsmall)) ///
    omitted baselevels ///
    bylabels("Assets - Cultivated land" ///
             "Biomass - Cultivated land" ///
             "Distance - Cultivated land" ///
             "Assets - Charcoal" ///
             "Biomass - Charcoal" ///
             "Distance - Charcoal") ///
    byopts(cols(3) graphregion(color(white)) ///
           note("") title("Heterogeneity in Land Use") ///
           xrescale) ///
    xlabel(1 "Q1" 2 "Q2" 3 "Q3" 4 "Q4", labsize(vsmall)) ///
    xscale(range(0.5 4.5)) ///
    name(heterogeneity_boxed, replace) ///
    ysize(30) xsize(50)
	



