# Data

Source: International Association for the Evaluation of Educational Achievement (IEA), *ICILS 2018 International Database*, [doi:10.58150/ICILS_2018_data](https://doi.org/10.58150/ICILS_2018_data). Use is subject to the [IEA terms and conditions](https://www.iea.nl/data-tools/repository/icils). Source data remain the property of IEA.

`raw/` contains student (`BSG`) and school (`BCG`) SPSS extracts for Germany, Denmark, Finland, France, Korea, Luxembourg, Portugal, and the United States. Only variables used in the analyses are included. All original rows, missing values, public-use identifiers, and variable labels are retained. These files are the direct inputs to the analysis scripts.

`analysis_data.csv` is a merged, compressed CSV of 29,832 student records, including incomplete records. School variables are joined by `IDSCHOOL` within country. Blank values denote missing data. `COUNTRY_CODE` identifies the national sample; `SCHOOL_UID` and `ROW_UID` combine the country code with the public-use school and student identifiers. `complete_case` indicates observed values for all 25 predictors and five CT plausible values. 

`variable_dictionary.csv` describes the 25 predictors and five CT outcomes. `source_variable_labels.csv` also covers identifiers, weights, demographic auxiliaries, and CIL plausible values for each country. `sample_counts.csv` records retention by national sample. PV1CIL names the CIL predictor slot; matched analyses substitute the corresponding CIL plausible value for each CT plausible value.

To rebuild the extracts from the official public-use release:

```bash
Rscript code/download_data.R data/official
Rscript code/prepare_data.R data/official/ICILS2018_IDB_SPSS/Data data
```

The full database and documentation are available from the [IEA repository](https://www.iea.nl/data-tools/repository/icils).
