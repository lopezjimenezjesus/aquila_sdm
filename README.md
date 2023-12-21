# Model 
Aquila adalberti Favourability and Distribution Modeling


# Instructions

1. Edit the config file to set up the output folder, VIF threshold and rest of parameters
2. Open the RStudio project.
3. Open render_report.R and run it
4. Check the report in the 04_RESULTS folder

# Folders

00_DOC: contain previous examples and monographies
01_DATA: primary species dataset and predictor variables
02_SCRIPTS: workflow and functions used to run the models 
03_BIBLIOGRAPHY: .bib format used to render bibliography (currently not used for rendering the provisional report)
04_RESULTS: contain the  intermediate results, figures and reports

# Note

- Most important scripts start with 02, 03 and 05. 01 and 04 are mean to just explore the relationship between variables and explore new modelling methods. Not of them are executed when running render_report.R.
- Depending on the value of the parameters the intermediate results, figures and report are saved in a different folder. 
- While metrics are similar for different species, the report should be adapted as differences in the interpretation are expected (obviously)

