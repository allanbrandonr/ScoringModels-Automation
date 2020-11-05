# ScoringModels-Automation
Providing functions for automate selection and fitting of scoring a& weighted regression based models

# Summary
The objective of the codes is to automate scoring model selection and fitting. The codes optimize fitting of a weighted regresion model 
on intervals constructed with a trend wize approach to every variable. The codes also provide the prediction procedure for new data.

Depending on the type of variable, different performance metrics and statistical significance tests are available.

# The codes
- trameado: Estimates the optimum weighted interval partition for each explanatory variable in the set. This partition is optimum on the performance metric chosen. It also provide the prediction procedure for new data.
- modelado: Fits a weighted regression model and makes an optimum variable backward selection conidering minimum count of variables, collinearity, variable weight thresholds and significance
- transformaciones: Apply transformations to re shape each variable or shifts to consider previous information.

