# WeightedPrincipalSVM
Implement principal weighted support vector machine by Shin et al. (2017). The function takes binary outcome and continuous covariates, and returns the central space of the covariates that predicts the outcome, estimated by the principal weighted svm. The central space can be used for dimension reduction.

# Dependency
Runs in R. Loads libraries `MASS` and `optimx`.

Shin, Seung Jun, et al. "Principal weighted support vector machines for sufficient dimension reduction in binary classification." Biometrika 104.1 (2017): 67-81.
