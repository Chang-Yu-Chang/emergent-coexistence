# Test the example code and data.
# Fit Sanger sequences of a mixed culture of four strains. This step may take a few seconds.

library(CASEU)
data('fourStrainExpt') # load an example dataset

results <- fitSangerMixture(mixture = fourStrainExpt$mixture, # mixture
                            components = fourStrainExpt$components, # Four individual 16S sequences
                            verbose = TRUE)

save(results, file = "~/Dropbox/lab/emergent-coexistence/data/temp/CASEU_test_results1.Rdata")


load("~/Dropbox/lab/emergent-coexistence/data/temp/CASEU_test_results1.Rdata")
print(results)                      # view the results
plot(results)                       # plot the fit
sangerFitDiagnosticPlot(results)    # plot a bunch more diagnostics
names(results)
