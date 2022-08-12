alpha1 = 0.05
maf = 0.001
alpha2 = 0.50
cov = 10
rm("test/test_1.syncx")
rm(string("test/test_1-FILTERED-alpha1_", alpha1, "-maf_", maf, "-alpha2_", alpha2, "-cov_", cov, ".pileup"))
rm(string("test/test_1-FILTERED-alpha1_", alpha1, "-maf_", maf, "-alpha2_", alpha2, "-cov_", cov, "-FROM_PILEUP.syncx"))
rm(string("test/test_1-FILTERED-alpha1_", alpha1, "-maf_", maf, "-alpha2_", alpha2, "-cov_", cov, ".syncx"))
# rm("test/test_1-IMPUTED-window_20-model_OLS-distance_true.syncx")

alpha1 = 0.05
maf = 0.0001
alpha2 = 0.50
cov = 1
rm("test/test_2.syncx")
rm(string("test/test_2-FILTERED-alpha1_", alpha1, "-maf_", maf, "-alpha2_", alpha2, "-cov_", cov, ".pileup"))
rm(string("test/test_2-FILTERED-alpha1_", alpha1, "-maf_", maf, "-alpha2_", alpha2, "-cov_", cov, "-FROM_PILEUP.syncx"))
rm(string("test/test_2-FILTERED-alpha1_", alpha1, "-maf_", maf, "-alpha2_", alpha2, "-cov_", cov, ".syncx"))
# rm("test/test_2-IMPUTED-window_20-model_OLS-distance_true.syncx")
