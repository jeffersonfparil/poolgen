alpha1 = 0.05
maf = 0.001
alpha2 = 0.50
cov = 10
window_size=20
model=["Mean", "OLS", "RR", "LASSO", "GLMNET"][2]
distance=true
rm("test/test_1.syncx")
rm(string("test/test_1-FILTERED-alpha1_", alpha1, "-maf_", maf, "-alpha2_", alpha2, "-cov_", cov, ".pileup"))
rm(string("test/test_1-FILTERED-alpha1_", alpha1, "-maf_", maf, "-alpha2_", alpha2, "-cov_", cov, "-FROM_PILEUP.syncx"))
rm(string("test/test_1-FILTERED-alpha1_", alpha1, "-maf_", maf, "-alpha2_", alpha2, "-cov_", cov, ".syncx"))
rm(string("test/test_1-IMPUTED-window_", window_size, "-model_", model, "-distance_", distance, ".syncx"))

alpha1 = 0.05
maf = 0.0001
alpha2 = 0.50
cov = 1
window_size=20
model=["Mean", "OLS", "RR", "LASSO", "GLMNET"][2]
distance=true
rm("test/test_2.syncx")
rm(string("test/test_2-FILTERED-alpha1_", alpha1, "-maf_", maf, "-alpha2_", alpha2, "-cov_", cov, ".pileup"))
rm(string("test/test_2-FILTERED-alpha1_", alpha1, "-maf_", maf, "-alpha2_", alpha2, "-cov_", cov, "-FROM_PILEUP.syncx"))
rm(string("test/test_2-FILTERED-alpha1_", alpha1, "-maf_", maf, "-alpha2_", alpha2, "-cov_", cov, ".syncx"))
rm(string("test/test_2-IMPUTED-window_", window_size, "-model_", model, "-distance_", distance, ".syncx"))
