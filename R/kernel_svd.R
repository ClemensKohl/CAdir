# https://bobrupakroy.medium.com/what-is-kernel-pca-using-r-python-4864c2471e62

# FIXME: WIP
kernel_ca <- function(data, dims = 2) {
    kpca = kernlab::kpca(~.,
                         data = data,
                         kernel = 'rbfdot',
                         features = dims)

    data_pca = as.data.frame(predict(kpca, data))
    head(training_set_pca)
    training_set_pca$Customer_Segment = training_set$Customer_Segmenttest_set_pca = as.data.frame(predict(kpca, test_set))
    test_set_pca$Customer_Segment = test_set$Customer_Segment
}
