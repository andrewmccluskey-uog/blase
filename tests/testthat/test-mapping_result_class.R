## Tests show
test_that("mapping_result_class show() prints what we expect", {
    mapping_result <- methods::new("MappingResult",
        bulk_name = "Test",
        best_bin = 1,
        best_correlation = 0.05,
        top_2_distance = 0.001,
        confident_mapping = TRUE,
        history = data.frame(bin = c(1), correlation = c(0.5), lower_bound = c(0.3), upper_bound = c(0.6)),
        bootstrap_iterations = 200
    )

    expect_output(
        print(mapping_result),
        "best_bin=1 correlation=0.05 top_2_distance=0.001\n\t Confident Result: TRUE (next max upper  -Inf )\n\t with history for scores against 1  bins\n\t Bootstrapped with 200 iterations",
        fixed = TRUE
    )
})

## Test getters
test_that("mapping_result_class show() getter(bulk_name) works", {
    mapping_result <- methods::new("MappingResult",
        bulk_name = "Test",
        best_bin = 1,
        best_correlation = 0.05,
        top_2_distance = 0.001,
        confident_mapping = TRUE,
        history = data.frame(bin = c(1), correlation = c(0.5), lower_bound = c(0.3), upper_bound = c(0.6)),
        bootstrap_iterations = 200
    )

    expect_equal(bulk_name(mapping_result), "Test")
})

test_that("mapping_result_class show() getter(best_bin) works", {
    mapping_result <- methods::new("MappingResult",
        bulk_name = "Test",
        best_bin = 1,
        best_correlation = 0.05,
        top_2_distance = 0.001,
        confident_mapping = TRUE,
        history = data.frame(bin = c(1), correlation = c(0.5), lower_bound = c(0.3), upper_bound = c(0.6)),
        bootstrap_iterations = 200
    )

    expect_equal(best_bin(mapping_result), 1)
})

test_that("mapping_result_class show() getter(best_correlation) works", {
    mapping_result <- methods::new("MappingResult",
        bulk_name = "Test",
        best_bin = 1,
        best_correlation = 0.05,
        top_2_distance = 0.001,
        confident_mapping = TRUE,
        history = data.frame(bin = c(1), correlation = c(0.5), lower_bound = c(0.3), upper_bound = c(0.6)),
        bootstrap_iterations = 200
    )

    expect_equal(best_correlation(mapping_result), 0.05)
})

test_that("mapping_result_class show() getter(top_2_distance) works", {
    mapping_result <- methods::new("MappingResult",
        bulk_name = "Test",
        best_bin = 1,
        best_correlation = 0.05,
        top_2_distance = 0.001,
        confident_mapping = TRUE,
        history = data.frame(bin = c(1), correlation = c(0.5), lower_bound = c(0.3), upper_bound = c(0.6)),
        bootstrap_iterations = 200
    )

    expect_equal(top_2_distance(mapping_result), 0.001)
})


test_that("mapping_result_class show() getter(confident_mapping) works", {
    mapping_result <- methods::new("MappingResult",
        bulk_name = "Test",
        best_bin = 1,
        best_correlation = 0.05,
        top_2_distance = 0.001,
        confident_mapping = TRUE,
        history = data.frame(bin = c(1), correlation = c(0.5), lower_bound = c(0.3), upper_bound = c(0.6)),
        bootstrap_iterations = 200
    )

    expect_equal(confident_mapping(mapping_result), TRUE)
})

test_that("mapping_result_class show() getter(mapping_history) works", {
    history <- data.frame(bin = c(1), correlation = c(0.5), lower_bound = c(0.3), upper_bound = c(0.6))
    mapping_result <- methods::new("MappingResult",
        bulk_name = "Test",
        best_bin = 1,
        best_correlation = 0.05,
        top_2_distance = 0.001,
        confident_mapping = TRUE,
        history = history,
        bootstrap_iterations = 200
    )

    expect_equal(mapping_history(mapping_result), history)
})

test_that("mapping_result_class show() getter(bootstrap_iterations) works", {
    mapping_result <- methods::new("MappingResult",
        bulk_name = "Test",
        best_bin = 1,
        best_correlation = 0.05,
        top_2_distance = 0.001,
        confident_mapping = TRUE,
        history = data.frame(bin = c(1), correlation = c(0.5), lower_bound = c(0.3), upper_bound = c(0.6)),
        bootstrap_iterations = 200
    )

    expect_equal(bootstrap_iterations(mapping_result), 200)
})
