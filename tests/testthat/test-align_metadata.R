test_that(".mbm_align_metadata reorders metadata by ID and says so", {
    meta <- data.frame(SampleID = c("S3", "S1", "S2"), Group = c("c", "a", "b"))

    expect_message(
        aligned <- .mbm_align_metadata(c("S1", "S2", "S3"), meta),
        "reordered"
    )
    expect_equal(aligned$SampleID, c("S1", "S2", "S3"))
    expect_equal(aligned$Group, c("a", "b", "c"))
})

test_that(".mbm_align_metadata is silent when metadata is already in order", {
    meta <- data.frame(SampleID = c("S1", "S2", "S3"), Group = c("a", "b", "c"))
    expect_silent(.mbm_align_metadata(c("S1", "S2", "S3"), meta))
})

test_that(".mbm_align_metadata never matches by position", {
    meta <- data.frame(ID = c("X1", "X2", "X3"), Group = c("a", "b", "c"))
    expect_error(.mbm_align_metadata(c("S1", "S2", "S3"), meta), "None of the samples")
})

test_that(".mbm_align_metadata leaves out table samples missing from metadata", {
    meta <- data.frame(SampleID = c("S1", "S3"), Group = c("a", "c"))
    expect_message(
        aligned <- .mbm_align_metadata(c("S1", "S2", "S3"), meta),
        "left out: S2"
    )
    expect_equal(aligned$SampleID, c("S1", "S3"))
})

test_that(".mbm_align_metadata rejects duplicated IDs", {
    meta <- data.frame(SampleID = c("S1", "S1", "S2"), Group = c("a", "b", "c"))
    expect_error(.mbm_align_metadata(c("S1", "S2"), meta), "Duplicated")
})
