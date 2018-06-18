context("plotting")

test_that("correlation plot", {
  data(OpenSWATH_data, package="SWATH2stats")
  data(Study_design, package="SWATH2stats")
  data.FDR <- sample_annotation(OpenSWATH_data, Study_design)


  correlations <- plot_correlation_between_samples(data.FDR)
  expect_that(correlations[correlations$Var1 == "Hela_Control_1" &
                           correlations$Var2 == "Hela_Control_1","value"], equals(1))
  expect_that(correlations[correlations$Var1 == "Hela_Treatment_1" &
                           correlations$Var2 == "Hela_Treatment_1","value"], equals(1))

  correlations <- plot_correlation_between_samples(data.FDR, label = FALSE)
  expect_that(correlations[correlations$Var1 == "Hela_Control_1" &
                           correlations$Var2 == "Hela_Control_1","value"], equals(1))


  data1 <- data.FDR[data.FDR$condition == "Hela_Control" &
                    data.FDR$bioreplicate == 1 &
                    data.FDR$decoy == 0, c("transition_group_id", "intensity")]
  data2 <- data.FDR[data.FDR$condition == "Hela_Treatment" &
                    data.FDR$bioreplicate == 1 &
                    data.FDR$decoy == 0, c("transition_group_id", "intensity")]
  data3 <- merge(data1, data2, by="transition_group_id")

  cor.p <- cor(data3$intensity.x, data3$intensity.y, method="pearson")
  cor.s <- cor(data3$intensity.x, data3$intensity.y, method="spearman")

  expect_that(correlations[correlations$Var1 == "Hela_Control_1" &
                           correlations$Var2 == "Hela_Treatment_1" &
                           correlations$method == "pearson","value"], equals(cor.p))
  expect_that(correlations[correlations$Var1 == "Hela_Treatment_1" &
                           correlations$Var2 == "Hela_Control_1" &
                           correlations$method == "spearman","value"], equals(cor.s))
})

test_that("variation plot", {
  data(OpenSWATH_data, package="SWATH2stats")
  data(Study_design, package="SWATH2stats")
  data.FDR <- sample_annotation(OpenSWATH_data, Study_design)

  variation <- plot_variation(data.FDR, label=FALSE)
  variation <- plot_variation(data.FDR)
  variation.data <- variation[[1]]

  # test extraction
  test1 <- variation[variation$transition_group_id == "124947_SGWVKPIIIGVLR_2_run0" &
                     variation$condition == "Hela_Control",]
  test2 <- variation[variation$transition_group_id == "214321_GTLNLDSYR_2_run0" &
                     variation$condition == "Hela_Treatment",]

  val1 <- data.FDR[data.FDR$transition_group_id == "124947_SGWVKPIIIGVLR_2_run0" &
                   data.FDR$condition == "Hela_Control", "intensity"]
  val2 <- data.FDR[data.FDR$transition_group_id == "214321_GTLNLDSYR_2_run0" &
                   data.FDR$condition == "Hela_Treatment", "intensity"]

  expect_true(sum(test1[, c("1", "2", "3")] %in% val1) == 3)
  expect_true(sum(test2[, c("1", "2", "3")] %in% val2) == 3)

  # test cv calculation

  mean.test1 <- mean(as.numeric(test1[, c("1", "2", "3")]))
  mean.test2 <- mean(as.numeric(test2[, c("1", "2", "3")]))

  expect_true(sd(as.numeric(test1[, c("1", "2", "3")])) /
              mean.test1 == as.numeric(test1[, "cv"]))
  expect_true(sd(as.numeric(test2[, c("1", "2", "3")])) /
              mean.test2 == as.numeric(test2[, "cv"]))
})

test_that("variation plot vs total", {
  data(OpenSWATH_data, package="SWATH2stats")
  data(Study_design, package="SWATH2stats")
  data.FDR <- sample_annotation(OpenSWATH_data, Study_design)

  variation <- plot_variation_vs_total(data.FDR)
  variation.data <- variation[[1]]

  # test extraction
  test1.total <- variation.data[variation.data$rep == "124947_SGWVKPIIIGVLR_2_run0",]
  test1.rep <- variation.data[variation.data$rep == "124947_SGWVKPIIIGVLR_2_run0 Hela_Control",]

  test2.total <- variation.data[variation.data$rep == "214321_GTLNLDSYR_2_run0",]
  test2.rep <- variation.data[variation.data$rep == "214321_GTLNLDSYR_2_run0 Hela_Treatment",]

  val1.total <- data.FDR[data.FDR$transition_group_id == "124947_SGWVKPIIIGVLR_2_run0", "intensity"]
  val2.total <- data.FDR[data.FDR$transition_group_id == "214321_GTLNLDSYR_2_run0", "intensity"]

  val1.rep <- data.FDR[data.FDR$transition_group_id == "124947_SGWVKPIIIGVLR_2_run0" &
                       data.FDR$condition == "Hela_Control", "intensity"]
  val2.rep <- data.FDR[data.FDR$transition_group_id == "214321_GTLNLDSYR_2_run0" &
                       data.FDR$condition == "Hela_Treatment", "intensity"]

  cv.val1.total <- sd(val1.total) / mean(val1.total)
  cv.val2.total <- sd(val2.total) / mean(val2.total)
  cv.val1.rep <- sd(val1.rep) / mean(val1.rep)
  cv.val2.rep <- sd(val2.rep) / mean(val2.rep)

  expect_true(cv.val1.total == as.numeric(test1.total[, "cv"]))
  expect_true(cv.val2.total == as.numeric(test2.total[, "cv"]))
  expect_true(cv.val1.rep == as.numeric(test1.rep[, "cv"]))
  expect_true(cv.val2.rep == as.numeric(test2.rep[, "cv"]))
})
