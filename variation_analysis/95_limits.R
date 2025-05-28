model <- "MP" # Choose from "MP", "CF"
target_sex <- 0.0 # Choose from 0.001, 0.0
sim <- read.delim(paste0("../stats/collected/", model, "_collected.tsv"), sep = "\t")
sim <- na.omit(sim)
# sim$SEX2 <- sim$SEX
# sim$SEX2[sim$SEX2 == 0] <- 1e-6
# sim$SEX2 <- log10(sim$SEX2)

# sim$REC2 <- sim$REC
# sim$REC2[sim$REC2 == 0] <- 1e-10
# sim$REC2 <- log10(sim$REC2)

recs <- unique(sim$REC)
secs <- unique(sim$SEX)
print(secs)
to.plot <- c("Hi", "Fis", "LD_0.1000", "S_1_0", "S_0_1", "S_SINGER", "S_IQ_TREE", "D_1_0", "D_0_1", "D_SINGER", "D_IQ_TREE")

p.sex <- array(NA, dim = c(length(recs), length(secs), length(to.plot)))

for (i in (1:length(recs))) {
    subs <- sim$REC == recs[i]
    subs1 <- (sim$REC == recs[i]) & (sim$SEX == target_sex)
    for (j in (1:length(to.plot))) {
        quants <- quantile(sim[, to.plot[j]][subs1], c(0.025, 0.975))
        for (k in (1:length(secs))) {
            subs2 <- (sim$REC == recs[i]) & (sim$SEX == secs[k])
            p.sex[i, k, j] <- sum((sim[, to.plot[j]][subs2] >= quants[1]) & (sim[, to.plot[j]][subs2] <= quants[2])) / sum(subs2)
        }
    }
}

for (j in seq_along(to.plot)) {
    df <- p.sex[, , j]
    rownames(df) <- recs
    colnames(df) <- secs
    write.table(
        df,
        file = paste0("./95_limits_", "SEX_", target_sex, "_", model, "/", to.plot[j], ".txt"),
        sep = "\t", quote = FALSE, col.names = NA
    )
}
