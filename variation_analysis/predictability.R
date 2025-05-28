### Predictability of sex rate from different statistics
model <- "MP" # Choose from "MP", "CF"
sim <- read.delim(paste0("../stats/collected/", model, "_collected.tsv"), sep = "\t")
# colSums(is.na(sim)) # check for NAs
sim <- na.omit(sim)

## make a continuous variable of log10(rate of sex), substituting 0 for 1e-6
sim$SEX2 <- sim$SEX
sim$SEX2[sim$SEX2 == 0] <- 1e-6
sim$SEX2 <- log10(sim$SEX2)

## same for log10(recombination rate) - not actually using that here
sim$REC2 <- sim$REC
sim$REC2[sim$REC2 == 0] <- 1e-10
sim$REC2 <- log10(sim$REC2)

## values of REC
recs <- rev(unique(sim$REC))
## metrics to plot - column headings
to.plot <- c("Hi", "Fis", "LD_0.1000", "S_1_0", "S_0_1", "S_SINGER", "S_IQ_TREE", "D_1_0", "D_0_1", "D_SINGER", "D_IQ_TREE")


# ## store Rsquare values
rsq.vals <- matrix(NA, nrow = length(recs), ncol = length(to.plot))
colnames(rsq.vals) <- to.plot
rownames(rsq.vals) <- paste0("REC=", recs)
# rsq.vals

# ## make plots
for (i in (length(recs):1)) {
    subs <- sim$REC == recs[i]
    for (j in (1:length(to.plot))) {
        m <- lm(sim$SEX2[subs] ~ poly(sim[, to.plot[j]][subs], 2))
        rsq.vals[i, j] <- round(summary(m)$r.squared, 2)
    }
}

## table of Rsq
rsq.vals
write.table(rsq.vals, file = paste0("R2_tables/", model, ".tsv"), sep = "\t", quote = FALSE, col.names = TRUE)
