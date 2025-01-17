suppressMessages({
	library(dplyr)
	library(data.table)
})

#import Olink data in wide and remove Nan
targets <- read.table("~/sean/INTRePID/STARNET/proteins/STARNET_Olink_20180383_20181235_cases_combined.tsv", sep='\t', header=TRUE)
targets[targets$status=="", ] <- NA
targets_na <- na.omit(targets)

targets_filt <- subset(targets_na, select = -c(starnet_id, status, Sex, Age, year))
targets_scale <- data.frame(scale(targets_filt))

dat_adj <- function(df) {
	fit <- lm(df ~ Age + Sex + year, data = targets_na)
	return(resid(fit))
}

adj_res <- lapply(targets_scale, dat_adj)
df_adj <- data.frame(adj_res)
cols <- subset(targets_na, select = c(starnet_id, status, Sex, Age, year))
out <- cbind(cols, df_adj)

write.table(out, "~/sean/INTRePID/STARNET/proteins/adjusted/STARNET_Olink_20180383_20181235_cases_combined_adj-age-sex-year.tsv", row.names=FALSE, sep="\t", quote = FALSE)
