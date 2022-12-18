cancer.df = read.table(system.file("extdata", "CANCER11a.DAT", package = "multimix"),
                       header = FALSE)
names(cancer.df) = c("Age", 
                     "Wt", 
                     "PF", 
                     "HX", 
                     "SBP",
                     "DBP",
                     "EKG",
                     "HG",
                     "SZ",
                     "SG",
                     "AP",
                     "BM")

library(dplyr)

cancer.df = cancer.df |> 
  mutate(PF = recode_factor(as.factor(PF),
                     `0` = "Active",
                     `1` = "Bed49",
                     `2` = "Bed51")) |> 
  mutate(HX = recode_factor(as.factor(HX),
                     `0` = "No_hist",
                     `1` = "Hist")) |> 
  mutate(EKG = recode(as.factor(EKG),
                      `0` = "Normal",
                      `1` = "Benign",
                      `2` = "Rythmic",
                      `3` = "H_blocks",
                      `4` = "H_strain",
                      `5` = "Old_MC1",
                      `6` = "New_MC1")
                      )

save(cancer.df, file = "data/cancer.df.rda")
