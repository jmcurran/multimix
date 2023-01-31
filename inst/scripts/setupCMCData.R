cmc.df = read.csv(system.file("extdata", "cmc.csv", package = "multimix"),
                       header = FALSE)

names(cmc.df) =  c("age",
                   "edu",
                   "eduh",
                   "nborn",
                   "islam",
                   "working",
                   "husocc",
                   "sol",
                   "medex",
                   "method")

library(dplyr)
cmc.df = cmc.df |> 
  mutate(edu = recode_factor(as.factor(edu),
                `1` = "low",
                `2` = "below",
                `3` = "above",
                `4` = "high"), 
        eduh = recode_factor(as.factor(eduh),
                `1` = "low",
                `2` = "below",
                `3` = "above",
                `4` = "high"),
        sol = recode_factor(as.factor(sol),
                             `1` = "low",
                             `2` = "below",
                             `3` = "above",
                             `4` = "high"),
        method = recode_factor(as.factor(method),
                                         `1` = "None",
                                         `2` = "Long.term",
                                         `3` = "Short.term")) |> 
  mutate(islam = factor(case_when(islam == 0 ~ "Non.Islam",
                                  TRUE ~ "Islam"))) |> 
  mutate(working = factor(case_when(working == 0 ~ "Yes",
                                    TRUE ~ "No"))) |> 
  mutate(husocc = factor(paste0("ho", husocc))) |> 
  mutate(medex = factor(case_when(medex == 0 ~ "Good",
                                  TRUE ~ "Not.Good")))


save(cmc.df, file = "cmc.df.rda")
save(cmc.df, file = "src/R/data/cmc.df.rda")

