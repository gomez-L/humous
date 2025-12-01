
# CELL SELECTION FOR HUMAN RAW PERMUTED HUMOUS - all pipeline since integration

library(Seurat) ; library(readr) ; library(dplyr) ; library(VennDiagram) ; library(gridExtra)


# load and preprocess data
# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- #
# load integrated object - with cell selection only per age
IntH <- readRDS("~/humous/humous_v3/out/IntH/IntH.rds")
# re-annotate age and diff groups as by esther annotations post integration
DiffTypes_Hm_Int <- read_csv("humous_v3/out/IntH/DiffTypes_Hm_Int_new.csv")
IntH$age_ek <- DiffTypes_Hm_Int$age[match(colnames(IntH),DiffTypes_Hm_Int$cells)]
IntH$diff_ek <- DiffTypes_Hm_Int$DiffTypes_final[match(colnames(IntH),DiffTypes_Hm_Int$cells)]
# how many cells we have per age - diff group?
table(IntH$age_ek,IntH$diff_ek)

# how many cells we have to select? as in LandsS_H
# load dataset with cell selection for original humous
LandsS_H <- readRDS("~/humous/humous_v4/out/landsH/LandsS_H.rds")
table(LandsS_H$age_ek,LandsS_H$diff_ek)
# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- #



# perform cell selection for permutations 1 and 2
# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- #

# annotate DiffTypes_Hm_Int with original selection
DiffTypes_Hm_Int$peOrig <- ifelse(DiffTypes_Hm_Int$cells %in% colnames(LandsS_H),TRUE,FALSE)

# how many cells we have for selecting that where not in original selection, per age-diff group?
table(DiffTypes_Hm_Int$age,DiffTypes_Hm_Int$DiffTypes_final,DiffTypes_Hm_Int$peOrig)

# random selection for pe1, with 490 cells for all age-diff combinations except for 12-IPC, where we will take 320

# add stable row id so we can mark rows after sampling
DiffTypes_Hm_Int <- DiffTypes_Hm_Int %>% mutate(.rowid = row_number())

# 1) Check group sizes vs required sizes
group_req <- DiffTypes_Hm_Int %>%
  count(age, DiffTypes_final, name = "group_n") %>%
  mutate(required = if_else(age == 12 & DiffTypes_final == "IPC", 320L, 490L))

# stop with informative message if any group is too small
bad <- filter(group_req, group_n < required)
if (nrow(bad) > 0) {
  stop(
    "Some (age, DiffTypes_final) groups don't have enough rows to sample.\nOffending groups:\n",
    paste0(
      "  age=", bad$age,
      ", DiffTypes_final='", bad$DiffTypes_final,
      "' — have ", bad$group_n,
      ", need ", bad$required,
      collapse = "\n"
    )
  )
}

# 2) Do first stratified selection (per-group sizes). Use group_modify so we can compute a scalar n per group
set.seed(123)   # reproducible draw 1
sel1_ids <- DiffTypes_Hm_Int %>%
  group_by(age, DiffTypes_final) %>%
  group_modify(~ {
    n_to_sample <- if (.y$age == 12 && .y$DiffTypes_final == "IPC") 320L else 490L
    dplyr::slice_sample(.x, n = n_to_sample)
  }) %>%
  ungroup() %>%
  pull(.rowid)

# 3) Do second independent stratified selection (different seed so draws differ)
set.seed(456)   # reproducible draw 2 (different seed => independent draw)
sel2_ids <- DiffTypes_Hm_Int %>%
  group_by(age, DiffTypes_final) %>%
  group_modify(~ {
    n_to_sample <- if (.y$age == 12 && .y$DiffTypes_final == "IPC") 320L else 490L
    dplyr::slice_sample(.x, n = n_to_sample)
  }) %>%
  ungroup() %>%
  pull(.rowid)

# 4) Annotate the original dataframe with two logical selection flags
DiffTypes_Hm_Int_annotated <- DiffTypes_Hm_Int %>%
  mutate(
    sel1 = .rowid %in% sel1_ids,
    sel2 = .rowid %in% sel2_ids
  ) %>%
  select(-.rowid)   # drop helper id if you want

#write.csv(DiffTypes_Hm_Int_annotated,"humous_v4/out/raw_permuted_humous/cellselection_RPH/DiffTypes_Hm_Int_annotated.csv")

# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- #


# Evaluate cell selection results
# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- #


# check selection cell numbers and how many cells are the same or different across selections
table(DiffTypes_Hm_Int_annotated$age,DiffTypes_Hm_Int_annotated$DiffTypes_final,DiffTypes_Hm_Int_annotated$peOrig)
table(DiffTypes_Hm_Int_annotated$age,DiffTypes_Hm_Int_annotated$DiffTypes_final,DiffTypes_Hm_Int_annotated$sel1)
table(DiffTypes_Hm_Int_annotated$age,DiffTypes_Hm_Int_annotated$DiffTypes_final,DiffTypes_Hm_Int_annotated$sel2)

intersection_summary <- DiffTypes_Hm_Int_annotated %>%
  mutate(
    intersection = case_when(
      peOrig & !sel1 & !sel2 ~ "peOrig only",
      !peOrig & sel1 & !sel2 ~ "sel1 only",
      !peOrig & !sel1 & sel2 ~ "sel2 only",
      peOrig & sel1 & !sel2 ~ "peOrig + sel1",
      peOrig & !sel1 & sel2 ~ "peOrig + sel2",
      !peOrig & sel1 & sel2 ~ "sel1 + sel2",
      peOrig & sel1 & sel2 ~ "peOrig + sel1 + sel2",
      TRUE ~ "none"
    )
  ) %>%
  group_by(intersection, age, DiffTypes_final) %>%
  summarise(n = n(), .groups = "drop") %>%
  arrange(intersection, desc(n))
#write.csv(intersection_summary,"humous_v4/out/raw_permuted_humous/cellselection_RPH/intersection_summary.csv")


# venn diagram with all age-diff groups together
venn.plot <- draw.triple.venn(
  area1 = sum(DiffTypes_Hm_Int_annotated$peOrig),
  area2 = sum(DiffTypes_Hm_Int_annotated$sel1),
  area3 = sum(DiffTypes_Hm_Int_annotated$sel2),
  n12 = sum(DiffTypes_Hm_Int_annotated$peOrig & DiffTypes_Hm_Int_annotated$sel1),
  n23 = sum(DiffTypes_Hm_Int_annotated$sel1 & DiffTypes_Hm_Int_annotated$sel2),
  n13 = sum(DiffTypes_Hm_Int_annotated$peOrig & DiffTypes_Hm_Int_annotated$sel2),
  n123 = sum(DiffTypes_Hm_Int_annotated$peOrig & DiffTypes_Hm_Int_annotated$sel1 & DiffTypes_Hm_Int_annotated$sel2),
  category = c("peOrig", "sel1", "sel2"),
  fill = c("red", "blue", "green"),
  alpha = 0.5,
  cex = 1.5,
  cat.cex = 1.5
)
grid.draw(venn.plot)


# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- #
