rm(list = ls())

library(readxl)
library(dplyr)

setwd("/Users/xiyu/Desktop/network")
GDP <- read_excel("GDP.xlsx")

GDP_clean <- GDP %>%
  select(地区, 年份, `地区生产总值(亿元)`, `人均地区生产总值(元/人)`)
GDP_clean <- GDP_clean %>%
  filter(年份 %in% c(1995:2000,2005:2010, 2015:2020))

GDP_clean <- GDP_clean %>%
  mutate(YearGroup = case_when(
    年份 >= 1995 & 年份 <= 2000 ~ "1995-2000",
    年份 >= 2005 & 年份 <= 2010 ~ "2005-2010",
    年份 >= 2015 & 年份 <= 2020 ~ "2015-2020",
    TRUE ~ NA_character_)) %>%
  filter(!is.na(YearGroup)) %>%
  group_by(地区, YearGroup) %>%
  mutate(pGDP = mean(`人均地区生产总值(元/人)`, na.rm = TRUE)) %>%
  ungroup()

GDP_summary <- GDP_clean %>%
  select(地区, YearGroup, pGDP) %>%
  distinct(地区, YearGroup, .keep_all = TRUE)

GDP_summary <- GDP_summary %>%
  arrange(desc(YearGroup))


write.xlsx(GDP_summary, "GDP_summary.xlsx")


