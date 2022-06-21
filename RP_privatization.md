Replication Package for ‘Does homeownership hinder labor market
activity? Evidence from housing privatization and restitution’
================
Stepan Mikula and Josef Montag
June 18, 2019

Contacts:

-   <stepan.mikula@econ.muni.cz> \|
    [https://sites.google.com/site/stepanmikula](https://sites.google.com/site/stepanmikula/).
-   <josef.montag@gmail.com> \|
    [https://sites.google.com/site/josefmontag](https://sites.google.com/site/josefmontag/).

This package file contains R script and outputs reported in Mikula,
Štěpán, and Josef Montag. 2019. Does homeownership hinder labor market
activity? Evidence from housing privatization and restitution. MUNI ECON
Working Paper No. 2019-06. Brno, Czech Republic: Department of
Economics, Masaryk University.

The manuscript is available at
<https://ideas.repec.org/p/mub/wpaper/2019-06.html>.

#### Data availability

*Privatized_houses.xls* is the database of privatized houses in Brno.
For each house it contains address and the year of privatization. The
dataset analyzed in the paper is built up on this database as described
in the paper.

The individual-level census data used in this paper were provided by the
Czech Statistical Office (CSO) under the condition of confidentiality
and thus cannot be provided as part of this package. However, the data
can be requested from the CSO and the authors can provide help with
that.

``` r
#Load libraries and functions, set constants and load data.

library(tidyverse)
library(rgdal)
library(broom)
library(ggmap)
library(stargazer)
library(lfe)
library(lmtest)
library(sandwich)
library(AER)
library(ivpack)
library(Formula)
library(geosphere)
library(DescTools)
library(readxl)
library(sf)
library(ggsn)
library(ggspatial)
library(ggmap)
library(cowplot)
library(sf)
library(knitr)

rm(list = ls())

# Load functions
source("CODE/functions.R")

# Constants

## Filtering
age_upper <- 61 #56
age_lower <- 19 #25
house_size <- 6 #6 minimum house size

## Matching
d_caliper <- 0.7
d_distance <- 700

## Robust clustered SE
cluster_by <- "house"

## Output format
table_type <- "text"

# Load data
load("DATA_RP/RPdata.RData")

# Number of digits
dig <- 3
```

Variables in `RPdata.RData`:

-   `brno`…Map of Brno (sf polygon)

-   `rzsj`…Map of neighborhoods (sf polygon)

-   `rzsj_body`…Centroides of neighborhoods (sf point)

-   `rzsj_houses`…No. of houses per group (type) and neighborhood (df)

-   `r91_data_raw`…Census data from 1991 (no filtering applied)

-   `r91_data`…1991 Census data (houses matching with 2001 census)

-   `desc_91`…1991 Census data pre-processed for descriptives rendering

-   `r91_houses`…1991 Census data aggregated by house

-   `r01_data_raw`…2001 Census data (no filtering applied)

-   `r01_data_basic`…2001 Census data (estimation sample)

-   `r01_data`…2001 Census data (houses matching with 1991 census)

Other variables stored in separate files:

-   `IDOB`…table for IDs conversion
-   `privyear`…exact year of privatization
-   `zsjsf`…equivalent of `rzsj` from `RPdata.RData`
-   `rzsj` from `rzsj_centroides.RData`…centroides of neighborhoods
    (different format, simple df)

Key variables:

-   `ownership`…owner/renter
-   `RZSJ`…neighborhood code
-   `EDUC`…education level
-   `LPOHLAV`/`MALE`…gender
-   `htype`…household type
-   `BornBrno`/`BB`…individual born in Brno
-   `LVEK`…age
-   `LEKAKTI`…economic activity (1 = employed, 5 = unemployed)
-   `JPOCBYT`/`byty`…house size
-   `sample`…group of houses (ATT = privatized before 2001, NTR =
    restituted, NTP = privatized in 2001 or later, NTC = city-owned and
    never privatized)

## Formulas

``` r
# Individual data models
model <- list(
  I(LEKAKTI == 5) ~ ownership,
  I(LEKAKTI == 5) ~ ownership + factor(RZSJ),
  I(LEKAKTI == 5) ~ ownership + EDUC,
  I(LEKAKTI == 5) ~ ownership + EDUC + LPOHLAV  + htype + BornBrno + factor(LVEK),
  I(LEKAKTI == 5) ~ ownership + EDUC + LPOHLAV  + htype + BornBrno + factor(RZSJ) + factor(LVEK)
) %>% lapply(as.Formula)

# Change dependent variable from unemployment to economic activity
model <- c(
  model,
  #model %>% lapply(update, I(LEKAKTI == 1) ~ .),
  model %>% lapply(update, I(LEKAKTI %in% c(1,5)) ~ .)
)

model_individual <- model

# Individual data models, gender-specific subsamples (not used)

model_individual_gender <- model_individual %>% lapply(update, .~. - LPOHLAV)

# Individual data models (IV)

modeliv_individual <- list(
  I(LEKAKTI == 5) ~ ownership | sample_dummy,
  I(LEKAKTI == 5) ~ ownership + factor(RZSJ) | sample_dummy + factor(RZSJ),
  I(LEKAKTI == 5) ~ ownership + EDUC | sample_dummy + EDUC,
    I(LEKAKTI == 5) ~ ownership + EDUC + LPOHLAV  + htype + BornBrno + factor(LVEK) | sample_dummy + EDUC + LPOHLAV  + htype + BornBrno + factor(LVEK),
  I(LEKAKTI == 5)  ~ ownership + EDUC + LPOHLAV  + htype + BornBrno + factor(RZSJ) + factor(LVEK) | sample_dummy + EDUC + LPOHLAV + htype + BornBrno + factor(RZSJ) + factor(LVEK)
) %>% lapply(as.Formula)

modeliv_1ststage <- list(
  ownership ~ sample_dummy,
  ownership ~ sample_dummy + factor(RZSJ),
  ownership ~ sample_dummy + EDUC,
  ownership ~ sample_dummy + EDUC + LPOHLAV + factor(LVEK) + htype + BornBrno,
  ownership ~ sample_dummy + EDUC + LPOHLAV + factor(LVEK) + htype + BornBrno + factor(RZSJ)
) %>% lapply(as.Formula)

# Change dependent variable from unemployment to economic activity

modeliv_individual <- c(
  modeliv_individual,
  modeliv_individual %>% lapply(update, I(LEKAKTI %in% c(1,5)) ~ .)
)

# Individual data models, gender-specific subsamples

modeliv_individual_gender <- modeliv_individual %>% lapply(update, .~. - LPOHLAV | . - LPOHLAV)

# Decision to privatize

model_fs <- list(
  I(sample != "NTC") ~ LVEK + MALE + EDUCmiddle + EDUChigh + BornBrno + byty,
  I(sample != "NTC") ~ LVEK + MALE + EDUCmiddle + EDUChigh + BornBrno + byty + factor(RZSJ)
) %>% lapply(as.Formula)

model_fs <- c(
  model_fs,
  model_fs %>% lapply(update, I(sample == "ATT") ~ .)
)

model_propscore <- as.Formula(X ~ EDUChigh + EDUCmiddle + BornBrno + LVEK + MALE + JPOCBYT + U + EA)
```

## Summary statistics

### Maps

#### Figure 3a

``` r
# Calculate no. of houses per neighborhood (ZSJ)
rzsj_houses <- r01_data_basic %>%
  select(adresado_01_klic,sample,RZSJ) %>% 
  distinct() %>% 
  #group_by(sample,RZSJ) %>%
  group_by(RZSJ) %>% 
  summarise(
    obs = n()
  ) %>% 
  ungroup() %>% 
  #spread(sample,obs) %>% 
  mutate(
    KOD_ZSJ = str_sub(RZSJ,1,6)
  )

rzsj %>% 
  rename(geometry = SHAPE) %>% 
  left_join(rzsj_houses) %>% 
  ggplot() +
  geom_sf(
    aes(fill = obs),
    size = 0.3
  ) +
  geom_sf(
    data = brno,
    fill = NA,
    size = 1,
    color = "black"
  ) +
  annotation_scale() +
  coord_sf(
    datum = NA
  ) +
  scale_fill_distiller(
    name = "Number of houses",
    palette = "YlOrRd",
    na.value = "white",
    guide = guide_colorbar(
      direction = "horizontal",
      title.position = "top" 
        )
  ) +
  theme_void(
    base_family = "times"
  ) +
  theme(
    legend.position = "bottom"
  )
```

![](RP_privatization_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

#### Figure 3b

``` r
bbox <- rzsj %>% 
  rename(geometry = SHAPE) %>% 
  left_join(rzsj_houses) %>% 
  filter(!is.na(obs)) %>% 
  st_bbox()
  
rzsj_houses <- r01_data_basic %>%
  filter(sample %in% c("ATT","NTR")) %>% 
  select(adresado_01_klic,sample,RZSJ) %>% 
  distinct() %>% 
  group_by(sample,RZSJ) %>%
  #group_by(RZSJ) %>% 
  summarise(
    obs = n()
  ) %>% 
  ungroup() %>% 
  spread(sample,obs) %>% 
  mutate(
    KOD_ZSJ = str_sub(RZSJ,1,6)
  )

ATTbody <- left_join(rzsj_body,rzsj_houses) %>% 
  rename(geomtry = SHAPE) %>% 
  filter(!is.na(ATT))

NTRshape <- rzsj %>% 
  rename(geometry = SHAPE) %>% 
  left_join(rzsj_houses)

NTRshape %>% 
  ggplot() +
  geom_sf(
    aes(fill = NTR),
    size = 0.3
  ) +
  geom_sf(
    data = ATTbody,
    aes(size = ATT),
    show.legend = "point",
    alpha = 0.6
  ) +
  geom_sf(
    data = brno,
    fill = NA,
    size = 1,
    color = "black"
  ) +
  annotation_scale() +
  coord_sf(
    datum = NA,
    xlim = bbox[c(1,3)],
    ylim = bbox[c(2,4)]
  ) +
  scale_fill_distiller(
    name = "Number of restituted\nhouses",
    palette = "YlOrRd",
    na.value = "white",
    guide = guide_colorbar(
      direction = "horizontal",
      title.position = "top" 
        )
  ) +
  scale_size(
    name = "Number of houses\nprivatized before 2001",
    guide = guide_legend(
      direction = "horizontal",
      title.position = "top",
      byrow = TRUE
        )
  ) +
  theme_void(
    base_family = "times"
  ) +
  theme(
    legend.direction = "vertical"
  ) +
  theme(
    legend.position = "bottom"
  )
```

![](RP_privatization_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

#### Table 1: Apartment buildings in Brno in 1991, means by privatization and restitution status

``` r
desc_91 %>% 
  mutate(
    EDUClow = as.integer(EDUC == "low"),
    EDUCmiddle = as.integer(EDUC == "middle"),
    EDUChigh = as.integer(EDUC == "high")
  ) %>% 
  select(-EDUC) %>% 
  group_by(adresado_01_klic,sample) %>%
  summarise_if(is.numeric, mean, na.rm = TRUE) %>% 
  ungroup()  -> desc_91_houses_klic

desc_91_houses_klic %>% 
  select(-adresado_01_klic) -> desc_91_houses

desc_91_houses <- r91_houses %>% 
  select(-adresado_01_klic)

## City-owned
CityOwned <- desc_91_houses %>% 
  select(-Act) %>% 
  filter(sample != "NTR") %>% 
  desctable() %>% 
  format_desctable() %>% 
  rename(CityOwned = value)

## Restituted
Restituted <- desc_91_houses %>% 
  select(-Act) %>% 
  filter(sample == "NTR") %>% 
  desctable() %>% 
  format_desctable() %>% 
  rename(Restituted = value)

# Diff (1) and (2)
desc_91_houses %>% 
  select(-Act) %>% 
  mutate(
    sample = case_when(
      sample != "NTR" ~ "CO",
      TRUE ~ "NTR"
    )
  ) %>% 
  gather(variable,value,-sample) %>%
  mutate(
    sample = factor(sample) %>% relevel(ref = "CO")
  ) %>% 
  #drop_na(value) %>% 
  mutate(value = as.double(value)) %>% 
  group_by(variable) %>% 
  do(
    estmlm = lm(
      value ~ sample,
      data = .
    ) %>% 
      coeftest(vcov. = vcovHC) %>% tidy() %>% 
      filter(str_detect(term,"sample"))
  ) %>% 
  ungroup() %>% 
  unnest(estmlm) %>% 
  rowwise() %>% 
  mutate(
    mean = estimate %>% format(digits = 1, nsmall = dig, trim = TRUE, scientific = FALSE) %>% 
      str_c(.,add_stars(p.value, latex = TRUE)),
    se = std.error %>% format(digits = 1, nsmall = dig, trim = TRUE, scientific = FALSE) %>%
      str_c("(",.,")")
  ) %>% 
  ungroup() %>% 
  select(variable,mean,se) %>% 
  gather(stat,dif1,-variable) -> dif1


## Subgroups
SelectedForP <- desc_91_houses %>%
  filter(sample != "NTR") %>% 
  filter(!is.na(Act)) %>% 
  select(-Act) %>% 
  desctable() %>% 
  format_desctable() %>% 
  rename(SelectedForP = value)

SelectedForP_priv <- desc_91_houses %>%
  filter(sample != "NTR") %>% 
  filter(!is.na(Act)) %>% 
  filter(sample == "ATT") %>% 
  select(-Act) %>% 
  desctable() %>% 
  format_desctable() %>% 
  rename(SelectedForP_priv = value)


SelectedForP_notpriv <- desc_91_houses %>%
  filter(sample != "NTR") %>% 
  filter(!is.na(Act)) %>% 
  filter(sample != "ATT") %>% 
  select(-Act) %>% 
  desctable() %>% 
  format_desctable() %>% 
  rename(SelectedForP_notpriv = value)

NotSelected <- desc_91_houses %>% 
  filter(sample != "NTR") %>% 
  filter(is.na(Act)) %>% 
  select(-Act) %>% 
  desctable() %>% 
  format_desctable() %>% 
  rename(NotSelected = value)


desc_91_houses %>% 
  filter(sample != "NTR") %>% 
  mutate(
    sample = case_when(
      is.na(Act) ~ "NDp",
      TRUE ~ "Dp"
    )
  ) %>% 
  select(-Act) %>% 
  gather(variable,value,-sample) %>%
  mutate(
    sample = factor(sample) %>% relevel(ref = "Dp")
  ) %>% 
  #drop_na() %>% 
  mutate(value = as.double(value)) %>% 
  group_by(variable) %>%
  do(
    estmlm = lm(
      value ~ sample,
      data = .
    ) %>% 
      coeftest(vcov. = vcovHC) %>% tidy() %>% 
      filter(str_detect(term,"sample"))
  ) %>% 
  ungroup() %>% 
  unnest(estmlm) %>% 
  rowwise() %>% 
  mutate(
    mean = estimate %>% format(digits = 1, nsmall = dig, trim = TRUE, scientific = FALSE) %>% 
      str_c(.,add_stars(p.value, latex = TRUE)),
    se = std.error %>% format(digits = 1, nsmall = dig, trim = TRUE, scientific = FALSE) %>%
      str_c("(",.,")")
  ) %>% 
  ungroup() %>% 
  select(variable,mean,se) %>% 
  gather(stat,dif2,-variable) -> dif2


FullSample <- desc_91_houses %>% 
  #filter(sample != "NTR") %>% 
  #filter(is.na(Act)) %>% 
  select(-Act) %>% 
  desctable() %>% 
  format_desctable() %>% 
  rename(FullSample = value)

house_tab91 <- list(
  #FullSample,
  Restituted,
  CityOwned,
  NotSelected,
  SelectedForP,
  SelectedForP_priv,
  SelectedForP_notpriv,
  dif1,
  dif2
) %>% 
  reduce(full_join) %>%
  mutate(
    variable = ifelse(stat == "n","obs",variable)
  ) %>% 
  filter(
    !(variable %in% c("EXP","E"))
  ) %>% 
  gather(sample,value,-variable,-stat) %>% 
  mutate(
    value = ifelse(stat == "n",as.character(as.integer(value)),value)
  ) %>% 
  spread(sample,value) %>% 
  mutate(
    variable = factor(
      variable,
      levels = c("EA","E","U","MALE","LVEK","EDUClow","EDUCmiddle","EDUChigh",
                                             "EXP","BB",
                                             "htypeComplete_family","htypeIncomplete_family",
                                             "htypeNon_family_household","htypeOne_person_household","POC_BYTU","obs"),
                       labels = c("Economically active (share)",
                                  "Employed (share)",
                                  "Unemployed (share)",
                                  "Male (share)", "Age",
                                  "Primary education (share)", 
                                  "Secondary education (share)", 
                                  "Tertiary education (share)",
                                  "Potential experience","Born in Brno (share)",
                                  "Household type: Complete family",
                                  "Household type: Incomplete family",
                                  "Household type: Non-family household",
                                  "Household type: One-person household",
                                  "Number of apartments",
                                  "Houses")
    )
  ) %>% 
  arrange(variable,stat) %>%
  mutate_all(as.character) %>% 
  mutate(
    variable = ifelse(stat == "se","",variable)
  ) %>% 
  select(
    variable,
  #FullSample,
  Restituted,
  CityOwned,
  NotSelected,
  SelectedForP,
  SelectedForP_priv,
  SelectedForP_notpriv,
  dif1,
  dif2
  )

house_tab91 %>% 
  mutate_all(replace_na,"") %>% 
  mutate_all(str_remove_all,"\\\\c") %>% 
      mutate_all(str_remove_all,"\\{") %>% 
      mutate_all(str_remove_all,"\\}") %>% 
      mutate_all(str_remove_all,"\\^") %>% 
  kable(
    col.names = c(
      "Variable",
      str_c("(",1:8,")")
    )
  )
```

| Variable                    | \(1\)   | \(2\)   | \(3\)   | \(4\)   | \(5\)   | \(6\)   | \(7\)       | \(8\)        |
|:----------------------------|:--------|:--------|:--------|:--------|:--------|:--------|:------------|:-------------|
| Economically active (share) | 0.526   | 0.506   | 0.499   | 0.562   | 0.580   | 0.553   | 0.021\*\*\* | -0.063\*\*\* |
|                             | (0.006) | (0.004) | (0.005) | (0.011) | (0.020) | (0.014) | (0.008)     | (0.012)      |
| Unemployed (share)          | 0.027   | 0.023   | 0.023   | 0.024   | 0.022   | 0.025   | 0.004\*     | -0.001       |
|                             | (0.002) | (0.001) | (0.001) | (0.003) | (0.005) | (0.004) | (0.002)     | (0.003)      |
| Male (share)                | 0.443   | 0.442   | 0.441   | 0.451   | 0.461   | 0.446   | 0.0008      | -0.010       |
|                             | (0.003) | (0.002) | (0.002) | (0.006) | (0.011) | (0.007) | (0.004)     | (0.006)      |
| Age                         | 49.590  | 50.036  | 50.290  | 47.874  | 47.288  | 48.149  | -0.446      | 2.416\*\*\*  |
|                             | (0.257) | (0.172) | (0.184) | (0.432) | (0.760) | (0.527) | (0.309)     | (0.469)      |
| Primary education (share)   | 0.556   | 0.573   | 0.586   | 0.466   | 0.444   | 0.476   | -0.017\*    | 0.120\*\*\*  |
|                             | (0.007) | (0.005) | (0.005) | (0.015) | (0.026) | (0.018) | (0.009)     | (0.016)      |
| Secondary education (share) | 0.286   | 0.282   | 0.278   | 0.319   | 0.317   | 0.321   | 0.003       | -0.041\*\*\* |
|                             | (0.005) | (0.003) | (0.004) | (0.010) | (0.015) | (0.012) | (0.006)     | (0.010)      |
| Tertiary education (share)  | 0.159   | 0.145   | 0.136   | 0.215   | 0.239   | 0.203   | 0.014\*\*   | -0.078\*\*\* |
|                             | (0.005) | (0.003) | (0.003) | (0.011) | (0.021) | (0.012) | (0.006)     | (0.011)      |
| Houses                      | 605     | 1371    | 1227    | 144     | 46      | 98      |             |              |
|                             | 0.534   | 0.499   | 0.495   | 0.528   | 0.558   | 0.513   | 0.036\*\*\* | -0.033\*\*\* |
|                             | 11.398  | 12.002  | 11.944  | 12.500  | 11.435  | 13.000  | -0.604\*\*  | -0.556       |
|                             | 13.140  | 13.917  | 13.560  | 16.958  | 16.304  | 17.265  | -0.776\*\*  | -3.398\*\*\* |
|                             | (0.006) | (0.004) | (0.004) | (0.011) | (0.014) | (0.014) | (0.007)     | (0.012)      |
|                             | (0.199) | (0.134) | (0.143) | (0.397) | (0.561) | (0.517) | (0.240)     | (0.422)      |
|                             | (0.279) | (0.186) | (0.194) | (0.563) | (0.951) | (0.701) | (0.335)     | (0.596)      |

#### Table 2: Working age (18 to 60 years) individuals in Brno in 2001, means by privatization and restitution status of the house

``` r
r01_data_basic %>% 
  select(one_of(all.vars(model[[5]])),adresado_01_klic,flat_ID,sample,priv_year,Act,LVEK) %>% 
  drop_na(-priv_year,-Act) %>% 
  mutate(
    U = as.integer(LEKAKTI == 5),
    #E = as.integer(LEKAKTI == 1),
    EA = as.integer(LEKAKTI %in% c(1,5)),
    owner = as.integer(ownership == "owner"),
    EDUClow = as.integer(EDUC == "low"),
    EDUCmiddle = as.integer(EDUC == "middle"),
    EDUChigh = as.integer(EDUC == "high"),
    MALE = as.integer(LPOHLAV == "male"),
    htypeIncomplete_family = as.integer(htype == "Incomplete_family"),
    htypeComplete_family = as.integer(htype == "Complete_family"),
    htypeOne_person_household = as.integer(htype == "One_person_household")
  ) %>% 
  #select(-EXP) %>% 
  select(-LEKAKTI,-ownership,-EDUC,-LPOHLAV,-RZSJ,-htype) -> desc01_ind


## Restituted

Ntab <- desc01_ind %>% 
  filter(sample == "NTR") %>% 
  summarise(
    individuals = n(),
    houses = unique(adresado_01_klic) %>% length(),
    households = unique(flat_ID) %>% length()
  ) %>% mutate_all(as.character) %>% 
  gather(variable,value) %>% 
  mutate(stat = "n")

NTRstat <- desc01_ind %>% 
  filter(sample == "NTR") %>% 
  select(-Act,-adresado_01_klic,-flat_ID) %>% 
  desctable() %>%
  format_desctable() %>% 
  bind_rows(.,Ntab) %>% 
  rename(NTR = value)

## City-owned

Ntab <- desc01_ind %>% 
  filter(sample != "NTR") %>% 
  summarise(
    individuals = n(),
    houses = unique(adresado_01_klic) %>% length(),
    households = unique(flat_ID) %>% length()
  ) %>% mutate_all(as.character) %>% 
  gather(variable,value) %>% 
  mutate(stat = "n")

COstat <- desc01_ind %>% 
  filter(sample != "NTR") %>% 
  select(-Act,-adresado_01_klic,-flat_ID) %>% 
  desctable() %>%
  format_desctable() %>% 
  bind_rows(.,Ntab) %>% 
  rename(CO = value)

## Privatized
Ntab <- desc01_ind %>% 
  filter(sample == "ATT") %>% 
  summarise(
    individuals = n(),
    houses = unique(adresado_01_klic) %>% length(),
    households = unique(flat_ID) %>% length()
  ) %>% mutate_all(as.character) %>% 
  gather(variable,value) %>% 
  mutate(stat = "n")

ATTstat <- desc01_ind %>% 
  filter(sample == "ATT") %>% 
  select(-Act,-adresado_01_klic,-flat_ID) %>% 
  desctable() %>%
  format_desctable() %>% 
  bind_rows(.,Ntab) %>% 
  rename(ATT = value)

## City-owned/NTP

Ntab <- desc01_ind %>% 
  filter(sample == "NTP") %>% 
  summarise(
    individuals = n(),
    houses = unique(adresado_01_klic) %>% length(),
    households = unique(flat_ID) %>% length()
  ) %>% mutate_all(as.character) %>% 
  gather(variable,value) %>% 
  mutate(stat = "n")

NTPstat <- desc01_ind %>% 
  filter(sample == "NTP") %>% 
  select(-Act,-adresado_01_klic,-flat_ID) %>% 
  desctable() %>%
  format_desctable() %>% 
  bind_rows(.,Ntab) %>% 
  rename(NTP = value)

## City-owned/NTC

Ntab <- desc01_ind %>% 
  filter(sample == "NTC") %>% 
  summarise(
    individuals = n(),
    houses = unique(adresado_01_klic) %>% length(),
    households = unique(flat_ID) %>% length()
  ) %>% mutate_all(as.character) %>% 
  gather(variable,value) %>% 
  mutate(stat = "n")

NTCstat <- desc01_ind %>% 
  filter(sample == "NTC") %>% 
  select(-Act,-adresado_01_klic,-flat_ID) %>% 
  desctable() %>%
  format_desctable() %>% 
  bind_rows(.,Ntab) %>% 
  rename(NTC = value)


### dif1

desc01_ind %>% 
  filter(sample %in% c("NTR","ATT")) %>% 
  select(-Act,-adresado_01_klic,-flat_ID) %>% 
  select(-priv_year) %>% 
  gather(variable,value,-sample) %>% 
  mutate(
    sample = factor(sample) %>% relevel(ref = "ATT")
  ) %>% 
  group_by(variable) %>% 
  do(
    estmlm = lm(
      value ~ sample,
      data = .
    ) %>% 
      coeftest(vcov. = vcovHC) %>% tidy() %>% 
      filter(str_detect(term,"sample"))
  ) %>% 
  ungroup() %>% 
  unnest(estmlm) %>% 
  rowwise() %>% 
  mutate(
    mean = estimate %>% format(digits = 1, nsmall = 3, trim = TRUE, scientific = FALSE) %>% 
      str_c(.,add_stars(p.value, latex = TRUE)),
    se = std.error %>% format(digits = 1, nsmall = 3, trim = TRUE, scientific = FALSE) %>%
      str_c("(",.,")")
  ) %>% 
  ungroup() %>% 
  select(variable,mean,se) %>% 
  gather(stat,dif1,-variable) -> dif1
  


estsamp2001 <- list(
  NTRstat,
  ATTstat,
  NTPstat,
  NTCstat,
  COstat,
  dif1
) %>% 
  reduce(full_join) %>%
  # mutate(
  #   variable = ifelse(stat == "n","obs",variable)
  # ) %>% 
  filter(
    !(variable %in% c("EXP","E","Observations"))
  ) %>% 
  gather(
  sample,value,-variable,-stat
  ) %>% 
  mutate(
    value = ifelse(stat == "n", format(as.integer(value), big.mark = ",", trim = TRUE), value)
  ) %>% 
  spread(sample,value) %>% 
  mutate(
    variable = factor(
      variable,
      levels = c("EA","E","U","MALE","LVEK","EDUClow","EDUCmiddle","EDUChigh",
                                             "EXP","BornBrno",
                                             "htypeComplete_family","htypeIncomplete_family",
                                             "htypeNon_family_household","htypeOne_person_household","owner",
                 "individuals","households","houses"),
                       labels = c("Economically active (=1)",
                                  "Employed (=1)",
                                  "Unemployed (=1)",
                                  "Male (=1)", "Age",
                                  "Primary education (=1)", 
                                  "Secondary education (=1)",
                                  "Tertiary education (=1)",
                                  "Potential experience","Born in Brno (=1)",
                                  "Household type: Complete family (=1)",
                                  "Household type: Incomplete family (=1)",
                                  "Household type: Non-family household (=1)",
                                  "Household type: One-person household (=1)",
                                  "Owner (=1)",
                                  "Individuals",
                                  "Households",
                                  "Houses")
    )
  ) %>% 
  arrange(variable,stat) %>% 
  mutate_all(as.character) %>% 
  mutate(
    variable = ifelse(stat == "se","",variable)
  ) %>% 
  select(
    variable,
    NTR,
  ATT,
  NTP,
  NTC,
  CO,
  dif1
  )

estsamp2001  %>% 
  mutate_all(as.character) %>% 
  mutate_all(replace_na,"") %>% 
  mutate_all(str_remove_all,"\\\\c") %>% 
      mutate_all(str_remove_all,"\\{") %>% 
      mutate_all(str_remove_all,"\\}") %>% 
      mutate_all(str_remove_all,"\\^") %>% 
  kable(
    col.names = c(
      "Variable",
      str_c("(",1:6,")")
    )
  )
```

| Variable                                  | \(1\)   | \(2\)   | \(3\)   | \(4\)   | \(5\)   | \(6\)        |
|:------------------------------------------|:--------|:--------|:--------|:--------|:--------|:-------------|
| Economically active (=1)                  | 0.786   | 0.787   | 0.795   | 0.786   | 0.792   | -0.001       |
|                                           | (0.005) | (0.017) | (0.004) | (0.005) | (0.003) | (0.018)      |
| Unemployed (=1)                           | 0.099   | 0.060   | 0.097   | 0.116   | 0.103   | 0.039\*\*\*  |
|                                           | (0.004) | (0.010) | (0.003) | (0.004) | (0.002) | (0.010)      |
| Male (=1)                                 | 0.487   | 0.494   | 0.475   | 0.484   | 0.479   | -0.007       |
|                                           | (0.006) | (0.020) | (0.005) | (0.006) | (0.004) | (0.021)      |
| Age                                       | 39.453  | 39.476  | 39.109  | 39.032  | 39.092  | -0.023       |
|                                           | (0.155) | (0.492) | (0.122) | (0.146) | (0.092) | (0.516)      |
| Primary education (=1)                    | 0.409   | 0.307   | 0.399   | 0.473   | 0.425   | 0.102\*\*\*  |
|                                           | (0.006) | (0.019) | (0.005) | (0.006) | (0.004) | (0.020)      |
| Secondary education (=1)                  | 0.367   | 0.365   | 0.378   | 0.341   | 0.363   | 0.002        |
|                                           | (0.006) | (0.020) | (0.005) | (0.006) | (0.004) | (0.021)      |
| Tertiary education (=1)                   | 0.225   | 0.328   | 0.223   | 0.187   | 0.212   | -0.104\*\*\* |
|                                           | (0.005) | (0.019) | (0.004) | (0.005) | (0.003) | (0.020)      |
| Born in Brno (=1)                         | 0.646   | 0.668   | 0.633   | 0.634   | 0.635   | -0.023       |
|                                           | (0.006) | (0.019) | (0.005) | (0.006) | (0.004) | (0.020)      |
| Household type: Complete family (=1)      | 0.628   | 0.690   | 0.622   | 0.618   | 0.623   | -0.062\*\*\* |
|                                           | (0.006) | (0.019) | (0.005) | (0.006) | (0.004) | (0.020)      |
| Household type: Incomplete family (=1)    | 0.205   | 0.206   | 0.222   | 0.215   | 0.219   | -0.001       |
|                                           | (0.005) | (0.017) | (0.004) | (0.005) | (0.003) | (0.017)      |
| Household type: One-person household (=1) | 0.167   | 0.104   | 0.155   | 0.167   | 0.158   | 0.063\*\*\*  |
|                                           | (0.005) | (0.013) | (0.004) | (0.005) | (0.003) | (0.013)      |
| Owner (=1)                                | 0.000   | 0.913   | 0.000   | 0.000   | 0.033   | -0.913\*\*\* |
|                                           | (0.000) | (0.012) | (0.000) | (0.000) | (0.001) | (0.012)      |
| Individuals                               | 6,088   | 597     | 9,552   | 6,611   | 16,760  | NA           |
| Households                                | 3,556   | 329     | 5,725   | 4,007   | 10,061  | NA           |
| Houses                                    | 599     | 46      | 774     | 548     | 1,368   | NA           |

## Regressions

## OLS

### Table 3: Homeownership and labor market activity (OLS), individuals from houses privatized before 2001 and all restituted houses

``` r
model <- model_individual

# Remove incomplete observations
model %>% 
  lapply(all.vars) %>% 
  unlist() %>% 
  unique() %>% 
  c("flat_ID","sample","RZSJ","MALE","adresado_01_klic",.) -> all_vars

r01_data_x <- r01_data_basic %>%
  dplyr::select(one_of(all_vars)) %>% 
  drop_na()
```

``` r
r01_data_x %>% 
  filter(sample %in% c("ATT","NTR")) -> reg_data

reg_data_ea <- reg_data %>% filter(LEKAKTI %in% c(1,5))

#save(reg_data, reg_data_ea, file = "DATA/processed/fullNAIVE_data.RData")

c(
  model[1:5] %>% lapply(function(x) lm(x, data = reg_data_ea)),
  model[6:10] %>% lapply(function(x) lm(x, data = reg_data))
  ) -> estm

if(cluster_by == "flat"){
  c(
    estm[1:5] %>% lapply(function(x) coeftest(x, vcov. = vcovCL(x, reg_data_ea$flat_ID))),
    estm[6:10] %>% lapply(function(x) coeftest(x, vcov. = vcovCL(x, reg_data$flat_ID)))
  ) -> estm_ctest
}else{
  c(
    estm[1:5] %>% lapply(function(x) coeftest(x, vcov. = vcovCL(x, reg_data_ea$adresado_01_klic))),
    estm[6:10] %>% lapply(function(x) coeftest(x, vcov. = vcovCL(x, reg_data$adresado_01_klic)))
  ) -> estm_ctest
}

get_tables(estm, estm_ctest, reorder = c(6,7,8,9,10,1,2,3,4,5), space = FALSE) %>% 
  mutate_all(replace_na,"") %>% 
  kable(
    col.names = c(
      "Variable",
      str_c("(",1:10,")")
    )
  )
```

| Variable                                                                                    | \(1\)       | \(2\)   | \(3\)       | \(4\)       | \(5\)       | \(6\)        | \(7\)        | \(8\)        | \(9\)        | \(10\)       |
|:--------------------------------------------------------------------------------------------|:------------|:--------|:------------|:------------|:------------|:-------------|:-------------|:-------------|:-------------|:-------------|
| Homeowner (=1)                                                                              | 0.009       | 0.002   | -0.007      | -0.002      | -0.005      | -0.060\*\*\* | -0.052\*\*\* | -0.038\*\*\* | -0.036\*\*\* | -0.038\*\*\* |
|                                                                                             | (0.015)     | (0.016) | (0.015)     | (0.014)     | (0.016)     | (0.014)      | (0.014)      | (0.014)      | (0.013)      | (0.013)      |
| Secondary education (=1)                                                                    |             |         | -0.005      | 0.026\*\*   | 0.027\*\*   |              |              | -0.123\*\*\* | -0.123\*\*\* | -0.119\*\*\* |
|                                                                                             |             |         | (0.012)     | (0.011)     | (0.011)     |              |              | (0.012)      | (0.012)      | (0.012)      |
| Tertiary education (=1)                                                                     |             |         | 0.133\*\*\* | 0.095\*\*\* | 0.097\*\*\* |              |              | -0.178\*\*\* | -0.157\*\*\* | -0.153\*\*\* |
|                                                                                             |             |         | (0.012)     | (0.011)     | (0.011)     |              |              | (0.011)      | (0.011)      | (0.011)      |
| Male (=1)                                                                                   |             |         |             | 0.157\*\*\* | 0.157\*\*\* |              |              |              | 0.001        | 0.0001       |
|                                                                                             |             |         |             | (0.009)     | (0.009)     |              |              |              | (0.008)      | (0.008)      |
| Born in Brno (=1)                                                                           |             |         |             | -0.004      | -0.005      |              |              |              | 0.004        | 0.007        |
|                                                                                             |             |         |             | (0.009)     | (0.009)     |              |              |              | (0.010)      | (0.010)      |
| Incomplete family (=1)                                                                      |             |         |             | 0.028\*\*   | 0.028\*\*   |              |              |              | 0.066\*\*\*  | 0.063\*\*\*  |
|                                                                                             |             |         |             | (0.012)     | (0.012)     |              |              |              | (0.013)      | (0.013)      |
| One-person family (=1)                                                                      |             |         |             | 0.015       | 0.014       |              |              |              | 0.051\*\*\*  | 0.051\*\*\*  |
|                                                                                             |             |         |             | (0.012)     | (0.013)     |              |              |              | (0.013)      | (0.013)      |
| Constant                                                                                    | 0.786\*\*\* |         | 0.758\*\*\* |             |             | 0.127\*\*\*  |              | 0.215\*\*\*  |              |              |
|                                                                                             | (0.005)     |         | (0.009)     |             |             | (0.006)      |              | (0.011)      |              |              |
| Age dummies                                                                                 | –           | –       | –           | Yes         | Yes         | –            | –            | –            | Yes          | Yes          |
| Neighborhood dummies                                                                        | –           | Yes     | –           | –           | Yes         | –            | Yes          | –            | –            | Yes          |
| Observations                                                                                | 6,685       | 6,685   | 6,685       | 6,685       | 6,685       | 5,256        | 5,256        | 5,256        | 5,256        | 5,256        |
| ![R2](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;R2 "R2") | 0.00004     | 0.009   | 0.020       | 0.219       | 0.224       | 0.003        | 0.021        | 0.054        | 0.092        | 0.102        |

## Matching

``` r
#Add RZSJ coordinates.
houses <- r01_data %>% 
  group_by(adresado_01_klic) %>%
  dplyr::summarise(
    RCASTOBCE = first(RCASTOBCE),
    RZSJ = first(RZSJ)#,
    #JPOCBYT = first(JPOCBYT)
  ) %>% 
  mutate(
    RZSJ = str_sub(RZSJ, start = 1, end = 6)
  ) %>% 
  ungroup() %>% 
  distinct()

load("DATA_RP/rzsj_centroides.RData")

rzsj %>%
  as.data.frame() %>% 
  as_data_frame() %>% 
  dplyr::select(
    RZSJ = KOD_ZSJ, long = coords.x1, lat = coords.x2
  ) %>% 
  mutate(
    RZSJ = as.character(RZSJ)
  ) %>% 
  distinct() %>% 
  left_join(houses,., by = "RZSJ") %>% 
  drop_na() -> houses

desc_91_houses_klic %>% 
  left_join(houses) -> r91_houses
```

### Table A1: Balancing test (means), matched houses in 1991

``` r
# Propensity score
r91_houses %>%
  dplyr::select(-Act) %>%
  filter(sample %in% c("ATT","NTR")) %>%
  filter(adresado_01_klic %in% r01_data_basic$adresado_01_klic) %>%
  #filter(JPOCBYT >= 6) %>%
  rename(JPOCBYT = POC_BYTU) %>%
  mutate(
    byty = cut(JPOCBYT, c(0,10,20,40), right = FALSE)
  ) %>%
  drop_na() %>%
  rename(BornBrno = BB) -> r91_houses_sample1


att_selection <- glm(update(model_propscore, I(sample == "ATT") ~ .), 
                     data = r91_houses_sample1, family = binomial("logit"))

r91_houses_sample1 %>% 
  predict.glm(att_selection, type = "response", newdata = .) -> r91_houses_sample1$pscore


### Regression data
modeliv <- modeliv_individual #%>% lapply(matchingFE)
model_rob <- modeliv[[5]]

# Remove incomplete observations
model_rob %>% 
  lapply(all.vars) %>% 
  unlist() %>% 
  unique() %>% 
  c("flat_ID","sample","RZSJ","MALE","adresado_01_klic","person_ID",.) -> all_vars

r01_data_x <- r01_data_basic %>%
  dplyr::select(one_of(all_vars)) %>% 
  drop_na()

r01_data_x %>% 
  filter(sample %in% c("ATT","NTR")) %>% 
  mutate(
    ownership = ownership == "owner",
    sample_dummy = sample == "ATT"
    ) -> reg_data

reg_data_ea <- reg_data %>% filter(LEKAKTI %in% c(1,5))

## Unemployment

load("DATA_RP/zsjsf.RData")

crossing(
  C = seq(from=0.1, to = 0.9, by = 0.1),
  D = c(100,400,700,1000)
) %>% 
  #filter(D %in% c(100,400,700,1000)) %>%
  # filter(D == 700) %>% 
  # filter(C == 0.4) %>% 
  rowwise() %>% 
  mutate(
    match = get_matched.geo(r91_houses_sample1, C, D, geo = zsjsf) %>% list(),
    #match = matchO$DataX %>% list(),
    emod  = get_ownership.g0(match$match, xdata = reg_data_ea, xmodel = model_rob) %>% list(),
    owner = emod$coefs %>% filter(str_detect(term,"ownershipTRUE")) %>% list()
  ) -> outputs_U
```

    ## [1] "Cluster Robust Standard Errors"
    ## [1] "Cluster Robust Standard Errors"
    ## [1] "Cluster Robust Standard Errors"
    ## [1] "Cluster Robust Standard Errors"
    ## [1] "Cluster Robust Standard Errors"
    ## [1] "Cluster Robust Standard Errors"
    ## [1] "Cluster Robust Standard Errors"
    ## [1] "Cluster Robust Standard Errors"
    ## [1] "Cluster Robust Standard Errors"
    ## [1] "Cluster Robust Standard Errors"
    ## [1] "Cluster Robust Standard Errors"
    ## [1] "Cluster Robust Standard Errors"
    ## [1] "Cluster Robust Standard Errors"
    ## [1] "Cluster Robust Standard Errors"
    ## [1] "Cluster Robust Standard Errors"
    ## [1] "Cluster Robust Standard Errors"
    ## [1] "Cluster Robust Standard Errors"
    ## [1] "Cluster Robust Standard Errors"
    ## [1] "Cluster Robust Standard Errors"
    ## [1] "Cluster Robust Standard Errors"
    ## [1] "Cluster Robust Standard Errors"
    ## [1] "Cluster Robust Standard Errors"
    ## [1] "Cluster Robust Standard Errors"
    ## [1] "Cluster Robust Standard Errors"
    ## [1] "Cluster Robust Standard Errors"
    ## [1] "Cluster Robust Standard Errors"
    ## [1] "Cluster Robust Standard Errors"
    ## [1] "Cluster Robust Standard Errors"
    ## [1] "Cluster Robust Standard Errors"
    ## [1] "Cluster Robust Standard Errors"
    ## [1] "Cluster Robust Standard Errors"
    ## [1] "Cluster Robust Standard Errors"
    ## [1] "Cluster Robust Standard Errors"
    ## [1] "Cluster Robust Standard Errors"
    ## [1] "Cluster Robust Standard Errors"
    ## [1] "Cluster Robust Standard Errors"

``` r
outputs_U %>% 
  unnest(owner) %>% 
  rowwise() %>% 
  mutate(
    n_t = match$match %>% filter(ctg == "treatment") %>% distinct(adresado_01_klic) %>% nrow(),
    n_c = match$match %>% filter(ctg == "control") %>% distinct(adresado_01_klic) %>% nrow()
  ) -> outputs_U


# Balance test
outputs_U %>% 
  #filter(C==d_caliper,D==d_distance) %>%
  filter(is_equal(C,d_caliper)) %>% 
  filter(is_equal(D,d_distance)) %>% 
  pull(match) -> baseline_match

aux <- baseline_match[[1]]$match %>% 
  left_join(r91_houses_sample1, by = "adresado_01_klic") %>% 
  select(-adresado_01_klic,-sample,-RZSJ,-RCASTOBCE,-long,-lat,-byty) %>% 
  gather(variable,value,-ctg,-id,-W)

aux1 <- aux %>% 
  group_by(ctg,variable) %>% 
  do(
    mean = lm(value ~ 1, data = ., weights = W) %>% coeftest(vcov. = vcovHC) %>% tidy() %>% filter(str_detect(term,"Intercept"))
  ) %>% 
  ungroup() %>% 
  unnest(mean) %>%
  rowwise() %>% 
  mutate(
    estimate = format(estimate, digits = 1, nsmall=3, trim = TRUE, scientific = FALSE),
    std.error = format(std.error, digits = 1, nsmall=3, trim = TRUE, scientific = FALSE) %>% str_c("(",.,")")
  ) %>% 
  ungroup() %>% 
  select(ctg,variable,estimate,std.error) %>% 
  gather(stat,value,-ctg,-variable) %>% 
  spread(ctg,value)

aux %>% 
  mutate(
    ctg = factor(ctg) %>% relevel(ref = "control")
  ) %>% 
  group_by(variable) %>% 
  do(
    dif = lm(value ~ ctg, data = ., weights = W) %>% coeftest(vcov. = vcovHC) %>% tidy() %>% filter(str_detect(term,"ctg"))
  ) %>% 
  ungroup() %>% 
  unnest(dif) %>%
  rowwise() %>% 
  mutate(
    estimate = format(estimate, digits = 1, nsmall = 3, trim = TRUE, scientific = FALSE) %>% str_c(.,add_stars(p.value, latex = TRUE)),
    std.error = format(std.error, digits = 1, nsmall = 3, trim = TRUE, scientific = FALSE) %>% str_c("(",.,")")
  ) %>% 
  select(variable,estimate,std.error) %>%
  gather(stat,dif,-variable) %>% 
  left_join(aux1,.) %>% 
  filter(!(variable %in% c("EXP","E","POC_BYTU","W"))) %>% 
  #distinct(variable)
  mutate(
    variable = factor(
      variable,
      levels = c("EA","E","U","MALE","LVEK","EDUClow","EDUCmiddle","EDUChigh",
                                             "EXP","BornBrno",
                                             "htypeComplete_family","htypeIncomplete_family",
                                      "htypeNon_family_household","htypeOne_person_household","owner","JPOCBYT","pscore","obs"),
                       labels = c("Economically active (share)",
                                  "Employed (share)",
                                  "Unemployed (share)",
                                  "Male (share)", "Age",
                                  "Primary education (share)", 
                                  "Secondary education (share)", 
                                  "Tertiary education (share)",
                                  "Potential experience","Born in Brno (share)",
                                  "Household type: Complete family (=1)",
                                  "Household type: Incomplete family (=1)",
                                  "Household type: Non-family household (=1)",
                                  "Household type: One-person household (=1)",
                                  "Owner (=1)",
                                  "Number of apartments",
                                  "Propensity score",
                                  "Observations (houses)")
    )
  ) %>% 
  select(variable,treatment,control,dif,stat) %>% 
  arrange(variable,stat) %>% 
  mutate_all(as.character) %>% 
  mutate(
    variable = ifelse(stat == "std.error","",variable)
    #dif = ifelse(stat == "std.error","",dif)
  ) %>% 
  select(-stat) -> match_balance
 
tribble(
  ~variable,~treatment,~control,~dif,
  "Houses",
  str_c(
    as.character(length(unique(baseline_match[[1]]$match_obj$index.treated))) %>% str_c("\\c{",.,"}")
    # "/",
    #as.character(length(baseline_match[[1]]$match_obj$index.treated))
  ),
  str_c(
    as.character(length(unique(baseline_match[[1]]$match_obj$index.control))) %>% str_c("\\c{",.,"}")
    # "/",
    #as.character(length(baseline_match[[1]]$match_obj$index.control))
  ),
  ""
) %>% 
bind_rows(match_balance,.) %>% 
   mutate_all(as.character) %>% 
  mutate_all(replace_na,"") %>% 
  mutate_all(str_remove_all,"\\\\c") %>% 
      mutate_all(str_remove_all,"\\{") %>% 
      mutate_all(str_remove_all,"\\}") %>% 
      mutate_all(str_remove_all,"\\^") %>% 
  kable(
    col.names = c(
      "Variable",
      str_c("(",1:3,")")
    )
  )
```

| Variable                    | \(1\)   | \(2\)   | \(3\)   |
|:----------------------------|:--------|:--------|:--------|
| Economically active (share) | 0.585   | 0.596   | -0.011  |
|                             | (0.020) | (0.018) | (0.027) |
| Unemployed (share)          | 0.021   | 0.020   | 0.001   |
|                             | (0.005) | (0.004) | (0.006) |
| Male (share)                | 0.461   | 0.448   | 0.012   |
|                             | (0.010) | (0.011) | (0.015) |
| Age                         | 47.220  | 47.494  | -0.273  |
|                             | (0.784) | (0.728) | (1.070) |
| Primary education (share)   | 0.443   | 0.438   | 0.005   |
|                             | (0.027) | (0.029) | (0.039) |
| Secondary education (share) | 0.313   | 0.329   | -0.015  |
|                             | (0.015) | (0.018) | (0.024) |
| Tertiary education (share)  | 0.244   | 0.233   | 0.010   |
|                             | (0.021) | (0.021) | (0.029) |
| Born in Brno (share)        | 0.557   | 0.551   | 0.006   |
|                             | (0.015) | (0.021) | (0.026) |
| Number of apartments        | 11.564  | 12.949  | -1.385  |
|                             | (0.566) | (0.837) | (1.010) |
| Propensity score            | 0.114   | 0.112   | 0.002   |
|                             | (0.011) | (0.010) | (0.015) |
| Houses                      | 46      | 46      |         |

## Balance test W

``` r
aux <- baseline_match[[1]]$match %>% 
  left_join(r91_houses_sample1, by = "adresado_01_klic") %>% 
  distinct(adresado_01_klic, .keep_all = TRUE) %>% 
  select(-adresado_01_klic,-sample,-RZSJ,-RCASTOBCE,-long,-lat,-byty) %>% 
  gather(variable,value,-ctg,-id,-W)

aux1 <- aux %>% 
  group_by(ctg,variable) %>% 
  do(
    mean = lm(value ~ 1, data = .) %>% coeftest(vcov. = vcovHC) %>% tidy() %>% filter(str_detect(term,"Intercept"))
  ) %>% 
  ungroup() %>% 
  unnest(mean) %>%
  rowwise() %>% 
  mutate(
    estimate = format(estimate, digits = 1, nsmall=3, trim = TRUE, scientific = FALSE),
    std.error = format(std.error, digits = 1, nsmall=3, trim = TRUE, scientific = FALSE) %>% str_c("(",.,")")
  ) %>% 
  ungroup() %>% 
  select(ctg,variable,estimate,std.error) %>% 
  gather(stat,value,-ctg,-variable) %>% 
  spread(ctg,value)

aux %>% 
  mutate(
    ctg = factor(ctg) %>% relevel(ref = "control")
  ) %>% 
  group_by(variable) %>% 
  do(
    dif = lm(value ~ ctg, data = .) %>% coeftest(vcov. = vcovHC) %>% tidy() %>% filter(str_detect(term,"ctg"))
  ) %>% 
  ungroup() %>% 
  unnest(dif) %>%
  rowwise() %>% 
  mutate(
    estimate = format(estimate, digits = 1, nsmall = 3, trim = TRUE, scientific = FALSE) %>% str_c(.,add_stars(p.value, latex = TRUE)),
    std.error = format(std.error, digits = 1, nsmall = 3, trim = TRUE, scientific = FALSE) %>% str_c("(",.,")")
  ) %>% 
  select(variable,estimate,std.error) %>%
  gather(stat,dif,-variable) %>% 
  left_join(aux1,.) %>% 
  filter(!(variable %in% c("EXP","E","POC_BYTU","W"))) %>% 
  #distinct(variable)
  mutate(
    variable = factor(
      variable,
      levels = c("EA","E","U","MALE","LVEK","EDUClow","EDUCmiddle","EDUChigh",
                                             "EXP","BornBrno",
                                             "htypeComplete_family","htypeIncomplete_family",
                                      "htypeNon_family_household","htypeOne_person_household","owner","JPOCBYT","pscore","obs"),
                       labels = c("Economically active (share)",
                                  "Employed (share)",
                                  "Unemployed (share)",
                                  "Male (share)", "Age",
                                  "Primary education (share)", 
                                  "Secondary education (share)", 
                                  "Tertiary education (share)",
                                  "Potential experience","Born in Brno (share)",
                                  "Household type: Complete family (=1)",
                                  "Household type: Incomplete family (=1)",
                                  "Household type: Non-family household (=1)",
                                  "Household type: One-person household (=1)",
                                  "Owner (=1)",
                                  "Number of apartments",
                                  "Propensity score",
                                  "Observations (houses)")
    )
  ) %>% 
  select(variable,treatment,control,dif,stat) %>% 
  arrange(variable,stat) %>% 
  mutate_all(as.character) %>% 
  mutate(
    variable = ifelse(stat == "std.error","",variable)
    #dif = ifelse(stat == "std.error","",dif)
  ) %>% 
  select(-stat) -> match_balance
 
tribble(
  ~variable,~treatment,~control,~dif,
  "Houses",
  str_c(
    as.character(length(unique(baseline_match[[1]]$match_obj$index.treated))) %>% str_c("\\c{",.,"}")
    # "/",
    #as.character(length(baseline_match[[1]]$match_obj$index.treated))
  ),
  str_c(
    as.character(length(unique(baseline_match[[1]]$match_obj$index.control))) %>% str_c("\\c{",.,"}")
    # "/",
    #as.character(length(baseline_match[[1]]$match_obj$index.control))
  ),
  ""
) %>% 
bind_rows(match_balance,.) %>% 
   mutate_all(as.character) %>% 
  mutate_all(replace_na,"") %>% 
  mutate_all(str_remove_all,"\\\\c") %>% 
      mutate_all(str_remove_all,"\\{") %>% 
      mutate_all(str_remove_all,"\\}") %>% 
      mutate_all(str_remove_all,"\\^") -> Tab3print

Tab3print %>% print_latex("match_balance.tex")

Tab3print %>% 
  kable(
    col.names = c(
      "Variable",
      str_c("(",1:3,")")
    )
  )
```

| Variable                    | \(1\)   | \(2\)   | \(3\)   |
|:----------------------------|:--------|:--------|:--------|
| Economically active (share) | 0.580   | 0.600   | -0.020  |
|                             | (0.020) | (0.019) | (0.028) |
| Unemployed (share)          | 0.022   | 0.021   | 0.0004  |
|                             | (0.005) | (0.004) | (0.006) |
| Male (share)                | 0.461   | 0.448   | 0.012   |
|                             | (0.011) | (0.011) | (0.015) |
| Age                         | 47.288  | 47.529  | -0.240  |
|                             | (0.760) | (0.767) | (1.080) |
| Primary education (share)   | 0.444   | 0.459   | -0.015  |
|                             | (0.026) | (0.030) | (0.039) |
| Secondary education (share) | 0.317   | 0.322   | -0.005  |
|                             | (0.015) | (0.019) | (0.024) |
| Tertiary education (share)  | 0.239   | 0.219   | 0.020   |
|                             | (0.021) | (0.020) | (0.029) |
| Born in Brno (share)        | 0.558   | 0.551   | 0.007   |
|                             | (0.014) | (0.020) | (0.025) |
| Number of apartments        | 11.435  | 12.891  | -1.457  |
|                             | (0.561) | (0.846) | (1.015) |
| Propensity score            | 0.112   | 0.103   | 0.009   |
|                             | (0.011) | (0.010) | (0.015) |
| Houses                      | 46      | 46      |         |

``` r
## Economic activity
modeliv <- modeliv_individual #%>% lapply(matchingFE)
model_rob <- modeliv[[10]]

outputs_U %>% 
  select(C,D,match) %>% 
  rowwise() %>% 
  mutate(
    #match = get_matched.g0(r91_houses_sample1, C, D) %>% list(),
    emod  = get_ownership.g0(match$match, xdata = reg_data, xmodel = model_rob) %>% list(),
    owner = emod$coefs %>% filter(str_detect(term,"ownershipTRUE")) %>% list()
  ) -> outputs_EA
```

    ## [1] "Cluster Robust Standard Errors"
    ## [1] "Cluster Robust Standard Errors"
    ## [1] "Cluster Robust Standard Errors"
    ## [1] "Cluster Robust Standard Errors"
    ## [1] "Cluster Robust Standard Errors"
    ## [1] "Cluster Robust Standard Errors"
    ## [1] "Cluster Robust Standard Errors"
    ## [1] "Cluster Robust Standard Errors"
    ## [1] "Cluster Robust Standard Errors"
    ## [1] "Cluster Robust Standard Errors"
    ## [1] "Cluster Robust Standard Errors"
    ## [1] "Cluster Robust Standard Errors"
    ## [1] "Cluster Robust Standard Errors"
    ## [1] "Cluster Robust Standard Errors"
    ## [1] "Cluster Robust Standard Errors"
    ## [1] "Cluster Robust Standard Errors"
    ## [1] "Cluster Robust Standard Errors"
    ## [1] "Cluster Robust Standard Errors"
    ## [1] "Cluster Robust Standard Errors"
    ## [1] "Cluster Robust Standard Errors"
    ## [1] "Cluster Robust Standard Errors"
    ## [1] "Cluster Robust Standard Errors"
    ## [1] "Cluster Robust Standard Errors"
    ## [1] "Cluster Robust Standard Errors"
    ## [1] "Cluster Robust Standard Errors"
    ## [1] "Cluster Robust Standard Errors"
    ## [1] "Cluster Robust Standard Errors"
    ## [1] "Cluster Robust Standard Errors"
    ## [1] "Cluster Robust Standard Errors"
    ## [1] "Cluster Robust Standard Errors"
    ## [1] "Cluster Robust Standard Errors"
    ## [1] "Cluster Robust Standard Errors"
    ## [1] "Cluster Robust Standard Errors"
    ## [1] "Cluster Robust Standard Errors"
    ## [1] "Cluster Robust Standard Errors"
    ## [1] "Cluster Robust Standard Errors"

``` r
outputs_EA %>% 
  unnest(owner) %>% 
  rowwise() %>% 
  mutate(
    n_t = match$match %>% dplyr::filter(ctg == "treatment") %>% distinct(adresado_01_klic) %>% nrow(),
    n_c = match$match %>% dplyr::filter(ctg == "control") %>% distinct(adresado_01_klic) %>% nrow(),
  ) -> outputs_EA
```

### Table 8: Differences in characteristics of working age individuals living in matched restituted and privatized houses between 2001 and 1991

``` r
# Balance test
outputs_U %>% 
  #filter(C==d_caliper,D==d_distance) %>%
  filter(is_equal(C,d_caliper)) %>% 
  filter(is_equal(D,d_distance)) %>% 
  pull(match) -> baseline_match

### Balance in working population (1991)
r91_data_BT <- r91_data %>% 
  filter(LVEK < age_upper) %>% 
  left_join(baseline_match[[1]]$match,.) %>% 
  group_by(ctg) %>% 
  mutate(obs = n()) %>% 
  ungroup() %>% 
  select(W,ctg,one_of(all.vars(model_propscore)),EDUC,obs,adresado_01_klic) %>% 
  mutate(
    BornBrno = as.integer(BornBrno),
    EDUChigh = as.integer(EDUC == "high"),
    EDUClow = as.integer(EDUC == "low"),
    EDUCmiddle = as.integer(EDUC == "middle")
  ) %>% 
  select(-EDUC) %>% 
  gather(variable,value,-W,-ctg,-obs,-adresado_01_klic) %>% 
  mutate(
    census = "1991"
  )

r01_data_BT <- r01_data_raw %>% 
  filter(LVEK >= age_lower) %>% 
  filter(LVEK < age_upper) %>% 
 left_join(baseline_match[[1]]$match,.) %>% 
  group_by(ctg) %>% 
  mutate(obs = n()) %>% 
  ungroup() %>% 
  select(W,ctg,one_of(all.vars(model_propscore)),EDUC,LEKAKTI,obs,adresado_01_klic) %>% 
  mutate(
    U = as.integer(LEKAKTI == 5),
    EA = as.integer(LEKAKTI %in% c(1,5)),
    BornBrno = as.integer(BornBrno),
    EDUChigh = as.integer(EDUC == "high"),
    EDUClow = as.integer(EDUC == "low"),
    EDUCmiddle = as.integer(EDUC == "middle")
  ) %>% 
  select(-EDUC,-LEKAKTI,-JPOCBYT) %>% 
  gather(variable,value,-W,-ctg,-obs,-adresado_01_klic) %>% 
  mutate(
    census = "2001"
  )

TimeTravel <- bind_rows(
  r91_data_BT,
  r01_data_BT
) %>% 
  mutate(
    ctg = factor(ctg) %>% relevel(ref = "control"),
    census = factor(census) %>% relevel(ref = "1991")
  )

OBS <- TimeTravel %>% 
  group_by(ctg,census) %>% 
  summarise(
    obs = mean(obs) %>% as.integer() %>% format(big.mark = ",", trim = TRUE) %>% str_c("\\c{",.,"}")
  ) %>% 
  ungroup() %>% 
  spread(ctg,obs) %>% 
  mutate(
    variable = str_c("obs_",census),
    stat = "n"
  ) %>% 
  select(-census)

aux <- TimeTravel %>% 
  group_by(ctg,variable) %>% 
  do(
    lm = felm(value ~ census | 0 | 0 | adresado_01_klic, data = ., weights = .$W) %>% 
      coeftest() %>% 
      tidy() %>% 
      filter(str_detect(term,"census"))
  ) %>% 
  ungroup() %>% 
  unnest(lm) %>% 
  rowwise() %>% 
  mutate(
    estimate = format(estimate, digits = 1, nsmall = 3, trim=TRUE, scientific = FALSE) %>% 
      str_c(.,add_stars(p.value,latex = TRUE)),
    std.error = format(std.error, digits = 1, nsmall = 3, trim = TRUE, scientific = FALSE) %>% 
      str_c("(",.,")")
  ) %>% 
  ungroup() %>% 
  select(ctg,variable,estimate,std.error) %>% 
  gather(stat,value,-ctg,-variable) %>% 
  spread(ctg,value)
  
TimeTravel %>% 
  group_by(variable) %>% 
  do(
    lm = felm(value ~ ctg*census | 0 | 0 | adresado_01_klic, data = ., weights = .$W) %>% 
      coeftest() %>% 
      tidy() %>% 
      filter(str_detect(term,"treatment:census"))
  ) %>% 
  ungroup() %>% 
  unnest(lm) %>% 
  rowwise() %>% 
  mutate(
    estimate = format(estimate, digits = 1, nsmall = 3, trim=TRUE, scientific = FALSE) %>% 
      str_c(.,add_stars(p.value,latex = TRUE)),
    std.error = format(std.error, digits = 1, nsmall = 3, trim = TRUE, scientific = FALSE) %>% 
      str_c("(",.,")")
  ) %>% 
  select(variable,estimate,std.error) %>% 
  gather(stat,DID,-variable) %>% 
  left_join(aux,.) %>% 
  bind_rows(OBS) %>% 
  mutate(
    variable = factor(
      variable,
      levels = c("EA","E","U","MALE","LVEK","EDUClow","EDUCmiddle","EDUChigh",
                                             "EXP","BornBrno",
                                             "htypeComplete_family","htypeIncomplete_family",
                                      "htypeNon_family_household","htypeOne_person_household","owner","JPOCBYT","pscore",
                 "obs_2001","obs_1991"),
                       labels = c("Economically active",
                                  "Employed",
                                  "Unemployed",
                                  "Male", "Age",
                                  "Primary education", 
                                  "Secondary education", 
                                  "Tertiary education",
                                  "Potential experience","Born in Brno",
                                  "Household type: Complete family",
                                  "Household type: Incomplete family",
                                  "Household type: Non-family household",
                                  "Household type: One-person household",
                                  "Owner",
                                  "Number of apartments",
                                  "Propensity score",
                                  "Individuals in 2001 census",
                                  "Individuals in 1991 census")
    )
  ) %>% 
  #select(variable,treatment,control,dif,stat) %>% 
  arrange(variable,stat) %>% 
  mutate_all(as.character) %>% 
  mutate(
    variable = ifelse(stat == "std.error","",variable)
    #dif = ifelse(stat == "std.error","",dif)
  ) %>% 
  select(-stat) %>% 
   mutate_all(as.character) %>% 
  mutate_all(replace_na,"") %>% 
  mutate_all(str_remove_all,"\\\\c") %>% 
      mutate_all(str_remove_all,"\\{") %>% 
      mutate_all(str_remove_all,"\\}") %>% 
      mutate_all(str_remove_all,"\\^") %>% 
  kable(
    col.names = c(
      "Variable",
      str_c("(",1:3,")")
    )
  )
```

| Variable                   | \(1\)        | \(2\)        | \(3\)      |
|:---------------------------|:-------------|:-------------|:-----------|
| Economically active        | -0.067\*\*\* | -0.028       | 0.039      |
|                            | (0.018)      | (0.020)      | (0.027)    |
| Unemployed                 | 0.072\*\*\*  | 0.036\*\*\*  | -0.036\*\* |
|                            | (0.011)      | (0.011)      | (0.015)    |
| Male                       | 0.008        | 0.005        | -0.003     |
|                            | (0.012)      | (0.013)      | (0.018)    |
| Age                        | 1.198\*\*\*  | 0.712        | -0.486     |
|                            | (0.460)      | (0.537)      | (0.703)    |
| Primary education          | -0.025\*     | -0.069\*\*\* | -0.044\*   |
|                            | (0.015)      | (0.020)      | (0.025)    |
| Secondary education        | 0.034\*\*    | 0.048\*\*    | 0.015      |
|                            | (0.017)      | (0.022)      | (0.027)    |
| Tertiary education         | -0.008       | 0.021        | 0.029      |
|                            | (0.024)      | (0.017)      | (0.029)    |
| Born in Brno               | 0.061\*\*    | 0.060\*\*\*  | -0.0004    |
|                            | (0.027)      | (0.018)      | (0.032)    |
| Individuals in 2001 census | 864          | 744          |            |
| Individuals in 1991 census | 882          | 750          |            |

``` r
auxD <- TimeTravel %>% 
  group_by(adresado_01_klic, variable, census) %>% 
  summarise(
    mval = mean(value, na.rm = TRUE)
  ) %>% 
  spread(census,mval) %>% 
  mutate(
    dif = `2001` - `1991`
  ) %>% 
  select(adresado_01_klic,variable,dif) %>% 
  spread(variable,dif)

TimeTravel %>% 
  group_by(adresado_01_klic) %>% 
  summarise(
    W = mean(W),
    ctg = first(ctg)
  ) %>% 
  ungroup() %>% 
  left_join(auxD) %>% 
  felm(
    I(ctg == "treatment") ~ BornBrno + EDUCmiddle + EDUChigh + LVEK + MALE,
    data = .,
    weights = .$W
  ) %>% 
  summary(robust = TRUE)
```

    ## 
    ## Call:
    ##    felm(formula = I(ctg == "treatment") ~ BornBrno + EDUCmiddle +      EDUChigh + LVEK + MALE, data = ., weights = .$W) 
    ## 
    ## Weighted Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.60346 -0.50417 -0.00117  0.46859  0.58929 
    ## 
    ## Coefficients:
    ##              Estimate Robust s.e t value Pr(>|t|)    
    ## (Intercept)  0.515110   0.070496   7.307  1.3e-10 ***
    ## BornBrno    -0.027326   0.410226  -0.067    0.947    
    ## EDUCmiddle   0.080020   0.503253   0.159    0.874    
    ## EDUChigh     0.313396   0.482853   0.649    0.518    
    ## LVEK        -0.006053   0.015020  -0.403    0.688    
    ## MALE        -0.386525   0.552284  -0.700    0.486    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.4978 on 86 degrees of freedom
    ## Multiple R-squared(full model): 0.01176   Adjusted R-squared: -0.04569 
    ## Multiple R-squared(proj model): 0.01176   Adjusted R-squared: -0.04569 
    ## F-statistic(full model, *iid*):0.2048 on 5 and 86 DF, p-value: 0.9597 
    ## F-statistic(proj model): 0.2316 on 5 and 86 DF, p-value: 0.9477

F-statistic(full model, *iid*):0.2048 on 5 and 86 DF, p-value: 0.9597

``` r
balance_vars <- c(all.vars(model_propscore)[-1],"EDUClow","pscore")
balance_data <- r91_houses_sample1 %>% select(adresado_01_klic, one_of(balance_vars))
  
balance_test <- function(xmatch,xdata){
  aux <- left_join(xmatch$match, xdata) %>% 
    gather(variable,value,-id,-W,-ctg,-adresado_01_klic) %>% 
    mutate(ctg = factor(ctg) %>% relevel(ref = "control")) %>% 
    group_by(variable) %>% 
    do(
      lm = lm(value~ctg, data = ., weights = W) %>% coeftest(vcov. = vcovHC) %>% tidy() %>% filter(str_detect(term,"ctg"))
    ) %>% 
    unnest(lm)
  
  Pp <- aux %>% filter(variable == "pscore") %>% pull(p.value)
  Bp <- aux %>% filter(variable != "pscore") %>% pull(p.value) %>% min()
  
  tribble(
    ~Pp, ~Bp,
    Pp, Bp
  )
}

balance_df <- outputs_U %>% 
  select(C,D,match) %>% 
  rowwise() %>%
  mutate(
    test = balance_test(match, balance_data) %>% list()
  ) %>% 
  ungroup() %>% 
  unnest(test) %>% 
  select(-match)
```

#### Table A4: Estimates of the effect of homeownership on unemployment under alternative matching parameters

``` r
outputs_U %>% 
  ungroup() %>% 
  left_join(.,select(balance_df,C,D,Pp,Bp)) %>% 
  rowwise() %>% 
  mutate(
    Aestimate = format(estimate, digits = 1, nsmall = 3, trim = TRUE, scientific = FALSE) %>% str_c(.,add_stars(p.value, latex = TRUE, ds = TRUE)),
    Bstd.error = format(std.error, digits = 1, nsmall = 3, trim = TRUE, scientific = FALSE) %>% str_c("(",.,")"),
    DBp = format(Bp,digits = 1, nsmall = 3, trim = TRUE, scientific = FALSE),
    EPp = format(Pp,digits = 1, nsmall = 3, trim = TRUE, scientific = FALSE),
    Cn_t = str_c("Np = ",n_t)
  ) %>% 
  ungroup() %>% 
  mutate(
    Dtests = str_c("Bp = ",DBp,", Pp = ",EPp)
  ) %>% 
  select(C,D,Aestimate,Bstd.error,Cn_t,Dtests) %>% 
  gather(stat,value,-C,-D) %>% 
  mutate_if(is.numeric,as.character) %>% 
  filter(D %in% c(100,400,700,1000)) %>% 
  filter(C %in% c("0.1","0.2","0.3","0.4","0.5","0.6","0.7","0.8","0.9")) %>% 
  spread(D,value) %>% 
  arrange(C,stat) %>% 
  mutate(
    C = ifelse(str_detect(stat,"^A"),C,"")
  ) %>% 
  select(C,`100`,`400`,`700`,`1000`) %>% 
  kable()
```

| C   | 100                                                                                                                          | 400                                                                                                                    | 700                                                                                                                          | 1000                                                                                                                         |
|:----|:-----------------------------------------------------------------------------------------------------------------------------|:-----------------------------------------------------------------------------------------------------------------------|:-----------------------------------------------------------------------------------------------------------------------------|:-----------------------------------------------------------------------------------------------------------------------------|
| 0.1 | -0.050![^{\*\*}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5E%7B%2A%2A%7D "^{**}")       | -0.039![^{\*}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5E%7B%2A%7D "^{*}")       | -0.089![^{\*\*\*}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5E%7B%2A%2A%2A%7D "^{***}") | -0.093![^{\*\*\*}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5E%7B%2A%2A%2A%7D "^{***}") |
|     | (0.025)                                                                                                                      | (0.022)                                                                                                                | (0.023)                                                                                                                      | (0.030)                                                                                                                      |
|     | Np = 30                                                                                                                      | Np = 31                                                                                                                | Np = 39                                                                                                                      | Np = 41                                                                                                                      |
|     | Bp = 0.570, Pp = 0.967                                                                                                       | Bp = 0.283, Pp = 0.954                                                                                                 | Bp = 0.168, Pp = 0.896                                                                                                       | Bp = 0.201, Pp = 0.892                                                                                                       |
| 0.2 | -0.055![^{\*\*\*}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5E%7B%2A%2A%2A%7D "^{***}") | -0.044![^{\*\*}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5E%7B%2A%2A%7D "^{**}") | -0.064![^{\*\*\*}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5E%7B%2A%2A%2A%7D "^{***}") | -0.089![^{\*\*\*}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5E%7B%2A%2A%2A%7D "^{***}") |
|     | (0.021)                                                                                                                      | (0.018)                                                                                                                | (0.024)                                                                                                                      | (0.028)                                                                                                                      |
|     | Np = 36                                                                                                                      | Np = 36                                                                                                                | Np = 44                                                                                                                      | Np = 45                                                                                                                      |
|     | Bp = 0.437, Pp = 0.987                                                                                                       | Bp = 0.230, Pp = 0.986                                                                                                 | Bp = 0.234, Pp = 0.923                                                                                                       | Bp = 0.195, Pp = 0.926                                                                                                       |
| 0.3 | -0.053![^{\*\*}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5E%7B%2A%2A%7D "^{**}")       | -0.041![^{\*\*}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5E%7B%2A%2A%7D "^{**}") | -0.058![^{\*\*}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5E%7B%2A%2A%7D "^{**}")       | -0.089![^{\*\*\*}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5E%7B%2A%2A%2A%7D "^{***}") |
|     | (0.021)                                                                                                                      | (0.018)                                                                                                                | (0.023)                                                                                                                      | (0.028)                                                                                                                      |
|     | Np = 39                                                                                                                      | Np = 39                                                                                                                | Np = 45                                                                                                                      | Np = 45                                                                                                                      |
|     | Bp = 0.546, Pp = 0.985                                                                                                       | Bp = 0.311, Pp = 0.961                                                                                                 | Bp = 0.188, Pp = 0.939                                                                                                       | Bp = 0.195, Pp = 0.926                                                                                                       |
| 0.4 | -0.052![^{\*\*\*}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5E%7B%2A%2A%2A%7D "^{***}") | -0.042![^{\*\*}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5E%7B%2A%2A%7D "^{**}") | -0.058![^{\*\*}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5E%7B%2A%2A%7D "^{**}")       | -0.089![^{\*\*\*}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5E%7B%2A%2A%2A%7D "^{***}") |
|     | (0.020)                                                                                                                      | (0.017)                                                                                                                | (0.023)                                                                                                                      | (0.028)                                                                                                                      |
|     | Np = 40                                                                                                                      | Np = 40                                                                                                                | Np = 45                                                                                                                      | Np = 45                                                                                                                      |
|     | Bp = 0.507, Pp = 0.983                                                                                                       | Bp = 0.286, Pp = 0.994                                                                                                 | Bp = 0.188, Pp = 0.939                                                                                                       | Bp = 0.195, Pp = 0.926                                                                                                       |
| 0.5 | -0.051![^{\*\*}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5E%7B%2A%2A%7D "^{**}")       | -0.041![^{\*\*}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5E%7B%2A%2A%7D "^{**}") | -0.058![^{\*\*}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5E%7B%2A%2A%7D "^{**}")       | -0.089![^{\*\*\*}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5E%7B%2A%2A%2A%7D "^{***}") |
|     | (0.020)                                                                                                                      | (0.017)                                                                                                                | (0.023)                                                                                                                      | (0.028)                                                                                                                      |
|     | Np = 41                                                                                                                      | Np = 41                                                                                                                | Np = 45                                                                                                                      | Np = 45                                                                                                                      |
|     | Bp = 0.520, Pp = 0.982                                                                                                       | Bp = 0.295, Pp = 0.958                                                                                                 | Bp = 0.188, Pp = 0.939                                                                                                       | Bp = 0.195, Pp = 0.926                                                                                                       |
| 0.6 | -0.051![^{\*\*}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5E%7B%2A%2A%7D "^{**}")       | -0.041![^{\*\*}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5E%7B%2A%2A%7D "^{**}") | -0.058![^{\*\*}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5E%7B%2A%2A%7D "^{**}")       | -0.089![^{\*\*\*}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5E%7B%2A%2A%2A%7D "^{***}") |
|     | (0.020)                                                                                                                      | (0.017)                                                                                                                | (0.023)                                                                                                                      | (0.028)                                                                                                                      |
|     | Np = 41                                                                                                                      | Np = 41                                                                                                                | Np = 45                                                                                                                      | Np = 45                                                                                                                      |
|     | Bp = 0.520, Pp = 0.982                                                                                                       | Bp = 0.295, Pp = 0.958                                                                                                 | Bp = 0.188, Pp = 0.939                                                                                                       | Bp = 0.195, Pp = 0.926                                                                                                       |
| 0.7 | -0.051![^{\*\*\*}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5E%7B%2A%2A%2A%7D "^{***}") | -0.041![^{\*\*}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5E%7B%2A%2A%7D "^{**}") | -0.060![^{\*\*\*}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5E%7B%2A%2A%2A%7D "^{***}") | -0.090![^{\*\*\*}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5E%7B%2A%2A%2A%7D "^{***}") |
|     | (0.020)                                                                                                                      | (0.017)                                                                                                                | (0.022)                                                                                                                      | (0.027)                                                                                                                      |
|     | Np = 42                                                                                                                      | Np = 42                                                                                                                | Np = 46                                                                                                                      | Np = 46                                                                                                                      |
|     | Bp = 0.579, Pp = 0.966                                                                                                       | Bp = 0.338, Pp = 0.986                                                                                                 | Bp = 0.173, Pp = 0.892                                                                                                       | Bp = 0.163, Pp = 0.879                                                                                                       |
| 0.8 | -0.051![^{\*\*\*}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5E%7B%2A%2A%2A%7D "^{***}") | -0.041![^{\*\*}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5E%7B%2A%2A%7D "^{**}") | -0.060![^{\*\*\*}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5E%7B%2A%2A%2A%7D "^{***}") | -0.090![^{\*\*\*}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5E%7B%2A%2A%2A%7D "^{***}") |
|     | (0.020)                                                                                                                      | (0.017)                                                                                                                | (0.022)                                                                                                                      | (0.027)                                                                                                                      |
|     | Np = 42                                                                                                                      | Np = 42                                                                                                                | Np = 46                                                                                                                      | Np = 46                                                                                                                      |
|     | Bp = 0.579, Pp = 0.966                                                                                                       | Bp = 0.338, Pp = 0.986                                                                                                 | Bp = 0.173, Pp = 0.892                                                                                                       | Bp = 0.163, Pp = 0.879                                                                                                       |
| 0.9 | -0.051![^{\*\*\*}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5E%7B%2A%2A%2A%7D "^{***}") | -0.041![^{\*\*}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5E%7B%2A%2A%7D "^{**}") | -0.060![^{\*\*\*}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5E%7B%2A%2A%2A%7D "^{***}") | -0.090![^{\*\*\*}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5E%7B%2A%2A%2A%7D "^{***}") |
|     | (0.020)                                                                                                                      | (0.017)                                                                                                                | (0.022)                                                                                                                      | (0.027)                                                                                                                      |
|     | Np = 42                                                                                                                      | Np = 42                                                                                                                | Np = 46                                                                                                                      | Np = 46                                                                                                                      |
|     | Bp = 0.579, Pp = 0.966                                                                                                       | Bp = 0.338, Pp = 0.986                                                                                                 | Bp = 0.173, Pp = 0.892                                                                                                       | Bp = 0.163, Pp = 0.879                                                                                                       |

#### Table A3: Estimates of the effect of homeownership on labor force participation under alternative matching parameters

``` r
outputs_EA %>% 
  ungroup() %>% 
  left_join(.,select(balance_df,C,D,Pp,Bp)) %>% 
  #left_join(.,select(psttest,C,D,psttest)) %>% 
  rowwise() %>% 
  mutate(
    Aestimate = format(estimate, digits = 1, nsmall = 3, trim = TRUE, scientific = FALSE) %>% str_c(.,add_stars(p.value, latex = TRUE, ds = TRUE)),
    Bstd.error = format(std.error, digits = 1, nsmall = 3, trim = TRUE, scientific = FALSE) %>% str_c("(",.,")"),
    DBp = format(Bp,digits = 1, nsmall = 3, trim = TRUE, scientific = FALSE),
    EPp = format(Pp,digits = 1, nsmall = 3, trim = TRUE, scientific = FALSE),
    Cn_t = str_c("Np = ",n_t)
  ) %>% 
  ungroup() %>% 
  mutate(
    Dtests = str_c("Bp = ",DBp,", Pp = ",EPp)
  ) %>% 
  select(C,D,Aestimate,Bstd.error,Cn_t,Dtests) %>% 
  gather(stat,value,-C,-D) %>% 
  mutate_if(is.numeric,as.character) %>% 
  filter(D %in% c(100,400,700,1000)) %>% 
  filter(C %in% c("0.1","0.2","0.3","0.4","0.5","0.6","0.7","0.8","0.9")) %>% 
  spread(D,value) %>% 
  arrange(C,stat) %>% 
  mutate(
    C = ifelse(str_detect(stat,"^A"),C,"")
  ) %>% 
  select(C,`100`,`400`,`700`,`1000`) %>% 
  kable()
```

| C   | 100                    | 400                    | 700                    | 1000                   |
|:----|:-----------------------|:-----------------------|:-----------------------|:-----------------------|
| 0.1 | -0.024                 | 0.003                  | -0.013                 | -0.007                 |
|     | (0.028)                | (0.027)                | (0.028)                | (0.031)                |
|     | Np = 30                | Np = 31                | Np = 39                | Np = 41                |
|     | Bp = 0.570, Pp = 0.967 | Bp = 0.283, Pp = 0.954 | Bp = 0.168, Pp = 0.896 | Bp = 0.201, Pp = 0.892 |
| 0.2 | -0.028                 | -0.004                 | -0.004                 | -0.006                 |
|     | (0.024)                | (0.024)                | (0.027)                | (0.028)                |
|     | Np = 36                | Np = 36                | Np = 44                | Np = 45                |
|     | Bp = 0.437, Pp = 0.987 | Bp = 0.230, Pp = 0.986 | Bp = 0.234, Pp = 0.923 | Bp = 0.195, Pp = 0.926 |
| 0.3 | -0.024                 | -0.003                 | 0.003                  | -0.006                 |
|     | (0.025)                | (0.026)                | (0.026)                | (0.028)                |
|     | Np = 39                | Np = 39                | Np = 45                | Np = 45                |
|     | Bp = 0.546, Pp = 0.985 | Bp = 0.311, Pp = 0.961 | Bp = 0.188, Pp = 0.939 | Bp = 0.195, Pp = 0.926 |
| 0.4 | -0.026                 | -0.006                 | 0.003                  | -0.006                 |
|     | (0.025)                | (0.025)                | (0.026)                | (0.028)                |
|     | Np = 40                | Np = 40                | Np = 45                | Np = 45                |
|     | Bp = 0.507, Pp = 0.983 | Bp = 0.286, Pp = 0.994 | Bp = 0.188, Pp = 0.939 | Bp = 0.195, Pp = 0.926 |
| 0.5 | -0.022                 | -0.003                 | 0.003                  | -0.006                 |
|     | (0.025)                | (0.025)                | (0.026)                | (0.028)                |
|     | Np = 41                | Np = 41                | Np = 45                | Np = 45                |
|     | Bp = 0.520, Pp = 0.982 | Bp = 0.295, Pp = 0.958 | Bp = 0.188, Pp = 0.939 | Bp = 0.195, Pp = 0.926 |
| 0.6 | -0.022                 | -0.003                 | 0.003                  | -0.006                 |
|     | (0.025)                | (0.025)                | (0.026)                | (0.028)                |
|     | Np = 41                | Np = 41                | Np = 45                | Np = 45                |
|     | Bp = 0.520, Pp = 0.982 | Bp = 0.295, Pp = 0.958 | Bp = 0.188, Pp = 0.939 | Bp = 0.195, Pp = 0.926 |
| 0.7 | -0.020                 | -0.001                 | 0.003                  | -0.003                 |
|     | (0.024)                | (0.024)                | (0.025)                | (0.027)                |
|     | Np = 42                | Np = 42                | Np = 46                | Np = 46                |
|     | Bp = 0.579, Pp = 0.966 | Bp = 0.338, Pp = 0.986 | Bp = 0.173, Pp = 0.892 | Bp = 0.163, Pp = 0.879 |
| 0.8 | -0.020                 | -0.001                 | 0.003                  | -0.003                 |
|     | (0.024)                | (0.024)                | (0.025)                | (0.027)                |
|     | Np = 42                | Np = 42                | Np = 46                | Np = 46                |
|     | Bp = 0.579, Pp = 0.966 | Bp = 0.338, Pp = 0.986 | Bp = 0.173, Pp = 0.892 | Bp = 0.163, Pp = 0.879 |
| 0.9 | -0.020                 | -0.001                 | 0.003                  | -0.003                 |
|     | (0.024)                | (0.024)                | (0.025)                | (0.027)                |
|     | Np = 42                | Np = 42                | Np = 46                | Np = 46                |
|     | Bp = 0.579, Pp = 0.966 | Bp = 0.338, Pp = 0.986 | Bp = 0.173, Pp = 0.892 | Bp = 0.163, Pp = 0.879 |

#### Table 4: Homeownership and labor market activity, IV estimates for individuals living in matched houses

``` r
### Regression data
model_iv <- modeliv_individual #%>% lapply(matchingFE)

# Remove incomplete observations
model_iv %>% 
  lapply(all.vars) %>% 
  unlist() %>% 
  unique() %>% 
  c("flat_ID","sample","RZSJ","MALE","adresado_01_klic","person_ID",.) -> all_vars

r01_data_x <- r01_data_basic %>%
  dplyr::select(one_of(all_vars)) %>% 
  drop_na()

r01_data_x %>% 
  filter(sample %in% c("ATT","NTR")) %>% 
  mutate(
    ownership = ownership == "owner",
    sample_dummy = sample == "ATT"
    ) -> reg_data

reg_data_ea <- reg_data %>% filter(LEKAKTI %in% c(1,5))

reg_data_U <- outputs_U %>% 
  filter(isTRUE(all.equal(C, d_caliper)), isTRUE(all.equal(D, d_distance))) %>% 
  magrittr::extract2("match") %>% 
  magrittr::extract2(1) %>% 
  magrittr::extract2("match") %>% 
  left_join(reg_data_ea)

reg_data_EA <- outputs_U %>% 
  filter(isTRUE(all.equal(C, d_caliper)), isTRUE(all.equal(D, d_distance))) %>% 
  magrittr::extract2("match") %>% 
  magrittr::extract2(1) %>% 
  magrittr::extract2("match") %>% 
  left_join(reg_data)

## IV estimation

c(
  modeliv[1:5] %>% lapply(function(x) ivreg(x, data = reg_data_U, weights = W)),
  modeliv[6:10] %>% lapply(function(x) ivreg(x, data = reg_data_EA, weights = W))
  ) -> estm

if(cluster_by == "flat"){
  c(
    estm[1:5] %>% lapply(function(x) cluster.robust.se(x, reg_data_U$flat_ID)),
    estm[6:10] %>% lapply(function(x) cluster.robust.se(x, reg_data_EA$flat_ID))
  ) -> estm_ctest
}else{
  c(
    estm[1:5] %>% lapply(function(x) cluster.robust.se(x, reg_data_U$adresado_01_klic)),
    estm[6:10] %>% lapply(function(x) cluster.robust.se(x, reg_data_EA$adresado_01_klic))
  ) -> estm_ctest
}
```

    ## [1] "Cluster Robust Standard Errors"
    ## [1] "Cluster Robust Standard Errors"
    ## [1] "Cluster Robust Standard Errors"
    ## [1] "Cluster Robust Standard Errors"
    ## [1] "Cluster Robust Standard Errors"
    ## [1] "Cluster Robust Standard Errors"
    ## [1] "Cluster Robust Standard Errors"
    ## [1] "Cluster Robust Standard Errors"
    ## [1] "Cluster Robust Standard Errors"
    ## [1] "Cluster Robust Standard Errors"

``` r
get_tables(estm, estm_ctest, ivtest = TRUE, 
           reorder = c(6,7,8,9,10,1,2,3,4,5), space = FALSE) %>% 
  mutate_all(replace_na,"") %>% 
  kable(
    col.names = c(
      "Variable",
      str_c("(",1:10,")")
    )
  )
```

| Variable                                                                                    | \(1\)       | \(2\)   | \(3\)       | \(4\)       | \(5\)       | \(6\)       | \(7\)        | \(8\)        | \(9\)        | \(10\)       |
|:--------------------------------------------------------------------------------------------|:------------|:--------|:------------|:------------|:------------|:------------|:-------------|:-------------|:-------------|:-------------|
| Homeowner (=1)                                                                              | 0.022       | 0.020   | 0.013       | 0.006       | 0.003       | -0.057\*\*  | -0.061\*\*\* | -0.044\*     | -0.047\*\*   | -0.060\*\*\* |
|                                                                                             | (0.025)     | (0.027) | (0.026)     | (0.022)     | (0.025)     | (0.024)     | (0.021)      | (0.023)      | (0.022)      | (0.022)      |
| Secondary education (=1)                                                                    |             |         | -0.024      | 0.028       | 0.032       |             |              | -0.097\*\*\* | -0.106\*\*\* | -0.092\*\*\* |
|                                                                                             |             |         | (0.035)     | (0.031)     | (0.032)     |             |              | (0.028)      | (0.028)      | (0.029)      |
| Tertiary education (=1)                                                                     |             |         | 0.095\*\*\* | 0.077\*\*\* | 0.089\*\*\* |             |              | -0.162\*\*\* | -0.148\*\*\* | -0.135\*\*\* |
|                                                                                             |             |         | (0.034)     | (0.027)     | (0.030)     |             |              | (0.026)      | (0.026)      | (0.029)      |
| Male (=1)                                                                                   |             |         |             | 0.155\*\*\* | 0.156\*\*\* |             |              |              | 0.034        | 0.030        |
|                                                                                             |             |         |             | (0.021)     | (0.022)     |             |              |              | (0.021)      | (0.021)      |
| Born in Brno (=1)                                                                           |             |         |             | 0.018       | 0.016       |             |              |              | -0.011       | 0.000007     |
|                                                                                             |             |         |             | (0.024)     | (0.025)     |             |              |              | (0.023)      | (0.024)      |
| Incomplete family (=1)                                                                      |             |         |             | -0.016      | -0.007      |             |              |              | 0.070\*\*    | 0.068\*\*    |
|                                                                                             |             |         |             | (0.026)     | (0.027)     |             |              |              | (0.032)      | (0.030)      |
| One-person family (=1)                                                                      |             |         |             | 0.008       | 0.021       |             |              |              | 0.008        | 0.011        |
|                                                                                             |             |         |             | (0.032)     | (0.034)     |             |              |              | (0.028)      | (0.030)      |
| Constant                                                                                    | 0.768\*\*\* |         | 0.754\*\*\* |             |             | 0.130\*\*\* |              | 0.210\*\*\*  |              |              |
|                                                                                             | (0.017)     |         | (0.022)     |             |             | (0.017)     |              | (0.026)      |              |              |
| Age dummies                                                                                 | –           | –       | –           | Yes         | Yes         | –           | –            | –            | Yes          | Yes          |
| Neighborhood dummies                                                                        | –           | Yes     | –           | –           | Yes         | –           | Yes          | –            | –            | Yes          |
| Observations                                                                                | 1,238       | 1,238   | 1,238       | 1,238       | 1,238       | 960         | 960          | 960          | 960          | 960          |
| ![R2](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;R2 "R2") | 0.001       | 0.018   | 0.015       | 0.255       | 0.265       | 0.012       | 0.063        | 0.056        | 0.140        | 0.172        |

#### Table 5: Homeownership and labor market activity, intention-to-treat (reduced-form) estimates for individuals living in matched houses

``` r
### Reduced form
model_reduced <- model %>% 
  lapply(update, . ~ sample_dummy + . - ownership)

c(
  model_reduced[1:5] %>% lapply(function(x) lm(x, data = reg_data_U, weights = W)),
  model_reduced[6:10] %>% lapply(function(x) lm(x, data = reg_data_EA, weights = W))
  ) -> estm


if(cluster_by == "flat"){
  c(
    estm[1:5] %>% lapply(function(x) coeftest(x, vcov. = vcovCL(x, reg_data_U$adresado_01_klic))),
    estm[6:10] %>% lapply(function(x) coeftest(x, vcov. = vcovCL(x, reg_data_EA$adresado_01_klic)))
  ) -> estm_ctest
}else{
  c(
    estm[1:5] %>% lapply(function(x) coeftest(x, vcov. = vcovCL(x, reg_data_U$adresado_01_klic))),
    estm[6:10] %>% lapply(function(x) coeftest(x, vcov. = vcovCL(x, reg_data_EA$adresado_01_klic)))
  ) -> estm_ctest
}

get_tables(estm, estm_ctest, space = FALSE, 
           reorder = c(6,7,8,9,10,1,2,3,4,5)) %>% 
  mutate_all(replace_na,"") %>% 
  kable(
    col.names = c(
      "Variable",
      str_c("(",1:10,")")
    )
  )
```

| Variable                                                                                    | \(1\)       | \(2\)   | \(3\)       | \(4\)       | \(5\)       | \(6\)       | \(7\)        | \(8\)        | \(9\)        | \(10\)       |
|:--------------------------------------------------------------------------------------------|:------------|:--------|:------------|:------------|:------------|:------------|:-------------|:-------------|:-------------|:-------------|
| Living in privatized house (=1)                                                             | 0.020       | 0.018   | 0.012       | 0.005       | 0.003       | -0.052\*\*  | -0.056\*\*\* | -0.040\*     | -0.043\*\*   | -0.056\*\*\* |
|                                                                                             | (0.023)     | (0.025) | (0.023)     | (0.020)     | (0.023)     | (0.022)     | (0.020)      | (0.021)      | (0.020)      | (0.021)      |
| Secondary education (=1)                                                                    |             |         | -0.023      | 0.028       | 0.032       |             |              | -0.099\*\*\* | -0.109\*\*\* | -0.095\*\*\* |
|                                                                                             |             |         | (0.034)     | (0.031)     | (0.032)     |             |              | (0.028)      | (0.028)      | (0.029)      |
| Tertiary education (=1)                                                                     |             |         | 0.096\*\*\* | 0.078\*\*\* | 0.089\*\*\* |             |              | -0.164\*\*\* | -0.151\*\*\* | -0.138\*\*\* |
|                                                                                             |             |         | (0.033)     | (0.026)     | (0.030)     |             |              | (0.026)      | (0.027)      | (0.029)      |
| Male (=1)                                                                                   |             |         |             | 0.155\*\*\* | 0.156\*\*\* |             |              |              | 0.034        | 0.031        |
|                                                                                             |             |         |             | (0.021)     | (0.022)     |             |              |              | (0.021)      | (0.021)      |
| Born in Brno (=1)                                                                           |             |         |             | 0.018       | 0.016       |             |              |              | -0.011       | -0.0008      |
|                                                                                             |             |         |             | (0.024)     | (0.025)     |             |              |              | (0.023)      | (0.024)      |
| Incomplete family (=1)                                                                      |             |         |             | -0.016      | -0.007      |             |              |              | 0.069\*\*    | 0.068\*\*    |
|                                                                                             |             |         |             | (0.026)     | (0.027)     |             |              |              | (0.032)      | (0.030)      |
| One-person family (=1)                                                                      |             |         |             | 0.008       | 0.020       |             |              |              | 0.007        | 0.012        |
|                                                                                             |             |         |             | (0.032)     | (0.034)     |             |              |              | (0.028)      | (0.030)      |
| Constant                                                                                    | 0.768\*\*\* |         | 0.754\*\*\* |             |             | 0.130\*\*\* |              | 0.211\*\*\*  |              |              |
|                                                                                             | (0.017)     |         | (0.022)     |             |             | (0.017)     |              | (0.026)      |              |              |
| Age dummies                                                                                 | –           | –       | –           | Yes         | Yes         | –           | –            | –            | Yes          | Yes          |
| Neighborhood dummies                                                                        | –           | Yes     | –           | –           | Yes         | –           | Yes          | –            | –            | Yes          |
| Observations                                                                                | 1,238       | 1,238   | 1,238       | 1,238       | 1,238       | 960         | 960          | 960          | 960          | 960          |
| ![R2](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;R2 "R2") | 0.0006      | 0.017   | 0.015       | 0.255       | 0.265       | 0.007       | 0.061        | 0.054        | 0.138        | 0.171        |

#### Table A2: First-stage regressions results for IV estimates in Table 4, matched sample

``` r
### First stage
modeliv <- modeliv_1ststage %>% lapply(update, ownership ~ .)

modeliv[[1]] <- ownership ~ I(as.numeric(sample_dummy)) - 1

c(
  modeliv %>% lapply(function(x) lm(x, data = reg_data_U, weights = W)),
  modeliv %>% lapply(function(x) lm(x, data = reg_data_EA, weights = W))
  ) -> estm


if(cluster_by == "flat"){
  c(
    estm[1:5] %>% lapply(function(x) coeftest(x, vcov. = vcovCL(x, reg_data_U$adresado_01_klic))),
    estm[6:10] %>% lapply(function(x) coeftest(x, vcov. = vcovCL(x, reg_data_EA$adresado_01_klic)))
  ) -> estm_ctest
}else{
  c(
    estm[1:5] %>% lapply(function(x) coeftest(x, vcov. = vcovCL(x, reg_data_U$adresado_01_klic))),
    estm[6:10] %>% lapply(function(x) coeftest(x, vcov. = vcovCL(x, reg_data_EA$adresado_01_klic)))
  ) -> estm_ctest
}

get_tables(estm, estm_ctest, space = FALSE, reorder = c(6,7,8,9,10,1,2,3,4,5)) %>% 
  mutate_all(replace_na,"") %>% 
  kable(
    col.names = c(
      "Variable",
      str_c("(",1:10,")")
    )
  )
```

| Variable                                                                                    | \(1\)       | \(2\)       | \(3\)       | \(4\)       | \(5\)       | \(6\)       | \(7\)       | \(8\)       | \(9\)       | \(10\)      |
|:--------------------------------------------------------------------------------------------|:------------|:------------|:------------|:------------|:------------|:------------|:------------|:------------|:------------|:------------|
| Living in privatized house (=1)                                                             | 0.911\*\*\* | 0.919\*\*\* | 0.907\*\*\* | 0.908\*\*\* | 0.921\*\*\* | 0.920\*\*\* | 0.930\*\*\* | 0.916\*\*\* | 0.919\*\*\* | 0.935\*\*\* |
|                                                                                             | (0.028)     | (0.024)     | (0.029)     | (0.028)     | (0.022)     | (0.028)     | (0.023)     | (0.029)     | (0.028)     | (0.021)     |
| Secondary education (=1)                                                                    |             |             | 0.045\*     | 0.048\*\*   | 0.046\*\*   |             |             | 0.052\*\*   | 0.058\*\*   | 0.050\*\*   |
|                                                                                             |             |             | (0.023)     | (0.022)     | (0.023)     |             |             | (0.020)     | (0.022)     | (0.021)     |
| Tertiary education (=1)                                                                     |             |             | 0.066\*\*\* | 0.070\*\*\* | 0.057\*\*   |             |             | 0.060\*\*   | 0.067\*\*   | 0.052\*\*   |
|                                                                                             |             |             | (0.025)     | (0.026)     | (0.023)     |             |             | (0.028)     | (0.029)     | (0.025)     |
| Male (=1)                                                                                   |             |             |             | 0.002       | 0.003       |             |             |             | -0.0007     | -0.002      |
|                                                                                             |             |             |             | (0.008)     | (0.008)     |             |             |             | (0.008)     | (0.010)     |
| Born in Brno (=1)                                                                           |             |             |             | 0.004       | 0.011       |             |             |             | 0.007       | 0.014       |
|                                                                                             |             |             |             | (0.012)     | (0.012)     |             |             |             | (0.013)     | (0.011)     |
| Incomplete family (=1)                                                                      |             |             |             | 0.009       | 0.002       |             |             |             | 0.011       | 0.003       |
|                                                                                             |             |             |             | (0.015)     | (0.013)     |             |             |             | (0.017)     | (0.014)     |
| One-person family (=1)                                                                      |             |             |             | -0.019      | -0.030      |             |             |             | 0.0005      | -0.018      |
|                                                                                             |             |             |             | (0.019)     | (0.020)     |             |             |             | (0.018)     | (0.019)     |
| Constant                                                                                    |             |             | -0.034\*\*  |             |             |             |             | -0.035\*\*  |             |             |
|                                                                                             |             |             | (0.014)     |             |             |             |             | (0.015)     |             |             |
| Age dummies                                                                                 | –           | –           | –           | Yes         | Yes         | –           | –           | –           | Yes         | Yes         |
| Neighborhood dummies                                                                        | –           | Yes         | –           | –           | Yes         | –           | Yes         | –           | –           | Yes         |
| Observations                                                                                | 1,238       | 1,238       | 1,238       | 1,238       | 1,238       | 960         | 960         | 960         | 960         | 960         |
| ![R2](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;R2 "R2") | 0.911       | 0.864       | 0.839       | 0.847       | 0.873       | 0.920       | 0.876       | 0.852       | 0.860       | 0.884       |

## Additional analysis

``` r
model_robust <- model_individual %>% lapply(update, . ~ I(sample != "NTR") + . - ownership)

# Remove incomplete observations
model_robust %>% 
  lapply(all.vars) %>% 
  unlist() %>% 
  unique() %>% 
  c("flat_ID","sample","RZSJ","MALE",.) -> all_vars

r01_data_x <- r01_data_basic %>%
  dplyr::select(one_of(all_vars),adresado_01_klic) %>% 
  drop_na()
```

Columns 1 and 5

``` r
privX <- r01_data %>% filter(!is.na(Act), priv_year == "2001-2005") %>% distinct(adresado_01_klic) %>% pull(adresado_01_klic)

# Houses sample
r91_houses %>% 
  filter(adresado_01_klic %in% privX | sample == "NTR") %>% 
  mutate(
    sample = ifelse(sample == "NTR","NTR","ATT")
  ) %>% 
  dplyr::select(-Act) %>% 
  filter(adresado_01_klic %in% r01_data_basic$adresado_01_klic) %>% 
  rename(JPOCBYT = POC_BYTU) %>% 
  mutate(
    byty = cut(JPOCBYT, c(0,10,20,40), right = FALSE)
  ) %>% 
  drop_na() %>% 
  rename(BornBrno = BB) -> r91_houses_sample1

att_selection <- glm(update(model_propscore, I(sample == "ATT") ~ . ), data = r91_houses_sample1, family = binomial("logit"))

r91_houses_sample1 %>% 
  predict.glm(att_selection, type = "response", newdata = .) -> r91_houses_sample1$pscore

reg_data_EA <- get_matched.geo(r91_houses_sample1, d_caliper, d_distance) %>% 
  magrittr::extract2("match") %>% 
  left_join(r01_data_x) %>% drop_na()

reg_data_U <- reg_data_EA %>% filter(LEKAKTI %in% c(1,5))

c(
  model_robust[1:5] %>% lapply(function(x) lm(x, data = reg_data_U, weights = W)),
  model_robust[6:10] %>% lapply(function(x) lm(x, data = reg_data_EA, weights = W))
  ) -> estm

if(cluster_by == "flat"){
  c(
    estm[1:5] %>% lapply(function(x) coeftest(x, vcov. = vcovCL(x, reg_data_U$flat_ID))),
    estm[6:10] %>% lapply(function(x) coeftest(x, vcov. = vcovCL(x, reg_data_EA$flat_ID)))
  ) -> estm_ctest
}else{
  c(
    estm[1:5] %>% lapply(function(x) coeftest(x, vcov. = vcovCL(x, reg_data_U$adresado_01_klic))),
    estm[6:10] %>% lapply(function(x) coeftest(x, vcov. = vcovCL(x, reg_data_EA$adresado_01_klic)))
  ) -> estm_ctest
}

tab2_model <- estm[c(10,5)]
tab2_ctest <- estm_ctest[c(10,5)]
```

Columns 3 and 7

``` r
# Houses sample
r91_houses %>% 
  filter((!is.na(Act) & sample == "NTC") | sample == "NTR") %>% 
  mutate(
    sample = ifelse(sample == "NTR","NTR","ATT")
  ) %>% 
  dplyr::select(-Act) %>% 
  filter(adresado_01_klic %in% r01_data_basic$adresado_01_klic) %>% 
  rename(JPOCBYT = POC_BYTU) %>% 
  mutate(
    byty = cut(JPOCBYT, c(0,10,20,40), right = FALSE)
  ) %>% 
  drop_na() %>% 
  rename(BornBrno = BB) -> r91_houses_sample1

att_selection <- glm(update(model_propscore, I(sample == "ATT") ~ .), data = r91_houses_sample1, family = binomial("logit"))

r91_houses_sample1 %>% 
  predict.glm(att_selection, type = "response", newdata = .) -> r91_houses_sample1$pscore

reg_data_EA <- get_matched.geo(r91_houses_sample1, d_caliper, d_distance) %>% 
  magrittr::extract2("match") %>% 
  left_join(r01_data_x)

reg_data_U <- reg_data_EA %>% filter(LEKAKTI %in% c(1,5))

c(
  model_robust[1:5] %>% lapply(function(x) lm(x, data = reg_data_U, weights = W)),
  model_robust[6:10] %>% lapply(function(x) lm(x, data = reg_data_EA, weights = W))
  ) -> estm

if(cluster_by == "flat"){
  c(
    estm[1:5] %>% lapply(function(x) coeftest(x, vcov. = vcovCL(x, reg_data_U$flat_ID))),
    estm[6:10] %>% lapply(function(x) coeftest(x, vcov. = vcovCL(x, reg_data_EA$flat_ID)))
  ) -> estm_ctest
}else{
  c(
    estm[1:5] %>% lapply(function(x) coeftest(x, vcov. = vcovCL(x, reg_data_U$adresado_01_klic))),
    estm[6:10] %>% lapply(function(x) coeftest(x, vcov. = vcovCL(x, reg_data_EA$adresado_01_klic)))
  ) -> estm_ctest
}

tab3_model <- estm[c(10,5)]
tab3_ctest <- estm_ctest[c(10,5)]
```

Columns 2 and 6

``` r
privX <- r01_data %>% filter(!is.na(Act), priv_year == ">2005") %>% distinct(adresado_01_klic) %>% pull(adresado_01_klic)

# Houses sample
r91_houses %>% 
  filter(adresado_01_klic %in% privX | sample == "NTR") %>% 
  mutate(
    sample = ifelse(sample == "NTR","NTR","ATT")
  ) %>% 
  dplyr::select(-Act) %>% 
  filter(adresado_01_klic %in% r01_data_basic$adresado_01_klic) %>% 
  rename(JPOCBYT = POC_BYTU) %>% 
  mutate(
    byty = cut(JPOCBYT, c(0,10,20,40), right = FALSE)
  ) %>% 
  drop_na() %>% 
  rename(BornBrno = BB) -> r91_houses_sample1

att_selection <- glm(update(model_propscore, I(sample == "ATT") ~ .), data = r91_houses_sample1, family = binomial("logit"))

r91_houses_sample1 %>% 
  predict.glm(att_selection, type = "response", newdata = .) -> r91_houses_sample1$pscore

reg_data_EA <- get_matched.geo(r91_houses_sample1, d_caliper, d_distance) %>% 
  magrittr::extract2("match") %>% 
  left_join(r01_data_x)

reg_data_U <- reg_data_EA %>% filter(LEKAKTI %in% c(1,5))

c(
  model_robust[1:5] %>% lapply(function(x) lm(x, data = reg_data_U, weights = W)),
  model_robust[6:10] %>% lapply(function(x) lm(x, data = reg_data_EA, weights = W))
  ) -> estm

if(cluster_by == "flat"){
  c(
    estm[1:5] %>% lapply(function(x) coeftest(x, vcov. = vcovCL(x, reg_data_U$flat_ID))),
    estm[6:10] %>% lapply(function(x) coeftest(x, vcov. = vcovCL(x, reg_data_EA$flat_ID)))
  ) -> estm_ctest
}else{
  c(
    estm[1:5] %>% lapply(function(x) coeftest(x, vcov. = vcovCL(x, reg_data_U$adresado_01_klic))),
    estm[6:10] %>% lapply(function(x) coeftest(x, vcov. = vcovCL(x, reg_data_EA$adresado_01_klic)))
  ) -> estm_ctest
}

tab4_model <- estm[c(10,5)]
tab4_ctest <- estm_ctest[c(10,5)]
```

Columns 4 and 8

``` r
# Houses sample
r91_houses %>% 
  filter(is.na(Act) | sample == "NTR") %>% 
  mutate(
    sample = ifelse(sample == "NTR","NTR","ATT")
  ) %>% 
  dplyr::select(-Act) %>% 
  filter(adresado_01_klic %in% r01_data_basic$adresado_01_klic) %>% 
  rename(JPOCBYT = POC_BYTU) %>% 
  mutate(
    byty = cut(JPOCBYT, c(0,10,20,40), right = FALSE)
  ) %>% 
  drop_na() %>% 
  rename(BornBrno = BB) -> r91_houses_sample1

att_selection <- glm(update(model_propscore, I(sample == "ATT") ~ .), data = r91_houses_sample1, family = binomial("logit"))

r91_houses_sample1 %>% 
  predict.glm(att_selection, type = "response", newdata = .) -> r91_houses_sample1$pscore

reg_data_EA <- get_matched.geo(r91_houses_sample1, d_caliper, d_distance) %>% 
  magrittr::extract2("match") %>% 
  left_join(r01_data_x)

reg_data_U <- reg_data_EA %>% filter(LEKAKTI %in% c(1,5))

c(
  model_robust[1:5] %>% lapply(function(x) lm(x, data = reg_data_U, weights = W)),
  model_robust[6:10] %>% lapply(function(x) lm(x, data = reg_data_EA, weights = W))
  ) -> estm

if(cluster_by == "flat"){
  c(
    estm[1:5] %>% lapply(function(x) coeftest(x, vcov. = vcovCL(x, reg_data_U$flat_ID))),
    estm[6:10] %>% lapply(function(x) coeftest(x, vcov. = vcovCL(x, reg_data_EA$flat_ID)))
  ) -> estm_ctest
}else{
  c(
    estm[1:5] %>% lapply(function(x) coeftest(x, vcov. = vcovCL(x, reg_data_U$adresado_01_klic))),
    estm[6:10] %>% lapply(function(x) coeftest(x, vcov. = vcovCL(x, reg_data_EA$adresado_01_klic)))
  ) -> estm_ctest
}

tab5_model <- estm[c(10,5)]
tab5_ctest <- estm_ctest[c(10,5)]
```

#### Table 7: Labor market activity of individuals living in the public housing stock as of 2001, by designation and privatization status, and individuals in restituted houses, OLS estimates on matched samples

``` r
estm <- c(#tab1_model, 
          tab2_model, 
          tab4_model,
          tab3_model,
          tab5_model)
estm_ctest <- c(#tab1_ctest, 
                tab2_ctest,
                tab4_ctest,
                tab3_ctest,
                tab5_ctest)

get_tables(estm, estm_ctest, space = FALSE) %>% 
  mutate_all(replace_na,"") %>% 
  kable(
    col.names = c(
      "Variable",
      str_c("(",1:8,")")
    )
  )
```

| Variable                                                                                    | \(1\)       | \(2\)        | \(3\)       | \(4\)        | \(5\)       | \(6\)        | \(7\)       | \(8\)        |
|:--------------------------------------------------------------------------------------------|:------------|:-------------|:------------|:-------------|:------------|:-------------|:------------|:-------------|
| Living in city-owned house (=1)                                                             | -0.026      | -0.0009      | 0.038       | 0.020        | 0.005       | -0.007       | 0.008       | 0.0002       |
|                                                                                             | (0.021)     | (0.025)      | (0.029)     | (0.032)      | (0.035)     | (0.029)      | (0.008)     | (0.008)      |
| Secondary education (=1)                                                                    | -0.050\*    | -0.108\*\*\* | 0.054\*     | -0.055\*     | 0.032       | -0.062\*     | 0.028\*\*\* | -0.124\*\*\* |
|                                                                                             | (0.027)     | (0.034)      | (0.031)     | (0.032)      | (0.033)     | (0.033)      | (0.008)     | (0.010)      |
| Tertiary education (=1)                                                                     | 0.046\*     | -0.161\*\*\* | 0.103\*\*\* | -0.089\*\*\* | 0.074\*\*   | -0.086\*\*\* | 0.087\*\*\* | -0.159\*\*\* |
|                                                                                             | (0.024)     | (0.028)      | (0.033)     | (0.030)      | (0.035)     | (0.031)      | (0.008)     | (0.009)      |
| Male (=1)                                                                                   | 0.127\*\*\* | -0.030       | 0.177\*\*\* | 0.003        | 0.157\*\*\* | -0.008       | 0.161\*\*\* | -0.006       |
|                                                                                             | (0.020)     | (0.020)      | (0.032)     | (0.023)      | (0.025)     | (0.032)      | (0.007)     | (0.007)      |
| Born in Brno (=1)                                                                           | -0.026      | 0.032        | 0.028       | 0.016        | 0.062\*\*   | 0.016        | -0.003      | -0.001       |
|                                                                                             | (0.020)     | (0.023)      | (0.033)     | (0.031)      | (0.027)     | (0.029)      | (0.007)     | (0.007)      |
| Incomplete family (=1)                                                                      | 0.005       | 0.096\*\*\*  | 0.027       | 0.037        | 0.023       | 0.029        | 0.013       | 0.064\*\*\*  |
|                                                                                             | (0.033)     | (0.034)      | (0.050)     | (0.040)      | (0.034)     | (0.042)      | (0.009)     | (0.011)      |
| One-person family (=1)                                                                      | -0.002      | 0.066\*      | 0.081\*\*   | 0.050        | 0.035       | -0.038       | 0.011       | 0.017\*      |
|                                                                                             | (0.027)     | (0.035)      | (0.038)     | (0.054)      | (0.043)     | (0.038)      | (0.010)     | (0.010)      |
| Age dummies                                                                                 | Yes         | Yes          | Yes         | Yes          | Yes         | Yes          | Yes         | Yes          |
| Neighborhood dummies                                                                        | Yes         | Yes          | Yes         | Yes          | Yes         | Yes          | Yes         | Yes          |
| Observations                                                                                | 1,171       | 930          | 593         | 490          | 797         | 625          | 26,908      | 21,113       |
| ![R2](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;R2 "R2") | 0.301       | 0.128        | 0.299       | 0.210        | 0.273       | 0.181        | 0.212       | 0.105        |

#### Table 6: Year of privatization and labor market outcomes of individuals living in privatized houses (OLS, 2000 is the omitted category)

``` r
model_y <- model[c(5,10)] %>% 
  lapply(update, . ~ fyear + . - ownership)

all_vars <- model_y %>%  lapply(all.vars) %>% unlist() %>% unique() %>% c(.,"adresado_01_klic")

load("DATA_RP/IDOB.RData")
load("DATA_RP/privyear.RData")

selected_privatized <- left_join(IDOB,privyear) %>% 
  drop_na()

r01_data_x <- r01_data_basic %>% 
  filter(adresado_01_klic %in% selected_privatized$adresado_01_klic) %>% 
  left_join(selected_privatized) %>% 
  filter(year <= 2000) %>% 
  filter(!is.na(Act)) %>% 
  mutate(
    fyear = as.factor(year) %>% relevel(ref = "2000")
  ) %>%
  dplyr::select(one_of(all_vars)) %>% 
  drop_na()

reg_data_wa <- r01_data_x
reg_data_ea <- r01_data_x %>% filter(LEKAKTI %in% c(1,5))


c(
  model_y[1] %>% lapply(function(x) lm(x, data = reg_data_ea)),
  model_y[2] %>% lapply(function(x) lm(x, data = reg_data_wa))
  ) -> estm

if(cluster_by == "flat"){
  c(
    estm[1] %>% lapply(function(x) coeftest(x, vcov. = vcovCL(x, reg_data_ea$flat_ID))),
    estm[2] %>% lapply(function(x) coeftest(x, vcov. = vcovCL(x, reg_data_wa$flat_ID)))
  ) -> estm_ctest
}else{
  c(
    estm[1] %>% lapply(function(x) coeftest(x, vcov. = vcovCL(x, reg_data_ea$adresado_01_klic))),
    estm[2] %>% lapply(function(x) coeftest(x, vcov. = vcovCL(x, reg_data_wa$adresado_01_klic)))
  ) -> estm_ctest
}

    flev <- c(
      "fyear1998",
      "fyear1999",
      "fyear2001",
      "fyear2002",
      "fyear2003",
      "EDUCmiddle",
      "EDUChigh",
      "EXP",
      "I((EXP^2)/100)",
      "LPOHLAVmale",
      "BornBrno",
      "htypeIncomplete_family",
      "htypeOne_person_household",
      "(Intercept)",
      "LVEKfe",
      "RZSJfe"
    )
    
    flab <- c(
      "Living in a house privatized in 1998 (=1)",
      "Living in a house privatized in 1999 (=1)",
      "Living in a house privatized in 2001 (=1)",
      "Living in a house privatized in 2002 (=1)",
      "Living in a house privatized in 2003 (=1)",
      "Secondary education (=1)",
      "Tertiary education (=1)",
      "Potential experience",
      "Potential experience$^2/100$",
      "Male (=1)",
      "Born in Brno (=1)",
      "Incomplete family (=1)",
      "One-person family (=1)",
      "Constant",
      "Age dummies",
      "Neighborhood dummies"
    )
    
    ylabs <- list()
    ylabs$flev <- flev
    ylabs$flab <- flab

get_tables(estm, estm_ctest, labels = ylabs, space = FALSE, reorder = c(2,1)) %>% 
  mutate_all(replace_na,"") %>% 
  kable(
    col.names = c(
      "Variable",
      str_c("(",1:2,")")
    )
  )
```

| Variable                                                                                    | \(1\)       | \(2\)        |
|:--------------------------------------------------------------------------------------------|:------------|:-------------|
| Living in a house privatized in 1998 (=1)                                                   | 0.004       | -0.012       |
|                                                                                             | (0.043)     | (0.036)      |
| Living in a house privatized in 1999 (=1)                                                   | -0.026      | -0.033       |
|                                                                                             | (0.039)     | (0.039)      |
| Secondary education (=1)                                                                    | 0.030       | -0.065       |
|                                                                                             | (0.051)     | (0.045)      |
| Tertiary education (=1)                                                                     | 0.081\*     | -0.117\*\*\* |
|                                                                                             | (0.043)     | (0.043)      |
| Male (=1)                                                                                   | 0.146\*\*\* | 0.028        |
|                                                                                             | (0.032)     | (0.026)      |
| Born in Brno (=1)                                                                           | 0.019       | 0.011        |
|                                                                                             | (0.031)     | (0.028)      |
| Incomplete family (=1)                                                                      | 0.034       | 0.056        |
|                                                                                             | (0.040)     | (0.047)      |
| One-person family (=1)                                                                      | 0.051       | 0.004        |
|                                                                                             | (0.043)     | (0.044)      |
| Age dummies                                                                                 | Yes         | Yes          |
| Neighborhood dummies                                                                        | Yes         | Yes          |
| Observations                                                                                | 597         | 470          |
| ![R2](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;R2 "R2") | 0.267       | 0.190        |
