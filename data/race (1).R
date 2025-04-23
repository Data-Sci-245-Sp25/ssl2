

# calculating race data for later merging with ACS data

rawData<-read_csv("../multCoData/vr_mortality_multco_limited_2006_2023.csv")

createFips<- function(tract, suffix){
  
  tractFixed<-sprintf("%04g", tract)
  suffixFixed<-sprintf("%02g", suffix)
  
  fips<-paste("41051", tractFixed, suffixFixed, sep="")
  return (fips)
}
# generate fips code
dataWithFips<- rawData |>
  mutate(fips=createFips(ddrestract, ddrestractsuf))

# limit to cancer-only deaths
cancerOnlyData <- dataWithFips |>
  filter(startsWith(dmcaACME, "C") | startsWith(dmcaACME, "D"))



deaths <- cancerOnlyData %>%
  mutate(race_group = case_when(
    drace2e != "NaN" ~ "Two or More Races",
    drace1e >= 100 & drace1e <= 199 ~ "White alone",
    drace1e >= 200 & drace1e <= 299 ~ "Black or African American alone",
    (drace1e >= 300 & drace1e <= 399) | (drace1e >= "A01" & drace1e <= "R99") ~ "American Indian or Alaska Native alone",
    drace1e >= 400 & drace1e <= 499 ~ "Asian alone",
    drace1e >= 500 & drace1e <= 599 ~ "Native Hawaiian and Other Pacific Islander alone",
    drace1e >= 600 & drace1e <= 999 ~ "Some Other Race alone",
    TRUE ~ "Unknown"
  ))

# need to adjust for Hispanic
deaths<-deaths |>
  mutate(race_group = case_when(
    ddethnicmex=="H" | ddethnicpr=="H" | ddethniccuban=="H" | ddethnicoth=="H" ~ "Hispanic",
    TRUE ~ race_group
  )
  )


# classify cancer type:

icdCodes<-read_csv("../multCoData/icd10-codes.csv")

deaths <-deaths |>
  mutate(diseaseSite=icdCodes$Category[match(dmcaACME, icdCodes$`ICD-10 Code`)])

#limit to our years

deaths <- deaths |>
  filter(ddodyear>=2014 & ddodyear<=2023)

# get top 3 overall
top_cancers <- deaths %>%
  filter(!is.na(diseaseSite)) %>%
  count(diseaseSite, sort = TRUE) %>%
  slice_max(n, n = 3) %>%
  pull(diseaseSite)


# limit to just these
deaths_top3 <- deaths %>%
  filter(diseaseSite %in% top_cancers)


# generate aggregated data
summary <- deaths_top3 %>%
  group_by(ddodyear, race_group, diseaseSite) %>%
  summarize(deaths = n(), .groups = "drop")


write_csv(summary, "top3SummaryByRac.csv")



# is there a different top 3 for different races?


# Black or African American alone

blackTop3<- deaths |>
  filter(race_group=="Black or African American alone" & !is.na(diseaseSite)) |>
  count(diseaseSite, sort = TRUE) |>
  slice_max(n, n = 3) |>
  pull(diseaseSite)

blackSummary <- deaths |>
  filter(diseaseSite %in% blackTop3 & race_group=="Black or African American alone") |>
  group_by(ddodyear, diseaseSite) |>
  summarize(deaths = n(), .groups = "drop")

write_csv(blackSummary, "blackSummary.csv")


# Asian alone


asianTop3<- deaths |>
  filter(race_group=="Asian alone" & !is.na(diseaseSite)) |>
  count(diseaseSite, sort = TRUE) |>
  slice_max(n, n = 3) |>
  pull(diseaseSite)

asianSummary <- deaths |>
  filter(diseaseSite %in% asianTop3 & race_group=="Asian alone") |>
  group_by(ddodyear, diseaseSite) |>
  summarize(deaths = n(), .groups = "drop")
write_csv(asianSummary, "asianSummary.csv")

#Hispanic

hispanicTop3<- deaths |>
  filter(race_group=="Hispanic" & !is.na(diseaseSite)) |>
  count(diseaseSite, sort = TRUE) |>
  slice_max(n, n = 3) |>
  pull(diseaseSite)

hispanicSummary<-deaths |>
  filter(diseaseSite %in% hispanicTop3 & race_group=="Hispanic") |>
  group_by(ddodyear, diseaseSite) |>
  summarize(deaths = n(), .groups = "drop")

write_csv(hispanicSummary, "hispanicSummary.csv")



#White

whiteTop3<- deaths |>
  filter(race_group=="White alone" & !is.na(diseaseSite)) |>
  count(diseaseSite, sort = TRUE) |>
  slice_max(n, n = 3) |>
  pull(diseaseSite)


whiteSummary<-deaths |>
  filter(diseaseSite %in% whiteTop3 & race_group=="White alone") |>
  group_by(ddodyear, diseaseSite) |>
  summarize(deaths = n(), .groups = "drop")

write_csv(whiteSummary, "whiteSummary.csv")


# American Indian or Alaska Native alone

amIndianTop3<- deaths |>
  filter(race_group=="American Indian or Alaska Native alone" & !is.na(diseaseSite)) |>
  count(diseaseSite, sort = TRUE) |>
  slice_max(n, n = 3) |>
  pull(diseaseSite)


amIndianSummary<-deaths |>
  filter(diseaseSite %in% amIndianTop3 & race_group=="American Indian or Alaska Native alone") |>
  group_by(ddodyear, diseaseSite) |>
  summarize(deaths = n(), .groups = "drop")
write_csv(amIndianSummary, "amIndianSummary.csv")


# Two or More Races

twoOrMoreTop3<- deaths |>
  filter(race_group=="Two or More Races" & !is.na(diseaseSite)) |>
  count(diseaseSite, sort = TRUE) |>
  slice_max(n, n = 3) |>
  pull(diseaseSite)


twoOrMoreSummary<-deaths |>
  filter(diseaseSite %in% twoOrMoreTop3 & race_group=="Two or More Races") |>
  group_by(ddodyear, diseaseSite) |>
  summarize(deaths = n(), .groups = "drop")
write_csv(twoOrMoreSummary, "twoOrMoreSummary.csv")

# Native Hawaiian and Other Pacific Islander alone

hawTop3<- deaths |>
  filter(race_group=="Native Hawaiian and Other Pacific Islander alone" & !is.na(diseaseSite)) |>
  count(diseaseSite, sort = TRUE) |>
  slice_max(n, n = 3) |>
  pull(diseaseSite)


hawSummary<-deaths |>
  filter(diseaseSite %in% hawTop3 & race_group=="Native Hawaiian and Other Pacific Islander alone") |>
  group_by(ddodyear, diseaseSite) |>
  summarize(deaths = n(), .groups = "drop")
write_csv(hawSummary, "hawSummary.csv")


##### query and merge population data from tidycensus
install.packages("tidycensus")
#install.packages("tidyverse")
install.packages("sf")


#load libraries
library(tidycensus)
#library(tidyverse)
library(sf)

#set API key - enter your own, or use the key above if requesting a key didn't work
apiKey="69991184bb877aec0620418b2c2a16001401acf4"
census_api_key(apiKey, install = TRUE, overwrite=TRUE)

readRenviron("~/.Renviron")


# Define the years you want to analyze
years <- 2014:2023 # Adjust as needed

# Create empty lists to store results
pop_data_list <- list()

# Loop through years to get population data
for (year in years) {
  pop_data_list[[as.character(year)]] <- get_acs(
    geography = "county",
    variables = c(
      "Total Population" = "B03002_001",
      "Not Hispanic or Latino" = "B03002_002",
      "Hispanic or Latino" = "B03002_012",
      "White alone" = "B02001_002",
      "Black or African American alone" = "B02001_003",
      "American Indian or Alaska Native alone" = "B02001_004",
      "Asian alone" = "B02001_005",
      "Native Hawaiian and Other Pacific Islander alone" = "B02001_006",
      "Some Other Race alone" = "B02001_007",
      "Two or More Races" = "B02001_008"
    ),
    state = "OR",
    county = "Multnomah",
    year = year,
    survey = "acs5"
  ) |>
    mutate(year = year) # Add year column for easier joining later
}

pop_data_all_years <- bind_rows(pop_data_list)

# Reshape to get both estimates and MOEs
pop_data_wide <- pop_data_all_years %>%
  select(year, variable, estimate, moe) %>%
  pivot_wider(
    names_from = variable, 
    values_from = c(estimate, moe)
  )

### attempting again, to get sex vars

# Set up variables
years <- 2014:2023

race_vars <- c(
  "White alone" = "B02001_002",
  "Black or African American alone" = "B02001_003",
  "American Indian or Alaska Native alone" = "B02001_004",
  "Asian alone" = "B02001_005",
  "Native Hawaiian and Other Pacific Islander alone" = "B02001_006",
  "Some Other Race alone" = "B02001_007",
  "Two or More Races" = "B02001_008"
)

sex_vars <- c(
  "Total Male" = "B01001_002",
  "Total Female" = "B01001_026"
)

hispanic_vars <- c(
  "Total Population" = "B03002_001",
  "Not Hispanic or Latino" = "B03002_002",
  "Hispanic or Latino" = "B03002_012"
)

# Function to get data for one year
get_multnomah_1yr <- function(year) {
  # Helper that handles missing years
  safely_get <- function(vars, var_names, type_label) {
    df <- tryCatch(
      get_acs(
        geography = "county",
        variables = vars,
        state = "OR",
        county = "Multnomah",
        year = year,
        survey = "acs5"
      ),
      error = function(e) tibble()
    )
    
    if (nrow(df) == 0) return(tibble())
    
    df %>%
      mutate(
        category = names(vars)[match(variable, vars)],
        type = type_label
      )
  }
  
  race_data     <- safely_get(race_vars, race_vars, "race")
  sex_data      <- safely_get(sex_vars, sex_vars, "sex")
  hispanic_data <- safely_get(hispanic_vars, hispanic_vars, "ethnicity")
  
  bind_rows(race_data, sex_data, hispanic_data) %>%
    mutate(year = year)
}

# Loop through all years
acs_5yr_all <- map_dfr(years, get_multnomah_1yr)

write_csv(acs_5yr_all, "acs5_10yearData.csv")


# sample joining

#cancerRates2022 <- popData2022 |>
#  left_join(deathsByFips2022, by = c("GEOID" = "fips")) |>
#  mutate(rate_per_100k = (deaths / total_populationE) * 100000)


popAsian<-acs_5yr_all |>
  filter(variable=="Asian alone")


asianRates<- asianSummary |>
  left_join(popAsian, by=c("ddodyear"="year")) |>
  mutate(rate_per_100k = (deaths / estimate) * 100000)

ggplot(asianRates, mapping=aes(x=ddodyear, y=rate_per_100k, color=diseaseSite)) +
  geom_line()

ggsave("test.jpg")






