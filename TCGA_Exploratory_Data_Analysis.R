# ============================================================
# Exploratory Data Analysis: TCGA patient covariates by hospital
# ============================================================
# Loads patient-level TCGA data, cleans up inconsistent hospital
# ("site") names, pools patients into the 6 largest hospitals,
# and visualises how demographics, cancer type, and stage vary
# across those hospitals.

library(ggplot2)

data <- read.csv("tcga_patient_covariates.csv", stringsAsFactors = FALSE)

# Trim whitespace from the raw site names before any further cleaning
data$site_clean <- trimws(data$site)

# --- Check how many distinct cancer types each hospital treats ---
cancers_per_hospital <- table(data$site, data$cancer_type)
n_cancer_types <- rowSums(cancers_per_hospital > 0)

# Hospitals treating more than one cancer type
multi_cancer_hospitals <- n_cancer_types[n_cancer_types > 1]
print(sort(multi_cancer_hospitals, decreasing = TRUE))

# --- Map inconsistent/duplicate hospital name spellings onto one canonical name ---
site_map <- c(
  # MD Anderson
  "MD Anderson" = "MD Anderson Cancer Center",
  "University of Texas MD Anderson Cancer Center" = "MD Anderson Cancer Center",

  # Mayo Clinic
  "Mayo" = "Mayo Clinic",
  "Mayo Clinic Arizona" = "Mayo Clinic",
  "Mayo Clinic - Rochester" = "Mayo Clinic",
  "Mayo Clinic Rochester" = "Mayo Clinic",

  # ABS / BLN naming
  "ABS IUPUI" = "ABS - IUPUI",
  "BLN Baylor" = "BLN - Baylor",
  "BLN UT Southwestern Medical Center at Dallas" = "BLN - UT Southwestern Medical Center at Dallas",

  # Misc typos / capitalization
  "University Of Michigan" = "University of Michigan",
  "Michigan University" = "University of Michigan",
  "International Genomics Conosrtium" = "International Genomics Consortium",
  "International Genomics Consortium" = "Int'l. Genomics Consortium",
  "Greenville Health Systems" = "Greenville Health System",

  # St. Joseph's
  "St. Joseph's Hospital Arizona" = "St. Joseph's Hospital (AZ)",
  "St. Joseph's Hospital AZ" = "St. Joseph's Hospital (AZ)",
  "St. Joseph's Medical Center-(MD)" = "St. Joseph's Medical Center (MD)",
  "St Joseph's Medical Center (MD)" = "St. Joseph's Medical Center (MD)",

  # Washington University (core site only — NOT the "- X" consortium cluster)
  "Washington University St. Louis" = "Washington University",
  "Washington University - St. Louis" = "Washington University",

  # Memorial Sloan Kettering
  "Memorial Sloan Kettering Cancer Centre" = "Memorial Sloan Kettering",
  "Memorial Sloan Kettering Cancer Center" = "Memorial Sloan Kettering",
  "MSKCC" = "Memorial Sloan Kettering",

  # Other formatting duplicates found in data$site_clean
  "University Health Network, Toronto" = "University Health Network",
  "Case Western - St Joes" = "Case Western - St Joes",  # kept distinct — likely a sub-site, verify
  "Hartford" = "Hartford Hospital",
  "UCSF" = "University of California San Francisco",
  "Duke University" = "Duke",
  "Vanderbilt University" = "Vanderbilt",
  "Yale University" = "Yale",
  "Roswell" = "Roswell Park",
  "Cleveland Clinic Foundation" = "Cleveland Clinic",
  "Christiana Care" = "Christiana Healthcare",
  "Gundersen Lutheran Health System" = "Gundersen Lutheran",
  "The Ohio State University" = "Ohio State University",
  "UNC" = "University of North Carolina",
  "Fox Chase Cancer Center" = "Fox Chase",
  "ILSBIO" = "ILSBio",
  "ILSbio" = "ILSBio",
  "Proteogenex, Inc" = "Proteogenex, Inc.",
  "Ontario Institute for Cancer Research (OICR)" = "Ontario Institute for Cancer Research",
  "Ontario Institute for Cancer Research (OICR)/Ottawa" = "Ontario Institute for Cancer Research",
  "Thoraxklinik at University Hospital Heidelberg" = "Thoraxklinik",
  "BC Cancer Agency" = "British Columbia Cancer Agency",
  "Emory University - Winship Cancer Inst." = "Emory University",
  "St. University of Colorado Denver" = "University of Colorado Denver",
  "University of Kansas Medical Center" = "University of Kansas"
)

# Apply the name mapping, leaving any name not in site_map unchanged
data$site_clean <- ifelse(data$site_clean %in% names(site_map),
                          site_map[data$site_clean], data$site_clean)
data$site_clean[data$site_clean == ""] <- NA

# Count patients per (cleaned) hospital and pick the 6 largest
site_totals <- as.data.frame(table(data$site_clean))
colnames(site_totals) <- c("site_clean", "n_patients")
site_totals <- site_totals[order(-site_totals$n_patients), ]
print(site_totals)

top6_hospitals <- head(site_totals, 6)
top6_names <- top6_hospitals$site_clean

# Restrict to patients from the 6 largest hospitals, dropping rows with no site
subsample_top6 <- data[data$site_clean %in% top6_names, ]
grouped_data <- subsample_top6[!is.na(subsample_top6$site), ]

print(sort(table(grouped_data$site_clean), decreasing = TRUE))

# Recode missing/unknown race and sex values to explicit categories 
grouped_data$race[grouped_data$race %in% c("not reported", "Unknown", "")] <- "Unreported"

grouped_data$race[grouped_data$race %in% c("american indian or alaska native", "native hawaiian or other pacific islander")] <- "Other"
grouped_data$sex[grouped_data$sex %in% c("", "unknown")] <- "Unreported"

# Keep only the columns needed for the analysis below (stage is needed later)
grouped_data <- grouped_data[, c("site_clean", "cancer_type", "age_at_diagnosis", "sex", "race", "stage", "year_of_diagnosis")]

par(mfrow = c(1, 1))

# Age at Diagnosis by hospital
ggplot(grouped_data, aes(x = site_clean, y = age_at_diagnosis)) +
  geom_boxplot(fill = "lightblue") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        plot.title = element_text(hjust = 0.5)) +
  labs(x = "", y = "Age at Diagnosis", title = "Age at Diagnosis by Hospital")

# Year of Diagnosis by hospital
ggplot(grouped_data, aes(x = site_clean, y = year_of_diagnosis)) +
  geom_boxplot(fill = "lightblue") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        plot.title = element_text(hjust = 0.5)) +
  labs(x = "", y = "Year of Diagnosis", title = "Year of Diagnosis by Hospital")

# Collapse detailed AJCC substages into I/II/III/IV, or "Undisclosed" 
stage_map <- c(
  "Stage IA" = "Stage I",
  "Stage IB" = "Stage I",
  "Stage IB1" = "Stage I",
  "Stage IB2" = "Stage I",
  "Stage IC" = "Stage I",
  "Stage I" = "Stage I",

  "Stage IIA" = "Stage II",
  "Stage IIA1" = "Stage II",
  "Stage IIA2" = "Stage II",
  "Stage IIa" = "Stage II",
  "Stage IIB" = "Stage II",
  "Stage IIb" = "Stage II",
  "Stage IIC" = "Stage II",
  "Stage II" = "Stage II",

  "Stage IIIA" = "Stage III",
  "Stage IIIB" = "Stage III",
  "Stage IIIC" = "Stage III",
  "Stage IIIC1" = "Stage III",
  "Stage IIIC2" = "Stage III",
  "Stage III" = "Stage III",

  "Stage IVA" = "Stage IV",
  "Stage IVB" = "Stage IV",
  "Stage IVC" = "Stage IV",
  "Stage IV" = "Stage IV",

  "Stage X" = "Undisclosed",
  "Stage 0" = "Undisclosed"
)

# Check for any stage values not covered by stage_map before trusting the mapping
unmapped <- unique(grouped_data$stage[!grouped_data$stage %in% names(stage_map) &
                                        !(grouped_data$stage == "" | is.na(grouped_data$stage))])
if (length(unmapped) > 0) {
  cat("Warning: unmapped stage values found:\n")
  print(unmapped)
}

grouped_data$stage_grouped <- stage_map[grouped_data$stage]

# Catch both "" and NA stage values as "Undisclosed"
grouped_data$stage_grouped[grouped_data$stage == "" | is.na(grouped_data$stage)] <- "Undisclosed"

table(grouped_data$stage_grouped, useNA = "always")

# Ordinal encoding of stage, for any downstream analysis that needs it as a number
stage_order <- c(
  "Stage I" = 1,
  "Stage II" = 2,
  "Stage III" = 3,
  "Stage IV" = 4
)
grouped_data$stage_ordinal <- stage_order[grouped_data$stage_grouped]

# Drop patients with unreported sex or undisclosed/unmapped stage
grouped_data <- grouped_data[!(grouped_data$sex == "Unreported" |
                                 grouped_data$stage_grouped == "Undisclosed" |
                                 is.na(grouped_data$stage_grouped)), ]

# Patient counts per hospital, used to label the proportion plots below with "n="
n_labels <- aggregate(list(n = grouped_data$site_clean),
                      by = list(site_clean = grouped_data$site_clean),
                      FUN = length)

# Sex distribution by hospital 
# Order hospitals by proportion female (base R, no dplyr)
prop_female_by_site <- tapply(grouped_data$sex == "female", grouped_data$site_clean, mean, na.rm = TRUE)
sex_order <- names(sort(prop_female_by_site))
grouped_data$site_clean <- factor(grouped_data$site_clean, levels = sex_order)

ggplot(grouped_data, aes(x = site_clean, fill = sex)) +
  geom_bar(position = "fill") +
  geom_text(data = n_labels, aes(x = site_clean, y = 1.02, label = paste0("n=", n)),
            inherit.aes = FALSE, size = 2.2) +
  scale_fill_manual(values = c("pink", "lightblue", "grey40")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 8),
        legend.key.size = unit(0.35, "cm")) +
  labs(x = "", y = "Proportion", fill = "Sex",
       title = "Sex Distribution by Hospital (Proportions)") +
  ylim(0, 1.08)

# Race distribution by hospital
# Order hospitals by proportion white (base R, no dplyr)
prop_white_by_site <- tapply(grouped_data$race == "white", grouped_data$site_clean, mean, na.rm = TRUE)
race_order <- names(sort(prop_white_by_site))
grouped_data$site_clean <- factor(grouped_data$site_clean, levels = race_order)

ggplot(grouped_data, aes(x = site_clean, fill = race)) +
  geom_bar(position = "fill") +
  geom_text(data = n_labels, aes(x = site_clean, y = 1.02, label = paste0("n=", n)),
            inherit.aes = FALSE, size = 2.2) +
  scale_fill_manual(values = c("steelblue", "lightskyblue","midnightblue", "powderblue","cadetblue" )) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 8),
        legend.key.size = unit(0.35, "cm")) +
  labs(x = "", y = "Proportion", fill = "Race",
       title = "Ethnicity Distribution by Hospital (Proportions)") +
  ylim(0, 1.08)

# Stage distribution by hospital
ggplot(grouped_data, aes(x = site_clean, fill = stage_grouped)) +
  geom_bar(position = "fill") +
  geom_text(data = n_labels, aes(x = site_clean, y = 1.02, label = paste0("n=", n)),
            inherit.aes = FALSE, size = 2.2) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 8),
        legend.key.size = unit(0.35, "cm")) +
  labs(x = "", y = "Proportion", fill = "Stage",
       title = "Overall Stage Distribution by Hospital") +
  ylim(0, 1.08)

# Cancer type composition by hospital
cancer_colors <- c(
  "Breast Invasive Carcinoma" = "steelblue",

  "Colon Adenocarcinoma" = "salmon",
  "Rectum Adenocarcinoma" = "olivedrab",
  "Stomach Adenocarcinoma" = "sienna",

  "Head and Neck Squamous Cell Carcinoma" = "firebrick",

  "Kidney Chromophobe" = "lightpink",
  "Kidney Renal Clear Cell Carcinoma" = "hotpink",
  "Kidney Renal Papillary Cell Carcinoma" = "lightblue",

  "Lung Squamous Cell Carcinoma" = "darkorchid",
  "Lung Adenocarcinoma" = "turquoise",

  "Pancreatic Adenocarcinoma" = "lightgreen",

  "Liver Hepatocellular Carcinoma" = "peru",

  "Thyroid Carcinoma" = "gold",

  "Ovarian Serous Cystadenocarcinoma" = "orchid",
  "Uterine Corpus Endometrial Carcinoma" = "mediumpurple"
)

# Keep the 11 most common cancer types as their own category; lump the rest into "Other"
top_n <- 11
top_cancers <- names(sort(table(grouped_data$cancer_type), decreasing = TRUE))[1:top_n]

grouped_data$cancer_type_grouped <- ifelse(
  grouped_data$cancer_type %in% top_cancers,
  grouped_data$cancer_type,
  "Other"
)

ggplot(grouped_data, aes(x = site_clean, fill = cancer_type_grouped)) +
  geom_bar(position = "fill") +
  geom_text(data = n_labels, aes(x = site_clean, y = 1.02, label = paste0("n=", n)),
            inherit.aes = FALSE, size = 2.2) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 8),
        legend.key.size = unit(0.35, "cm")) +
  labs(x = "", y = "Proportion", fill = "Cancer Type",
       title = "Cancer Type Composition by Hospital") +
  ylim(0, 1.08)

# Outlier detection (age/year of diagnosis) via the IQR rule
# Flags values more than k*IQR below Q1 or above Q3 for a numeric vector
detect_outliers_iqr <- function(x, k = 1.5) {
  qs <- quantile(x, probs = c(0.25, 0.75), na.rm = TRUE)
  iqr <- diff(qs)
  lower <- qs[1] - k * iqr
  upper <- qs[2] + k * iqr
  which(x < lower | x > upper)
}

# Age outliers, computed separately within each hospital
outlier_flags <- tapply(seq_len(nrow(grouped_data)), grouped_data$site_clean, function(idx) {
  age_out <- detect_outliers_iqr(grouped_data$age_at_diagnosis[idx])
  idx[age_out]
})
outlier_rows_age <- unlist(outlier_flags)
grouped_data[outlier_rows_age, c("site_clean", "cancer_type", "age_at_diagnosis")]

# Year-of-diagnosis outliers, computed separately within each hospital
outlier_flags_year <- tapply(seq_len(nrow(grouped_data)), grouped_data$site_clean, function(idx) {
  year_out <- detect_outliers_iqr(grouped_data$year_of_diagnosis[idx])
  idx[year_out]
})
outlier_rows_year <- unlist(outlier_flags_year)
grouped_data[outlier_rows_year, c("site_clean", "cancer_type", "year_of_diagnosis")]

# Age outliers, computed separately within each (hospital, cancer type) group,
# skipping any group with fewer than 5 patients (too small for a stable IQR)
outlier_flags_age_by_cancer <- tapply(seq_len(nrow(grouped_data)),
                         list(grouped_data$site_clean, grouped_data$cancer_type),
                         function(idx) {
                           if (length(idx) < 5) return(NULL)  # skip tiny groups
                           age_out <- detect_outliers_iqr(grouped_data$age_at_diagnosis[idx])
                           idx[age_out]
                         })
outlier_rows_age_by_cancer <- unlist(outlier_flags_age_by_cancer)
