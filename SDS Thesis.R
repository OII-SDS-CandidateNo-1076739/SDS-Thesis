
# Load Required Libraries
library(tidyverse)  # Includes ggplot2, dplyr, tidyr, and more
library(MASS)       # GLM Negative Binomial
library(tscount)    # Time series count data
library(plm)        # Panel data regression
library(Hmisc)      # Harrell Miscellaneous
library(zoo)        # For RowMean replace
library(mgcv)       # Generalized Additive Models
library(stargazer)  # Plot tables in LaTeX
library(RColorBrewer) # Color palettes
library(scales)     # Graphical scales
library(magrittr)   # Pipe operators
library(gridExtra)  # Arrange ggplot objects
library(GGally)     # Extension to ggplot
library(MuMIn)      # Multi-model inference
library(glmmTMB)    # Generalized linear mixed models
library(lme4)       # Linear mixed-effects models
library(lmerTest)   # Tests for random and fixed effects
library(lmtest)     # Testing linear regression models
library(sandwich)   # Robust covariance matrix estimators
library(nlme)       # Linear and Nonlinear Mixed Effects Models
library(tseries)    # Time-series analysis
library(reshape2)   # Data reshaping
library(ggrepel)    # Automatically position ggplot labels
library(purrr)      # Functional programming tools
library(flextable)  # Table formatting
library(itsadug)    # Tools for visualising mixed-effects models
library(randomForest) # Random Forests
library(missForest)  # Nonparametric missing value imputation
library(gratia)     # ggplot-based graphics for GAMs
library(patchwork)   # Arrange ggplot objects
library(mgcViz)     # Visualisations for GAMs
library(DHARMa)     # Residual diagnostics for hierarchical models




# --------------------------- Data Reprocessing ---------------------------.

# Step 1 - Time series Value Imputation
data_simple_imputed <- df_orig %>%
  group_by(NUTS_3) %>%
  mutate(Int_Access_missing = is.na(Int_Access) & !is.na(dplyr::lead(Int_Access))) %>%
  fill(Int_Access, .direction = "down") %>%
  ungroup() %>%
  select(-Int_Access_missing)

# Step 2 - Country unconditional average Impution
data_country_imputed <- data_simple_imputed %>%
  group_by(COUNTRY) %>%
  mutate(Int_Access = ifelse(is.na(Int_Access), mean(Int_Access, na.rm = TRUE), Int_Access)) %>%
  ungroup()

# Step 3 - MICE Imputation
data_names_fixed <- df_orig
names(data_names_fixed)[names(data_names_fixed) == "Int_Access"] <- "Int_Access"
mice_imputed <- mice(data_names_fixed, m=1, maxit=50, method='cart', seed=500)
data_mice_imputed <- complete(mice_imputed, 1)

# Apply Imputed Values to Original Dataframe
df_orig$Int_Access <- data_mice_imputed$Int_Access

# Create a DataFrame with No Missing Values
df_no_missing <- df_orig %>%
  mutate(all_vars = rowSums(!is.na(.)) == ncol(.)) %>%
  filter(all_vars)

# --------------------------- Data Collection ---------------------------.

# Read Project Collection Data
project_collection <- read_csv("data/final_platform_data.csv")
project_collection$project_start_date <- as.Date(project_collection$project_start_date, format="%d/%m/%Y")
project_collection$year_month <- format(project_collection$project_start_date, "%Y-%m")
count_by_year_month <- project_collection %>%
  group_by(year_month) %>%
  summarise(count = n()) %>%
  arrange(year_month)




# --------------------------- Data Visualization ---------------------------.



###### Data collection Types Visualization ######.

project_collection <- read_csv("data/final_platform_data.csv")

# Convert 'project_start_date' to Date format with the correct format
project_collection$project_start_date <- as.Date(project_collection$project_start_date, format="%d/%m/%Y")

# Check for missing or malformed dates again
malformed_dates <- which(is.na(project_collection$project_start_date))
print(paste("Number of malformed dates:", length(malformed_dates)))

# Assuming no malformed dates, proceed to create the 'year_month' and 'count_by_year_month'
project_collection$year_month <- format(project_collection$project_start_date, "%Y-%m")

count_by_year_month <- project_collection %>% 
  group_by(year_month) %>% 
  summarise(count = n()) %>% 
  arrange(year_month)

# Convert 'year_month' back to character for plotting
count_by_year_month$year_month <- as.character(count_by_year_month$year_month)

# Calculate the number of green bars
green_bars_count <- sum(as.Date(paste0(count_by_year_month$year_month, "-01"), format="%Y-%m-%d") <= as.Date("2016-03-01"))

# Define bar colors based on conditions
bar_colors_updated <- c(rep("#77b9cf", 6), 
                        rep("#348fe8", green_bars_count - 6), 
                        rep("#3142aa", nrow(count_by_year_month) - green_bars_count))

# Find the first instance of each year to use as breaks
year_breaks <- unique(substr(count_by_year_month$year_month, 1, 4))
year_breaks <- sapply(year_breaks, function(y) paste0(y, "-01"))

# Generate the barplot
projects_collection_plot <- ggplot(count_by_year_month, aes(x=factor(year_month), y=count, fill=factor(bar_colors_updated))) +
  geom_bar(stat="identity") +
  scale_fill_manual(values=c("#77b9cf", "#348fe8", "#3142aa"), 
                    labels=c("Projects Provided Directly by the Platform", "Projects Collected via the API – Before the OLI", "Projects Collected via the API – After the OLI"),
                    breaks=c("#77b9cf", "#348fe8", "#3142aa")) +
  labs(title="Monthly Distribution of Projects by Collection Method",
       x="Year",
       y="Number of Projects",
       fill=NULL) +
  theme_minimal() +
  theme(legend.position = c(0.3, 0.8)) +  # Adjust legend position slightly down
  scale_x_discrete(breaks = year_breaks, labels = unique(substr(year_breaks, 1, 4)))  # Show only full years on the x-axis

projects_collection_plot

ggsave(filename = "correlation_plot_collection.png", plot = projects_collection_plot, width = 12, height = 6, dpi = 300, bg = "white")








##### Visualizing Original Distributions #####.

# Define visualization parameters
bar_fill_color <- "#77b9cf"
bar_border_color <- "white"
density_line_color <- "#3142aa"
bar_alpha <- 0.5

# Custom function to format numbers as millions or thousands
format_millions <- function(x) {
  if (x >= 1e6) {
    return(paste0(round(x / 1e6, 1), "M"))
  } else {
    return(paste0(round(x / 1e3, 1), "K"))
  }
}

# Custom breaks and labels for Population and GDP
breaks_pop <- seq(0, max(df_orig$Population, na.rm = TRUE), by = 2e6)
labels_pop <- sapply(breaks_pop, format_millions)
breaks_gdp <- seq(0, max(df_orig$GDP, na.rm = TRUE), by = 5e4)
labels_gdp <- sapply(breaks_gdp, format_millions)

## Plot1 - Population
plot1 <- ggplot(df_orig, aes(x = Population)) +
  geom_histogram(aes(y = ..density..), bins = 30, fill = bar_fill_color, color = bar_border_color, alpha = bar_alpha) +
  geom_density(colour = density_line_color, adjust = 2) +
  labs(x = "Population", y = "Density") +
  scale_x_continuous(breaks = breaks_pop, labels = labels_pop) +
  scale_y_continuous(labels = scales::comma) +  # Add this line
  theme_minimal()
## Plot2 - Int_Access
plot2 <- ggplot(df_orig, aes(x = Int_Access)) +
  geom_histogram(aes(y = ..density..), bins = 30, fill = bar_fill_color, color = bar_border_color, alpha = bar_alpha) +
  geom_density(colour = density_line_color, adjust = 2) +
  labs(x = "Int_Access", y = "Density") +
  theme_minimal()
## Plot3 - GDP
plot3 <- ggplot(df_orig, aes(x = GDP)) +
  geom_histogram(aes(y = ..density..), bins = 30, fill = bar_fill_color, color = bar_border_color, alpha = bar_alpha) +
  geom_density(colour = density_line_color, adjust = 2) +
  labs(x = "GDP", y = "Density") +
  scale_x_continuous(breaks = breaks_gdp, labels = labels_gdp) +
  scale_y_continuous(labels = scales::comma) +  # Add this line
  theme_minimal() 
## Plot4 - Rural
plot4 <- ggplot(df_orig, aes(x = factor(Rural))) +
  geom_bar(fill = bar_fill_color, color = bar_border_color, alpha = bar_alpha) +
  labs(x = "Urban/Rural", y = "Count") +
  theme_minimal() 
## Plot5 - Unemp
plot5 <- ggplot(df_orig, aes(x = Unemp)) +
  geom_histogram(aes(y = ..density..), bins = 30, fill = bar_fill_color, color = bar_border_color, alpha = bar_alpha) +
  geom_density(colour = density_line_color, adjust = 2) +
  labs(x = "Unemploument", y = "Density") +
  theme_minimal()
## Plot6 - Education
plot6 <- ggplot(df_orig, aes(x = Education)) +
  geom_histogram(aes(y = ..density..), bins = 30, fill = bar_fill_color, color = bar_border_color, alpha = bar_alpha) +
  geom_density(colour = density_line_color, adjust = 2) +
  labs(x = "Education", y = "Density") +
  theme_minimal()
## Plot7 - Projects
plot7 <- ggplot(df_orig, aes(x = Projects)) +
  geom_histogram(aes(y = ..density..), bins = 30, fill = bar_fill_color, color = bar_border_color, alpha = bar_alpha) +
  geom_density(colour = density_line_color, adjust = 2) +
  labs(x = "Projects", y = "Density") +
  theme_minimal()
## Combine and display the plots with a main title
combined_plots <- grid.arrange(plot1, plot2, plot3, plot4, plot5, plot6, plot7, ncol = 3, top = "Original Distributions")
print(combined_plots)
ggsave(filename = "distribution_orig.png", plot = combined_plots, width = 10, height = 6, dpi = 300, bg = "white")





##### Visualizing Transformed Distributions ######.

## Plot1 - Population
plot1 <- ggplot(df_orig, aes(x = log(Population))) +
  geom_histogram(aes(y = ..density..), bins = 30, fill = bar_fill_color, color = bar_border_color, alpha = bar_alpha) +
  geom_density(colour = density_line_color, adjust = 2) +
  labs(x = "Log Population", y = "Density") +
  theme_minimal() 
## Plot2 - Int_Access
plot2 <- ggplot(df_orig, aes(x = Int_Access)) +
  geom_histogram(aes(y = ..density..), bins = 30, fill = bar_fill_color, color = bar_border_color, alpha = bar_alpha) +
  geom_density(colour = density_line_color, adjust = 2) +
  labs(x = "Int_Access", y = "Density") +
  theme_minimal()
## Plot3 - GDP
plot3 <- ggplot(df_orig, aes(x = log(GDP))) +
  geom_histogram(aes(y = ..density..), bins = 30, fill = bar_fill_color, color = bar_border_color, alpha = bar_alpha) +
  geom_density(colour = density_line_color, adjust = 2) +
  labs(x = "Log GDP", y = "Density") +
  theme_minimal() 
## Plot4 - Rural
plot4 <- ggplot(df_orig, aes(x = factor(Rural))) +
  geom_bar(fill = bar_fill_color, color = bar_border_color, alpha = bar_alpha) +
  labs(x = "Urban/Rural", y = "Count") +
  theme_minimal() 
## Plot5 - Unemp
plot5 <- ggplot(df_orig, aes(x = log(Unemp))) +
  geom_histogram(aes(y = ..density..), bins = 30, fill = bar_fill_color, color = bar_border_color, alpha = bar_alpha) +
  geom_density(colour = density_line_color, adjust = 2) +
  labs(x = "Log Unemploument", y = "Density") +
  theme_minimal()
## Plot6 - Education
plot6 <- ggplot(df_orig, aes(x = Education)) +
  geom_histogram(aes(y = ..density..), bins = 30, fill = bar_fill_color, color = bar_border_color, alpha = bar_alpha) +
  geom_density(colour = density_line_color, adjust = 2) +
  labs(x = "Education", y = "Density") +
  theme_minimal()
## Plot7 - Projects
plot7 <- ggplot(df_orig, aes(x = log(Projects))) +
  geom_histogram(aes(y = ..density..), bins = 30, fill = bar_fill_color, color = bar_border_color, alpha = bar_alpha) +
  geom_density(colour = density_line_color, adjust = 2) +
  labs(x = "Log Projects", y = "Density") +
  theme_minimal()
## Plot8 - Projects
plot8 <- ggplot(df_orig, aes(x = log1p(Projects))) +
  geom_histogram(aes(y = ..density..), bins = 30, fill = bar_fill_color, color = bar_border_color, alpha = bar_alpha) +
  geom_density(colour = density_line_color, adjust = 2) +
  labs(x = "Log1p Projects", y = "Density") +
  theme_minimal()
## Plot9 - Projects
plot9 <- ggplot(df_orig, aes(x = ihs(Projects))) +
  geom_histogram(aes(y = ..density..), bins = 30, fill = bar_fill_color, color = bar_border_color, alpha = bar_alpha) +
  geom_density(colour = density_line_color, adjust = 2) +
  labs(x = "IHS Projects", y = "Density") +
  theme_minimal()
## Combine and display the plots with a main title
combined_plots <- grid.arrange(plot1, plot2, plot3, plot4, plot5, plot6, plot7, plot8, plot9, ncol = 3, top = "Transformed Distributions")
print(combined_plots)
ggsave(filename = "distribution_transf.png", plot = combined_plots, width = 10, height = 6, dpi = 300, bg = "white")








##### Distribution and Dispersion of the Dependent: Project Counts ######.

# Distribution and Dispersion of Project Counts
## Check data dispersion ( we have overdispersion -> mean << variance)
dist_mean <- round(mean(df_orig$Projects, na.rm=T),2) # 18.94
dist_median<- median(df_orig$Projects, na.rm=T) # 3
dist_var <- round(var(df_orig$Projects, na.rm=T),0) # 3896
dist_zero <- round(length(df_orig$Projects[df_orig$Projects == 0])/length(df_orig$Projects == 0)*100,2) #41.23

## Distribution with mean and median
dist1 <- df_orig %>% ggplot(aes(x=Projects)) + geom_density(alpha=.25, fill="#df2644") + 
  geom_vline(aes(xintercept=mean(Projects, na.rm=T)), color="#a80000", linetype="dashed", linewidth=0.7) +
  geom_vline(aes(xintercept=median(Projects, na.rm=T)), color="#0a63bf", linetype="dashed", linewidth=0.7) + 
  theme_light() + 
  annotate(geom="text", x=1400, y=0.05, label = paste("Median: ",dist_median),color="#0a63bf") + 
  annotate(geom="text", x=1400, y=0.045, label=paste("Mean: ",dist_mean),color="#a80000") +
  annotate(geom="text", x=1400, y=0.04, label=paste("Variance: ",dist_var),color="black") +
  annotate(geom="text", x=1400, y=0.035, label=paste("% Zeros: ",dist_zero),color="black")

## Distribution on log scale
dist2 <- df_orig %>% ggplot(aes(x=Projects)) + geom_density(alpha=.25, fill="#df2644") + 
  geom_vline(aes(xintercept=mean(Projects, na.rm=T)), color="#a80000", linetype="dashed", linewidth=0.7) +
  geom_vline(aes(xintercept=median(Projects, na.rm=T)), color="#0a63bf", linetype="dashed", linewidth=0.7) + 
  theme_light() + scale_x_log10()

## Distribution as boxplot
dist3 <- df_orig %>% ggplot(aes(x=Projects)) + geom_boxplot() + scale_x_continuous(trans='log2') + coord_flip() +
  theme_light() + theme(axis.text.x=element_blank())

## Arrange all distribution charts
final_dist <- ggarrange(dist1, ggarrange(dist2, dist3, ncol = 2, labels = c("B", "C")), nrow = 2,labels = "A") 







##### VIOLIN PLOTS of Socioeconomic Variables By Region Type #######.

# Define a colour palette
colour_palette <- c("1" = "#df2644", "2" = "#ffd631", "3" = "#1898c2")
# Function to set the labels for the 'Rural' factor
set_rural_labels <- function() {
  scale_fill_manual(values = colour_palette, labels = c("Urban", "Intermediate", "Rural"))
}

## Prepare Data
### 1. Select necessary columns and process data
violin_subset <- df_orig %>% 
  dplyr::select(NUTS_3, Projects, Rural, Population, `Int_Access`, Unemp, Education, GDP) %>% 
  dplyr::group_by(NUTS_3) %>%
  dplyr::summarise(
    Population = median(Population, na.rm = TRUE),
    Projects = sum(Projects, na.rm = TRUE),
    `Int_Access` = median(`Int_Access`, na.rm = TRUE),
    Unemp = median(Unemp, na.rm = TRUE),
    Education = median(Education, na.rm = TRUE),
    GDP = median(GDP, na.rm = TRUE),
    Rural = first(Rural) # Take the first Rural value for each group
  ) %>%
  dplyr::mutate(projects_per_10000 = Projects / Population * 10000) # Create a variable for number of projects per 10,000 people

### 2. Handle factors and missing values
violin_subset$Rural <- factor(violin_subset$Rural, levels = c(1, 2, 3), labels = c('Urban', 'Intermediate', 'Rural'))
violin_subset <- na.omit(violin_subset) # Remove rows with NA values

### 3. Convert data to long format
violin_long <- violin_subset %>% gather(key = "Variable", value = "Value", -NUTS_3, -Rural, -Population, -Projects)

### 4. Standardize data
violin_long <- violin_long %>% group_by(Variable) %>% mutate(ValueStd = scale(Value))

## Plotting
violin_plot <- ggplot(violin_long, aes(x = Rural, y = Value, fill = Rural)) +
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.4, method = "histodot", colour = "white") +
  geom_boxplot(width = 0.5, alpha = 0, outlier.shape = NA, lwd = 0.6) +
  scale_y_log10() +
  labs(x = "Type of Region", y = "Metrics Value (Log Scale)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, hjust = 1), legend.position = "none") +
  scale_color_brewer(palette = "Dark2") +
  facet_wrap(~ Variable, scales = "free_y", ncol = 3) # Create separate panels for each variable
print(violin_plot)








###### TIME SERIES PLOTS of Socioeconomic Variables By Region Type #######.

## Prepare Data
### 1. Aggregate and process data for projects
data_agg_pr <- aggregate(cbind(Projects, Population) ~ year + Rural,  
                         data = df_orig,
                         FUN = sum, 
                         na.rm = TRUE)
names(data_agg_pr)[3] <- "sum_projects" 
names(data_agg_pr)[4] <- "sum_population"
data_agg_pr <- data_agg_pr %>%
  mutate(projects_per_10000 = (sum_projects / sum_population) * 10000)
totals <- data_agg_pr %>%
  group_by(year) %>%
  summarise(total_projects = sum(projects_per_10000, na.rm = TRUE))
data_agg_pr <- left_join(data_agg_pr, totals, by = "year") %>%
  mutate(prop_projects = projects_per_10000 / total_projects)


### 2. Set titles for the plot
titles <- c(
  "Unemp" = "Unemployment Rate",
  "Int_Access" = "Internet Access Rate",
  "Education" = "Education Level",
  "Population" = "Population Count in thousands",
  "GDP" = "GDP per capita",
  "prop_projects" = "Proportion of Projects per 10k people"
)

### 3. Aggregate other variables
data_agg_temp <- list(
  Unemp = df_orig %>% aggregate(Unemp ~ year + Rural, ., median),
  `Int_Access` = df_orig %>% aggregate(`Int_Access` ~ year + Rural, ., median),
  Education = df_orig %>% aggregate(Education ~ year + Rural, ., median),
  Population = df_orig %>% aggregate(Population ~ year + Rural, ., median),
  GDP = df_orig %>% aggregate(GDP ~ year + Rural, ., median)
) %>%
  Reduce(function(x, y) merge(x, y, by=c("year", "Rural")), .)

data_agg_temp$Population <- data_agg_temp$Population / 1e3 # Adjust the population to 10k units

data_agg <- data_agg_temp %>%
  left_join(data_agg_pr, by = c("year", "Rural")) %>%
  gather(key = "Variable", value = "Value", -year, -Rural) %>%
  filter(Variable %in% names(titles))

## Plotting
trends_plot <- ggplot(data_agg, aes(x = year, y = Value, color = factor(Rural))) +
  geom_line(size = 0.7, alpha = 0.7) +
  scale_color_manual(values = colour_palette, labels = c("Urban", "Intermediate", "Rural")) +  # Updated this line to use your custom colour palette
  labs(
    title = "Timeline: Various Metrics",
    x = "Year", 
    y = "Median Value",
    color = "Region Type"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "bottom",
    axis.text.x = element_text(angle = 0),
    panel.grid.major = element_line(color = "grey92"),
    panel.grid.minor = element_blank(),
    strip.text = element_text(size = 12, face = "bold")
  ) +
  scale_x_continuous(breaks = unique(data_agg$year)) +
  facet_wrap(~ Variable, scales = "free_y", ncol = 3, labeller = as_labeller(titles))

print(trends_plot)







####### PROJECTS Disparities by Region Type and Country 2020 ######.

# Filter out NaN in the Rural column
dataset <- subset(df_orig, !is.na(Rural))

# Filter the data for the year 2020
dataset_2020 <- subset(dataset, year == 2020)

# Aggregate the data by 'COUNTRY' and 'Rural', summing 'Projects' and 'Population'
grouped_data <- dataset_2020 %>% 
  group_by(COUNTRY, Rural) %>%
  summarise(Projects = sum(Projects), Population = sum(Population)) %>%
  mutate(Projects_per_million = (Projects / Population) * 1e6)

# Calculate the median project count per 1 million population for each country
median_project_count <- grouped_data %>% 
  group_by(COUNTRY) %>%
  summarise(Median_Projects = median(Projects_per_million, na.rm = TRUE)) %>%
  arrange(Median_Projects)

# Merge the median project count back into the original grouped data
grouped_data <- merge(grouped_data, median_project_count, by = "COUNTRY")

# Create the bubble plot with adjustments
final_plot <- ggplot(grouped_data, aes(x = reorder(COUNTRY, Median_Projects), y = Projects_per_million, fill = factor(Rural), size = Projects)) +
  geom_jitter(alpha = 0.6, width = 0.3, shape = 21, color = "black", stroke = 1) +
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  annotation_logticks(sides = "l") +
  scale_size_continuous(breaks = c(1, 10, 100, 250, 500, 1000), range = c(1, 10)) +
  set_rural_labels() +
  guides(fill = guide_legend(override.aes = list(size = 5))) +
  labs(title = "Projects Count Per 1 Million Population by Country and Area Type (2020)",
       x = "Country",
       y = "Projects Count Per 1 Million Population",
       fill = "Region Type",
       size = "Total Projects") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.position = "right",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank())  # Remove dotted horizontal lines

# Show the plot
print(final_plot)

ggsave(filename = "project_circles.png", plot = final_plot, width = 10, height = 4, dpi = 300, bg = "white")









######## Imbalance between Urban/Rural Regions in OLM #########

# Read the full dataset
all_data <- df_orig

# Aggregate the data by 'COUNTRY' and 'Rural', summing 'Projects' and 'Population'
grouped_all_data <- all_data %>%
  group_by(COUNTRY, Rural) %>%
  summarise(Projects = sum(Projects), Population = sum(Population)) %>%
  mutate(Projects_per_million = (Projects / Population) * 1e6)

# Count the total number of projects for each country
total_projects <- all_data %>%
  group_by(COUNTRY) %>%
  summarise(Total_Projects = sum(Projects))

# Pivot the data to wide format
wide_all_data <- grouped_all_data %>%
  pivot_wider(id_cols = COUNTRY, names_from = Rural, values_from = Projects_per_million, names_prefix = "Area_")

# Calculate the imbalance ratio excluding intermediate areas
wide_all_data <- wide_all_data %>%
  mutate(Imbalance_Ratio = Area_1 / Area_3) %>%
  arrange(-Imbalance_Ratio)

# Merge the total number of projects back to the wide_all_data
wide_all_data <- merge(wide_all_data, total_projects, by = "COUNTRY")

# Create the bar plot with gradient colour
imbalance_plot <- ggplot(wide_all_data, aes(x = reorder(COUNTRY, -Imbalance_Ratio), y = Imbalance_Ratio)) +
  geom_bar(stat = "identity", aes(fill = Imbalance_Ratio)) +
  scale_fill_gradient(low = "#E68A9B", high = "#E0435F") +
  geom_text(aes(label = formatC(Total_Projects, big.mark = ",", format = "d")), 
            angle = 90, y=0.4, hjust=0.2, vjust=0.3, color = "#611D2A") +
  labs(title = "Online Labour Markets Imbalance Between Urban and Rural Areas in Online Labour Markets (2013-2020)",
       x = "Country",
       y = "Imbalance Ratio",
       fill = "Imbalance Ratio") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.position = "none")

# Show the plot
print(imbalance_plot)

ggsave(filename = "imbalance_plot.png", plot = imbalance_plot, width = 10, height = 4, dpi = 300, bg = "white")

wide_all_data$COUNTRY





###### Projects Imbalance Map #######.

# Group, Summarize, and Calculate Metrics 
data <- df_orig %>%
  group_by(NUTS_3, COUNTRY) %>%  # Group by both NUTS_3 and COUNTRY
  summarise(
    Projects = sum(Projects, na.rm = TRUE),
    Rural = first(Rural, order_by = NULL),
    Population = mean(Population, na.rm = TRUE)
  ) %>%
  mutate(
    projects_per_10k = (Projects / Population) * 10000,
    projects_per_10k_log = log1p(projects_per_10k)  # log(1+x)
  ) %>%
  na.omit()  # Remove all NAs


# Calculate the imbalance ratio for each country
country_imbalance_ratio <- data %>%
  group_by(COUNTRY, Rural) %>%
  summarise(mean_projects_per_10k = mean(projects_per_10k, na.rm = TRUE)) %>%
  spread(Rural, mean_projects_per_10k) %>%
  mutate(imbalance_ratio = `1` / `3`) %>%
  na.omit() # Remove NAs

# Join with country geospatial data, assuming "ADMIN" is the country name column in country_geo
country_geo_imbal <- country_geo %>%
  left_join(country_imbalance_ratio, by = c("ISO_A2_EH" = "COUNTRY"))

# Filter country_geo_imbal to only include countries present in country_imbalance_ratio
country_geo_imbal <- country_geo_imbal %>% 
  filter(ISO_A2_EH %in% country_imbalance_ratio$COUNTRY)



# Custom colour palette
custom_palette <- rev(c("#e0435f", "#e65461", "#ec6463", "#f88567", "#ece1b4", "#e0e0cd"))  # Reverse the palette

# Your plotting code, with two geom_sf layers
imbalance_map <- ggplot() +
  
  # 1. Draw the relief
  geom_raster(data = relief, inherit.aes = FALSE, aes(x = x, y = y, alpha = HYP_50M_SR_W_1)) +
  scale_alpha(name = "", range = c(0.6, 0), guide = F) +
  
  # 2. Draw all country borders, fill set to NA for no fill
  geom_sf(data = country_geo, fill = NA, color = "white", lwd = 0.5) +
  
  # 3. Draw the choropleth for countries of interest, colored by imbalance_ratio
  geom_sf(data = country_geo_imbal, aes(fill = imbalance_ratio), color = "white", lwd = 0.5) +
  
  # Aesthetics & Scales
  scale_fill_gradientn(
    colors = scales::alpha(custom_palette, 0.8),  # Custom palette with alpha for transparency
    name = "Imbalance Ratio",
    guide = guide_colourbar(
      title.position = "top",
      reverse = FALSE)  # Set to FALSE to not reverse the legend
  ) +
  
  # Labels
  labs(
    x = NULL, y = NULL,
    title = "Imbalance Ratio of Online Projects in European Countries",
    subtitle = "Ratio of Projects per 10k People in Urban vs Rural Areas",
    caption = default_caption
  ) +
  
  # Theme & Coordinates
  theme_map() +
  coord_sf(xlim = c(-25, 40), ylim = c(34, 71))  # Trim to focus on Europe

# Display the plot
print(imbalance_map)

# Save the plot to a PDF
ggsave("map_imbalance_ratio_custom_palette_reversed.pdf", imbalance_map, width = 8, height = 6, bg = default_background_color, device = cairo_pdf)














#----------------------------- Regression Models -----------------------------.


# Stepwise model selection (all vars are significant)
## Fit a full model first
full.model <- lm(ihs(Projects) ~ log(Population) + Int_Access + log(GDP) 
                 + factor(Rural) + log(Unemp) + Education, data = df_orig
                 %>% filter(all_vars == 1))

## Run stepwise regression
stepwise.model <- step(full.model, direction = "both")

## Print the summary of the final model
summary(stepwise.model)


## Generate the HTML table using tab_model
tab_model(stepwise.model, 
          p.style = c("stars"),
          file = "Stepwise_Model_Summary.html")



## OLS log1p
summary(m_ols_log1p <- lm(ihs(Projects) ~ log(Population) + factor(Rural) + Int_Access + log(GDP) + Education + log(Unemp), 
                          data = df_no_missing %>% filter(all_vars == 1)))

## OLS IHS
summary(m_ols_ihs <- lm(ihs(Projects) ~ log(Population) + factor(Rural) + Int_Access + log(GDP) + Education + log(Unemp),
                        data=df_no_missing %>% filter(all_vars==1)))

## OLS log1p with Country and Year Fixed Effects
summary(m_ols_log1p_fe <- plm(ihs(Projects) ~ log(Population) + factor(Rural) + Int_Access + log(GDP) + Education + log(Unemp),
                              data = df_no_missing %>% filter(all_vars == 1),
                              index = c("COUNTRY", "year"),
                              model = "within",
                              effect = "twoways",
                              effect.id = TRUE))

#### OLS log1p with Country and Year Random Effects: This is a type of panel model that considers both within-group 
summary(m_ols_log1p_re <- plm(ihs(Projects) ~ log(Population) + factor(Rural) + Int_Access + log(GDP) + Education + log(Unemp),
                              data = df_orig %>% filter(all_vars == 1),
                              index = c("COUNTRY", "year"),
                              model = "random",
                              random.method = "walhus",
                              effect = "twoways"))

#### OLS LOG with NAs instead of 0
summary(m_ols_log <- lm(ihs(Projects.NA) ~  log(Population) + factor(Rural) + Int_Access + log(GDP) + Education + log(Unemp),
                        data=df_orig %>% filter(all_vars==1)))


#### Test Fixed Effects versus OLS
pFtest(m_ols_log1p_fe, m_ols_log1p)
#### Test Fixed versus Random Effects
phtest(m_ols_log1p_fe, m_ols_log1p_re)

#### GLM Negative Binomial
summary(m_nbin <- glm.nb(Projects ~ log(Population) + factor(Rural) + Int_Access + log(GDP) + Education + log(Unemp),
                         data=df_orig %>% filter(all_vars==1)))


# Nested Mixed Effect Model
m_mult_mix_log1p <- lmer(ihs(Projects) ~ log(Population) + factor(Rural) + Int_Access + 
                           log(GDP) + Education + log(Unemp) +
                           (1 | COUNTRY) + (1 | COUNTRY:NUTS_3) + factor(year), 
                         data = df_no_missing)

summary(m_mult_mix_log1p)
AIC(m_mult_mix_log1p)
BIC(m_mult_mix_log1p)
r.squaredGLMM(m_mult_mix_log1p)



# Negative Binomial Nested Mixed Effect Model
m_nb_fe <- glmmTMB(Projects ~ log(Population) + factor(Rural) + Int_Access + 
                     log(GDP) + Education + log(Unemp) +
                     (1 | COUNTRY) + (1 | COUNTRY:NUTS_3) + (1 | year), 
                   data = df_no_missing, family = nbinom2)

summary(m_nb_fe)
AIC(m_nb_fe)
BIC(m_nb_fe)
logLik(m_nb_fe)
r.squaredGLMM(m_nb_fe)

# GAM Model
m_gam <- gam(log1p(Projects) ~ s(Population, bs = "ps") + 
               factor(Rural) + 
               s(Int_Access, bs = "ps") +
               s(GDP, bs = "ps") + 
               s(Education, bs = "ps") +
               s(Unemp, bs = "ps") + 
               factor(COUNTRY) + 
               factor(year),  
             data = df_no_missing,
             method = "REML")





#-------------------------- Regression Results ---------------------------------.



###### Main Stat Table ######.
tab_model(m_ols_ihs, m_ols_log1p, m_ols_log1p_fe, m_nb_fe, m_mult_mix_log1p, m_mult_mix_log1p, 
          file = "model_comparison_impute_final_no_imputation.html", 
          transform = NULL, # show original coefficients instead of Incidence Rate Ratios (negative binomial)
          p.style = c("stars"),
          show.ci = FALSE,
          pred.labels = c("(Intercept)", "Population [log]", "Urb/Rur: Intermediate", "Urb/Rur: Rural",
                          "Int_Access", "GDP per cap [log]", "Education", "Unemployment [log]"),
          dv.labels = c("OLS<br>(IHS)<br>(1)", 
                        "OLS<br>(Log+1)<br>(2)",
                        "OLS FE<br>(Log+1)<br>(3)", 
                        "NB Mixed<br>(No Transf)<br>(4)",
                        "Mixed FE<br>(Log+1)<br>(5)", 
                        "GAM<br>(Log+1)<br>(6)"))




###### Table: MAE	MSE	RMSE Correlation AIC	BIC ###### .

# Prepare a list of models and their names
models <- list(m_ols_log1p, m_ols_log1p_fe, m_ols_log1p_re, m_ols_log, m_nbin, m_mult_mix_log1p, m_nb_fe, m_gam)
model_names <- c("OLS\n(log1p)", "OLS FE\n(log1p)", "OLS RE\n(log1p)",
                 "OLS NAs\n(LOG)", "NB\n(No Transf)", "Mixed FE\n(log1p)",
                 "NB Mixed FE\n(No Transf)", "GAM Model\n(LOG1P)")

# Initialize an empty data frame to store the results
results <- data.frame(Model = character(), MAE = numeric(), MSE = numeric(), RMSE = numeric(), Correlation = numeric(), AIC = numeric(), BIC = numeric())

# Loop over the models
for (i in 1:length(models)) {
  # Get the predicted values
  preds <- predict(models[[i]], newdata = df_no_missing %>% filter(all_vars == 1))
  
  # Calculate the metrics
  actuals <- df_no_missing$Projects.log1p
  mae <- mean(abs(preds - actuals))
  mse <- mean((preds - actuals)^2)
  rmse <- sqrt(mse)
  cor <- cor(preds, actuals)
  aic <- AIC(models[[i]])
  bic <- BIC(models[[i]])
  
  # Add the results to the data frame
  results <- rbind(results, data.frame(Model = model_names[i], MAE = mae, MSE = mse, RMSE = rmse, Correlation = cor, AIC = aic, BIC = bic))
}

# Normalize the measures so they're all on the same scale
results$normalized_MAE <- scale(results$MAE)
results$normalized_MSE <- scale(results$MSE)
results$normalized_RMSE <- scale(results$RMSE)
results$normalized_Corr <- scale(results$Correlation)
results$normalized_AIC <- scale(results$AIC)
results$normalized_BIC <- scale(results$BIC)

# Create a score that is the average of the normalized measures
# For MAE, MSE, RMSE, AIC, and BIC, lower is better so we multiply by -1
# For Correlation, higher is better so we don't multiply by -1
results$Score <- rowMeans(data.frame(
  -results$normalized_MAE,
  -results$normalized_MSE,
  -results$normalized_RMSE,
  results$normalized_Corr,
  -results$normalized_AIC,
  -results$normalized_BIC
))

# Now we can sort by the Score
results <- results[order(results$Score, decreasing = TRUE),]

# Drop the normalized and score columns
results <- subset(results, select = -c(normalized_MAE, normalized_MSE, normalized_RMSE, 
                                       normalized_Corr, normalized_AIC, normalized_BIC, Score))


# Define color vectors based on relative performance
pastel_green <- "#66C776"
pastel_red <- "#E8646F"
pastel_yellow <- "#FFFFFF"

color_RMSE <- color_Corr <- c()

for(i in 1:nrow(results)){
  if(i == 1){
    color_RMSE <- c(color_RMSE, pastel_green) 
    color_Corr <- c(color_Corr, pastel_red) 
  } else if(i == nrow(results)){
    color_RMSE <- c(color_RMSE, pastel_red)
    color_Corr <- c(color_Corr, pastel_green)
  } else {
    color_RMSE <- c(color_RMSE, pastel_yellow)
    color_Corr <- c(color_Corr, pastel_yellow)
  }
}

gt(results) %>%
  fmt_markdown(columns = c("Model")) %>%
  tab_header(title = "Model Performance") %>%
  tab_style(style = cell_fill(color = "white"), 
            locations = cells_body(columns = "Model")) %>%
  cols_align(align = "center", columns = everything()) %>%
  cols_align(align = "left", columns = "Model") %>%
  tab_style(style = cell_borders(sides = "all", color = "gray", weight = px(1)), 
            locations = cells_body(columns = everything())) %>%
  data_color(columns = c("MAE", "MSE", "RMSE"), colors = color_RMSE) %>%
  data_color(columns = c("Correlation"), colors = color_Corr) %>%
  data_color(columns = c("AIC", "BIC"), colors = color_RMSE)






######  Cross-validated MAE / Correlation : 10-fold CV loops ######.

library(tidyverse)
library(rsample)
library(ggplot2)
library(mgcv)

#------------------- 1. Setup and Data Filtering -------------------#
df <- df_orig %>% filter(all_vars == 1)

threshold <- 50 
rare_countries <- df %>%
  group_by(COUNTRY) %>%
  tally() %>%
  filter(n < threshold) %>%
  pull(COUNTRY)

df <- df %>%
  filter(!(COUNTRY %in% rare_countries))

#------------------- 2. Model Definition -------------------#

# List of model formulas
model_formulas <- list(
  ols_ihs = list(
    formula = ihs(Projects) ~ log(Population) + factor(Rural) + Int_Access + log(GDP) + Education + log(Unemp),
    func = lm,
    response_transformation = "ihs"
  ),
  ols_log1p = list(
    formula = log1p(Projects) ~ log(Population) + factor(Rural) + Int_Access + log(GDP) + Education + log(Unemp),
    func = lm,
    response_transformation = "log1p" 
  ),
  m_ols_log1p_fe = list(
    formula = log1p(Projects) ~ log(Population) + factor(Rural) + Int_Access + log(GDP) + Education + log(Unemp),
    func = plm,
    index = c("COUNTRY", "year"),
    model = "within",
    effect = "twoways",
    effect.id = TRUE,
    response_transformation = "log1p" 
  ),
  m_nb_fe = list(
    formula = Projects ~ log(Population) + factor(Rural) + Int_Access + log(GDP) + Education + log(Unemp) + (1 | COUNTRY) + (1 | COUNTRY:NUTS_3) + (1 | year),
    func = glmmTMB,
    family = nbinom2,
    response_transformation = "raw"
  ),
  gam_log1p = list(
    formula = log1p(Projects) ~ s(log(Population)) + factor(Rural) + s(Int_Access) + s(log(GDP)) + s(Education) + s(log(Unemp)) + factor(COUNTRY) + factor(year),
    func = gam,
    response_transformation = "log1p" 
  ),
  m_mult_mix_log1p = list(
    formula = log1p(Projects) ~ log(Population) + factor(Rural) + Int_Access + log(GDP) + Education + log(Unemp) + (1 | COUNTRY) + (1 | COUNTRY:NUTS_3) + (1 | year),
    func = lmer,
    response_transformation = "log1p"
  )
)

#------------------- 3. Cross-Validation Loop -------------------#
set.seed(17)
folds <- vfold_cv(df, v = 10)

results <- folds$splits %>%
  map_dfr(function(split) {
    train_data <- training(split)
    test_data <- testing(split)
    test_data <- pdata.frame(test_data, index = c("COUNTRY", "year"))
    
    imap_dfr(model_formulas, function(model_details, model_name) {
      args_list <- list(formula = model_details$formula, data = train_data)
      if ("index" %in% names(model_details)) {
        args_list$index <- model_details$index
      }
      model <- do.call(model_details$func, args_list)
      preds <- predict(model, newdata = test_data)
      
      if (model_details$response_transformation == "log") {
        preds <- exp(preds)
        true_values <- test_data$Projects.NA
      } else if (model_details$response_transformation == "log1p") {
        preds <- expm1(preds)
        true_values <- test_data$Projects
      } else if (model_details$response_transformation == "ihs") {
        preds <- sinh(preds)
        true_values <- test_data$Projects
      } else if (model_details$response_transformation == "raw") {
        true_values <- test_data$Projects
      }
      
      tibble(
        model = model_name,
        mae = mean(abs(true_values - preds), na.rm = TRUE),
        rho = cor(true_values, preds, use = "complete.obs")
      )
    })
  })

# make sure the models are in the correct order
desired_order <- names(model_formulas)
results$model <- factor(results$model, levels = desired_order)

#------------------- 4. Performance Evaluation and Plotting -------------------#
testPerformance <- results %>%
  pivot_longer(cols = c(mae, rho), names_to = "measure", values_to = "value")

cv.plot <- ggplot(testPerformance, aes(x = model, y = value, col = model, fill = model)) + 
  geom_point(position = position_jitter(0.1), shape = 21, col = "black", size = 3, stroke = 0.2) + 
  stat_summary(fun.y = mean, geom = "point", size = 2) + 
  stat_summary(fun.y = mean, geom = "point", size = 1, col = "white") + 
  facet_wrap(~ measure, scales = "free") +
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", fun.args = list(mult = 2), width = 0.2, lwd = 1) +
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", fun.args = list(mult = 2), width = 0.18, lwd = 0.1, col = "white") +
  scale_x_discrete(labels = c("OLS\n(ihs)", "OLS\n(log1p)", "FE\n(log1p)", "NB\nMixed", "GAM\n(log1p)", "Mult-L\nMixed log1p")) +
  labs(x = "Model", y = "Cross-validated MAE / Correlation") + 
  theme_bw() + theme(panel.grid = element_blank(), legend.position = "none", text = element_text(size = 14))

print(cv.plot)




######  GAM Smoothed curves partial effects Plots ######.


# Initialise list to store plots
draw_list <- list()

# Adjust y-axis limits to allow space for confidence intervals
y_limits_pop <- c(-1, 4.6)
y_limits_gdp <- c(-2, 1.5)
y_limits_unemp <- c(-0.5, 0.2)
y_limits_educ <- c(-0.5, 0.2)

# Custom labels for x-axis
custom_breaks = c(10000, 2000000, 4000000, 6000000)
custom_labels = c("1K", "2M", "4M", "6M")

# Set up common ggplot themes for all plots
common_theme <- list(
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    plot.title = element_text(size = 16, hjust = 0.5)
  ),
  geom_hline(yintercept = 0, color = "red", size = 0.3)
)


## For speeding up
# Reduce number of points for faster plotting
n_fewer_points = 50

# Use a subsample of the data for faster plotting (assuming your data frame is df_no_missing)
data_subsample = df_no_missing[sample(1:nrow(df_no_missing), 1000), ]

# Update the GAM model accordingly
m_gam_advanced_subsample <- update(m_gam, data = data_subsample)

# Modify draw() function to use fewer points
draw_list[['Population']] <- draw(m_gam_advanced_subsample, select = "s(Population)", n = n_fewer_points) + common_theme


# Effect of Population on Projects
draw_list[['Population']] <- draw(m_gam, select = "s(Population)") + common_theme +
  labs(title = NULL, y = "Projects") +
  scale_x_continuous(limits = c(1000, 6000000), breaks = custom_breaks, labels = custom_labels) +
  scale_y_continuous(limits = y_limits_pop)

# Effect of Internet Access on Projects
draw_list[['Int_Access']] <- draw(m_gam, select = "s(Int_Access)") + common_theme +
  labs(title = NULL, y = NULL)

# Effect of GDP on Projects
draw_list[['GDP']] <- draw(m_gam, select = "s(GDP)") + common_theme +
  labs(title = NULL, y = NULL) +
  scale_x_continuous(limits = c(0, 150000)) +
  scale_y_continuous(limits = y_limits_gdp)

# Effect of Education on Projects
draw_list[['Education']] <- draw(m_gam, select = "s(Education)") + common_theme +
  labs(title = NULL, y = "Projects") +
  scale_x_continuous(limits = c(11, 50)) + 
  scale_y_continuous(limits = y_limits_educ)

# Effect of Unemployment on Projects
draw_list[['Unemp']] <- draw(m_gam, select = "s(Unemp)") + common_theme +
  labs(title = NULL, y = NULL) +
  scale_x_continuous(limits = c(1, 30)) +
  scale_y_continuous( limits=y_limits_unemp)

# Combine all the plots into a single plot with 6 columns
final_plot <- wrap_plots(draw_list, ncol = 3, nrow = 2)


# Export the combined plots to a PNG
ggsave("combined_plots1111.png", final_plot, width = 18, height = 8, dpi = 300)
# Export the combined plots to a PDF
# ggsave("combined_plots.pdf", final_plot, width = 20, height = 4)


## Plot for Rural-Urban Divide factor term
# Convert the fitted GAM model to gamViz class
m_gam_viz <- getViz(m_gam)
# Plot
rural_plot <- plot(m_gam_viz, select = 6, allTerms = TRUE) +
  labs(title = NULL, y = NULL) +
  xlab("Region Type") + 
  ylab(NULL) +
  theme(
    axis.text = element_text(size = 13),
    plot.title = element_text(size = 16, hjust = 0.5),
    axis.title = element_text(size = 14, vjust = -2),
    panel.grid.major.y = element_line(color = "grey90", size = 0.4),  
    panel.grid.minor.y = element_line(color = "grey90", size = 0.4) 
  ) +
  scale_x_discrete(labels = c("1" = "Urban", "2" = "Intermediate", "3" = "Rural"))  # Custom x-axis labels


# Define the file name and dimensions
file_name <- "rural_plot1.png"
width_in_inches <- 6.3
height_in_inches <- 4
dpi <- 300

# Open a new PNG device
png(filename = file_name, width = width_in_inches * dpi, height = height_in_inches * dpi, res = dpi)

# Print the plot to the device
print(rural_plot)

# Close the device
dev.off()








#----------------------Models' Assumption Checks -----------------------------.



#######  OLS  #######.

## Linearity
### Fit the model
m_ols_log1p <- lm(Projects.log1p ~ log(Population) + factor(Rural) + Int_Access + log(GDP) + Education + log(Unemp),
                  data=df_no_missing %>% filter(all_vars==1))
### Predicted values
df_no_missing$predicted <- predict(m_ols_log1p, df_no_missing)
### Plot observed vs predicted values
ggplot(df_no_missing, aes(x = predicted, y = Projects.log1p)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "red")

## Independence
### Durbin-Watson test
library(car)
durbinWatsonTest(m_ols_log1p)

## Homoscedasticity
### Residuals vs Predicted values
ggplot(df_no_missing, aes(x = predicted, y = residuals(m_ols_log1p))) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "red")
### Breusch-Pagan test
bptest(m_ols_log1p)

## Normality
### QQ plot
qqnorm(residuals(m_ols_log1p))
qqline(residuals(m_ols_log1p))
### Shapiro-Wilk test
set.seed(123)  # For reproducibility
subsampled_residuals <- sample(residuals(m_ols_log1p), 5000)
shapiro.test(subsampled_residuals)

## No multicollinearity
### VIF
car::vif(m_ols_log1p)







####### Residuals Diagnostic for Mixed models #######.

# Extract residuals and standardized residuals from the model
residuals <- resid(m_mult_mix_log1p)
std_residuals <- residuals / sigma(m_mult_mix_log1p)

# Set up a 2x2 layout for plots
par(mfrow = c(2, 2), mar = c(4, 4, 2, 2))

# 1. Plot residuals vs fitted values
plot(fitted(m_mult_mix_log1p), residuals, 
     ylab = "Residuals", xlab = "Fitted Values", 
     main = "Residuals vs Fitted")
abline(h = 0, col = "red")

# 2. Plot normal Q-Q plot
qqnorm(std_residuals, main = "Normal Q-Q")
qqline(std_residuals, col = "red")

# 3. Plot scale-location (spread-location plot)
plot(fitted(m_mult_mix_log1p), sqrt(abs(std_residuals)), 
     ylab = "Square Root of Standardized Residuals", xlab = "Fitted Values",
     main = "Scale-Location")
abline(h = 0, col = "red")

# 4. Plot residuals by country with rotated labels
boxplot(residuals ~ df_no_missing$COUNTRY, ylab = "Residuals", xlab = "Country",
        main = "Residuals by Country", las=2)

# Add a main title for all plots
title("Residuals Diagnostic Plots", outer = TRUE, line = -21)






####### Statistical Tests: Multi-level vs Fixed-Effects #######.

# Run the fixed-effects model
m_fixed <- lm(log1p(Projects) ~ log(Population) + factor(Rural) + Int_Access + 
                log(GDP) + Education + log(Unemp) + factor(COUNTRY) + factor(year), 
              data = df_no_missing)


# 1. Likelihood Ratio Test
library(lmtest)
lrtest(m_fixed, m_mult_mix_log1p)

# 2. AIC and BIC
aic_fixed <- AIC(m_fixed)
bic_fixed <- BIC(m_fixed)

aic_mult <- AIC(m_mult_mix_log1p)
bic_mult <- BIC(m_mult_mix_log1p)
# Compare the AIC and BIC
print(paste("AIC Fixed: ", aic_fixed, " AIC Multi-level: ", aic_mult))
print(paste("BIC Fixed: ", bic_fixed, " BIC Multi-level: ", bic_mult))





###### Skewness Analysis by country ########.

# Group residuals by country
by_country <- tapply(residuals, df_no_missing$COUNTRY, function(x) x)

# Calculate and print skewness for each country
skewness_by_country <- sapply(by_country, skewness)
print(skewness_by_country)


### Residual Distribution Plots for the most skewed countries

# Create a dataframe with residuals and countries
df_res <- data.frame(residuals = residuals,
                     COUNTRY = df_no_missing$COUNTRY)

# Specify countries to inspect and initialize list for plots
countries_to_inspect <- c("LV", "EE", "LU")
plot_list <- list()

# Generate histogram plot for each specified country
for (country in countries_to_inspect) {
  plot_list[[country]] <- ggplot(df_res[df_res$COUNTRY == country,], aes(x=residuals)) +
    geom_histogram(fill = "blue", alpha = 0.5, bins = 30) +
    theme_minimal() +
    ggtitle(paste("Residual distribution for", country))
}

# Combine individual plots and display
combined_plot <- do.call(grid.arrange, c(plot_list, ncol = 3))
print(combined_plot)



#### inspection of Large Negative Residuals

# Add residuals to the df_no_missing dataframe
df_no_missing$residuals <- residuals

# Inspect and summarize observations with large negative residuals
df_no_missing %>%
  filter(residuals < quantile(residuals, 0.01)) %>%
  summary()




####### GAM Modle Fit and Residuals ########.


# Residuals vs Fitted: Check for non-linear patterns in residuals
plot(m_gam, residuals = TRUE, pch = 16, cex = 0.5)

# QQ Plot: Check if residuals are normally distributed
qqnorm(residuals(m_gam))
qqline(residuals(m_gam))

# gam.check: Provides several diagnostic plots and tests
gam.check(m_gam)

# Plot the Smoothed and Linear Terms: Visualise the effect of each predictor
plot(m_gam, pages=1, all.terms = TRUE)


## Assumption Checks

# Test for Autocorrelation: Check if residuals are independent
acf(residuals(m_gam))

# Check Concurvity: Check if smooth terms are explaining the same variance
concurvity(m_gam)

# Check Outliers: Identify any influential points
plot(resid(m_gam))
