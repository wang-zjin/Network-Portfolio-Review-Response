rm(list = ls())

setwd("~/Documents/GitHub/Network-Portfolio/Empirical_results_Involve_MV/Monthly_rebalance")

# Load Functions and other Files
source('./PackagesNetworkPortfolio.R')
source('./FunctionsNetworkPortfolio.R')

# load data
MV<-read_excel("Price_MV.xlsx",sheet="MV", na = "#N/A N/A")
# Identify columns with no NAs, excluding the 'Dates' column
non_na_cols <- sapply(MV[,-1], function(col) all(!is.na(col)))
# Include the 'Dates' column and columns with no NAs
MV_clean <- MV[, c(TRUE, non_na_cols)]
# Create the zoo object
ZOO <- zoo(MV_clean[,-1], order.by = as.Date(MV_clean$Dates, format = '%Y-%m-%d'))

# Prepare return data
MV <- ZOO
MV <- MV[-1, ]
MVstd<-xts(MV)
p=dim(MV)[2]
T.length=dim(MV)[1]

# Set labels
node.label=colnames(MV)
node.label<-gsub("Equity","",node.label)
node.label<-gsub("UN","",node.label)
node.label<-gsub("UW","",node.label)
node.label<-gsub(" ","",node.label)
names(MVstd) = node.label

# rolling window
W<-list()
for(t in 0: (floor((T.length-500)/22)-1)){
  W[[(t+1)]]=log(MVstd)[(1+t*22):(522+t*22),]
}
W_in<-list()
W_out<-list()
for(t in 0: (floor((T.length-500)/22)-1)){
  W_in[[(t+1)]]=W[[t+1]][c(1:500),]
  W_out[[(t+1)]]=W[[t+1]][c(501:522),]
}
T.windows<-length(W)

# Size list
Size_list<-list()
a=c()
for (t in 1: length(W_in)) {
  Size_list[[t]]=col_means(W_in[[t]])
  Size_list[[t]] <- (Size_list[[t]] - min(Size_list[[t]])) / (max(Size_list[[t]]) - min(Size_list[[t]]))
}

# Load centrality data
load("eigenvector_centrality_20241120.RData")
EC_in=eigenvector_centrality$eigenvector_absolute_value
EC_DS=eigenvector_centrality$eigenvector_centrality_Dantzig
ER_in=eigenvector_centrality$expected_return
C_in=eigenvector_centrality$correlation_matrix
COV_in=eigenvector_centrality$covariance_matrix

# Generate dates starting from MV_clean$Dates[500] with 22-day increments
start_date <- as.Date(MV_clean$Dates[500], format = '%Y-%m-%d')
date_sequence <- seq(from = start_date, by = 22, length.out = 88)

# Create data frame for size  and centrality
df_size <- data.frame(date_sequence, matrix(unlist(Size_list), nrow = length(Size_list), byrow = TRUE))
colnames(df_size) <- c("Date", names(Size_list[[1]]))

df_EC_DS <- data.frame(date_sequence, matrix(unlist(EC_DS), nrow = length(EC_DS), byrow = TRUE))
colnames(df_EC_DS) <- c("Date", names(EC_DS[[1]]))

# Reshape data frames to long format
size_long <- df_size %>%
  pivot_longer(
    cols = -Date,  # All columns except Date
    names_to = "Stock",  # Name of the new column for stock names
    values_to = "MV"  # Name of the new column for the values
  )

EC_DS_long <- df_EC_DS %>%
  pivot_longer(
    cols = -Date,  # All columns except Date
    names_to = "Stock",  # Name of the new column for stock names
    values_to = "EC_DS"  # Name of the new column for the values
  )

# Merge the two data frames by Date and Stock
merged_data <- merge(size_long, EC_DS_long, by = c("Date", "Stock"))

# Display the merged data
print("Merged Data:")
print(head(merged_data))

# Calculate correlation and perform a significance test
cor_test <- cor.test(merged_data$MV, merged_data$EC_DS)

# Print correlation result and p-value
print(paste("Correlation Coefficient:", cor_test$estimate))
print(paste("P-value:", cor_test$p.value))

# Check significance level
if (cor_test$p.value < 0.05) {
  print("The correlation is statistically significant (p < 0.05).")
} else {
  print("The correlation is not statistically significant (p >= 0.05).")
}


################################################################################
################## Calculate correlation for each date #########################
correlation_by_date <- merged_data %>%
  group_by(Date) %>%
  summarise(
    Correlation = cor(MV, EC_DS, use = "complete.obs"),
    P_value = cor.test(MV, EC_DS, use = "complete.obs")$p.value
  ) %>%
  mutate(Significant = ifelse(P_value < 0.05, "Significant", "Not Significant"))

# Display the correlation time series data
print("Correlation by Date:")
print(head(correlation_by_date))

# Plot the correlation time series without significance
plot_correlation_by_date <- ggplot(correlation_by_date, aes(x = Date, y = Correlation)) +
  geom_line(color = "blue", size = 1) +  # Time series line
  geom_point(aes(color = "black"), size = 2, alpha = 0.8) +  # Highlight data points
  scale_color_manual(values = c("Significant" = "green", "Not Significant" = "red")) +
  labs(
    x = "",
    y = "Correlation"
  ) +
  scale_x_date(
    breaks = seq(as.Date("2016-01-01"), as.Date("2022-01-01"), by="1 year"),
    date_labels = "%Y"
  ) +
  theme(
    panel.grid = element_blank(), # Remove all grid lines
    panel.background = element_rect(fill = "transparent", color = NA),  # Make panel background transparent
    plot.background = element_rect(fill = "transparent", color = NA),   # Make plot background transparent
    panel.border = element_rect(color="black", fill=NA, size=1),
    axis.title.y = element_text(size = 16)
  )
# Save the plot
ggsave(filename = "correlation_time_series_no_significance.png",
       plot = plot_correlation_by_date,
       width = 10,
       height = 6,
       dpi = 300)
