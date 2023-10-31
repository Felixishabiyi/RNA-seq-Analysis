#Loading the data
var <- read.csv("C:/Users/Administrator/Documents/R/Dataset/variants_long_table.csv")

#Printing the head
head(var)

#Checking the dimension
dim(var)

#The structure
str(var)

#summary of a POS column
summary(var$POS)

#Viewing in an interactive format
View(var)

#Loading packages
library('dplyr')
library('ggplot2')

#checking column names
colnames(var)

# Select columns 1, 4, and 5 with default display
select(var, SAMPLE, REF, ALT) %>% head(3)

# Select all columns except the column “CALLER” with selected display
select(var, -CALLER) %>% head(3)

# Transform the data frame into a tibble
var_tb <- as_tibble(var)

select(var_tb, SAMPLE, REF, ALT) %>% head(3)

# Select rows with selected display using base R code
var_tb[var_tb$SAMPLE == "SRR13500958",]

# Select rows with selected display using dplyr functions
filter(var_tb, SAMPLE == "SRR13500958") %>% head(3)

# Select sample type (rows) and variables (columns) with selected display
var_tb %>% filter(SAMPLE == "SRR13500958") %>% 
  select(CHROM, POS, REF, ALT) %>% head(3)

# To select all data related to the sample specified
var_tb %>% filter(SAMPLE == "SRR13500958") %>% 
  select(CHROM, POS, REF, ALT, DP)

# To select only values for which DP>=500 for the same sample
var_tb %>% filter(SAMPLE == "SRR13500958" & DP>=500) %>% 
  select(CHROM, POS, REF, ALT, DP)

# To select only values for which DP>=1000 for the same sample
var_tb %>% filter(SAMPLE == "SRR13500958" & DP>=1000) %>% 
  select(CHROM, POS, REF, ALT, DP)

# Count how many rows are associated with each sample in the data
var_tb %>% count(SAMPLE)

#Sorting the count
var_tb %>% count(SAMPLE, sort = TRUE)

# Distribution of genes per sample and counts
var_tb %>% count(SAMPLE, GENE, sort = TRUE) %>% head()

# Maximum value of column DP
max(var_tb$DP)

# Minimum value of column DP
min(var_tb$DP)

# Mean value of column DP
mean(var_tb$DP)

# Compute a LOG2 transformation on the DP values
var_tb_log <- var_tb %>% mutate(DP_log2 = log2(DP))

head(var_tb_log)

# View a selected content including the new column
select(var_tb_log, SAMPLE, REF, ALT, DP, DP_log2) %>% head()

# Show the maximum value of DP for each sample
var_tb %>% group_by(SAMPLE) %>% summarize(max(DP))

# Show the minimum value of DP for each sample
var_tb %>% group_by(SAMPLE) %>% summarize(min(DP))

var_tb %>% group_by(SAMPLE) %>% summarize(mean(DP))

var_tb %>% group_by(SAMPLE) %>% count(CHROM)

var_tb %>% mutate(DP100 = DP/100) %>% select(SAMPLE, DP, DP100) %>% str()

#Points
ggplot(data = var_tb, aes(x=SAMPLE, y=DP)) + geom_point() + ylim(0,10000)

#Boxplots
ggplot(data = var_tb, aes(x=SAMPLE, y=DP)) + geom_boxplot() + ylim(0,10000)

#Points
ggplot(data = var_tb, aes(x=SAMPLE, y=DP)) + geom_point() + 
  scale_y_continuous(name="dp", limits=c(0, 10000))

#Boxplot
ggplot(data = var_tb, aes(x=SAMPLE, y=DP)) + 
  geom_boxplot() + scale_y_continuous(name="dp", limits=c(0, 10000))

ggplot(data = var_tb, aes(x=SAMPLE, y=DP)) + 
  geom_point() + scale_y_continuous(trans='log10')

ggplot(data = var_tb, aes(x=SAMPLE, y=DP)) + geom_point() + scale_y_log10()

ggplot(data = var_tb, aes(x=SAMPLE, y=DP)) + geom_boxplot() + scale_y_log10()

# Colours of shapes
ggplot(data = var_tb, aes(x=SAMPLE, y=DP, colour = SAMPLE)) + 
  geom_boxplot() + ylim(0,10000)

# Colours for filling
ggplot(data = var_tb, aes(x=SAMPLE, y=DP, fill = SAMPLE)) + 
  geom_boxplot() + ylim(0,10000)

library(RColorBrewer)

ggplot(data = var_tb, aes(x=SAMPLE, y=DP, fill= SAMPLE)) + 
  geom_boxplot() + ylim(0,10000) + scale_fill_brewer(palette="RdYlBu")


ggplot(data = var_tb, aes(x=SAMPLE, y=DP, fill= SAMPLE)) + 
  geom_boxplot() + ylim(0,10000) + scale_fill_brewer(palette="RdYlBu") + 
  theme(legend.position="top")

ggplot(data = var_tb, aes(x=SAMPLE, y=DP, fill= SAMPLE)) + 
  geom_boxplot() + ylim(0,10000) + scale_fill_brewer(palette="RdYlBu") + 
  theme(legend.position="none")

ggplot(data = var_tb, aes(x=SAMPLE, y=DP, fill= SAMPLE)) + 
  geom_boxplot() + ylim(0,10000) + scale_fill_brewer(palette="RdYlBu") + 
  theme(legend.position="bottom") + 
  labs(title="DP_per_Sample", x="SampleID", y = "DP")

ggplot(data = var_tb, aes(x=SAMPLE, y=DP)) +  
  geom_point(shape = 21, fill = "#e4dbc1", color = "#b92e17", size = 6) + 
  ylim(0,10000)

ggplot(data = var_tb, aes(x=CHROM, y=DP, fill= SAMPLE)) + 
  geom_boxplot() + ylim(0,10000) + scale_fill_brewer(palette="RdYlBu") + 
  labs(title="DP_per_Chromosome") + facet_grid(. ~ SAMPLE)

# Define a variable with plotting options
p_DP_CHROM <- ggplot(data = var_tb, aes(x=CHROM, y=DP, fill= SAMPLE)) + 
  ylim(0,10000) + scale_fill_brewer(palette="RdYlBu") + 
  labs(title="DP_per_Chromosome") + theme(legend.position="bottom")

# Test boxplots with faceting
p_DP_CHROM + geom_boxplot() + facet_grid(. ~ SAMPLE)

# Combine violin plots and boxplots with faceting
p_DP_CHROM + geom_violin(trim=FALSE) + 
  facet_grid(. ~ SAMPLE) + 
  geom_boxplot(width=0.1)

#Count number of different effects per sample
p_EFFECT <- ggplot(data = var_tb, aes(x=EFFECT, fill= SAMPLE)) + 
  scale_fill_brewer(palette="RdBu") + labs(title="Effect_per_Sample") + 
  theme(legend.position="bottom")

p_EFFECT + geom_bar()

# Flip orientation
p_EFFECT_flip <- ggplot(data = var_tb, aes(y=EFFECT, fill= SAMPLE)) + 
  scale_fill_brewer(palette="RdBu") + labs(title="Effect_per_Sample") + 
  theme(legend.position="bottom")

p_EFFECT_flip + geom_bar()

# Count the number of different effects
var_tb %>% count(EFFECT)

# Count the number of different effects and link them to sample information
var_tb %>% count(EFFECT, SAMPLE, sort = TRUE)

# Counting the effects per gene
var_tb %>% count(EFFECT, GENE, sort = TRUE)

# Filtering option 1 to select for effect on stop
filter(var_tb, EFFECT == "stop_lost" | EFFECT == "stop_gained")

# Filtering on effect and selected columns
filter(var_tb, EFFECT %in% c("stop_lost", "stop_gained")) %>% 
  select(SAMPLE, CHROM, GENE, EFFECT)

# Define your variable
p_DP_POS <- ggplot(data = var_tb, aes(x=POS, y=DP, fill= SAMPLE)) + 
  scale_fill_brewer(palette="RdBu") + labs(title="DP_per_Position") + 
  theme(legend.position="bottom")

# Plot
p_DP_POS + geom_point(shape = 21, size = 5)

# Plot with transparency options
p_DP_POS + geom_point(shape = 21, size = 5, alpha = 0.7)

p_DP_POS <- ggplot(data = var_tb, aes(x=DP, y=ALT, fill= SAMPLE)) + 
  scale_fill_brewer(palette="RdBu") + labs(title="DP_VS_ALT") + 
  theme(legend.position="bottom")

p_DP_POS + geom_point(shape = 21, size = 5)





