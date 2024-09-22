library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(ggpubr)

#load data:
file <- file.choose(new = FALSE)
data <- read.csv(file,header=T,dec=".",sep=",")

data1<-subset(data, period="drought")
data2<-subset(data, period="rewetting")

# Reshape the data from wide to long format
long_data1 <- data1 %>%
  pivot_longer(cols = where(is.numeric) & matches("^(i|a|cy|X)"), names_to = "Category", values_to = "Value")

long_data1 <- long_data1 %>%
  mutate(Category = str_remove(Category, "^X"))


long_data2 <- data2 %>%
  pivot_longer(cols = where(is.numeric) & matches("^(i|a|cy|X)"), names_to = "Category", values_to = "Value")

long_data2 <- long_data2 %>%
  mutate(Category = str_remove(Category, "^X"))
# Calculate mean and standard error for each combination of Category, Group, and DeutLabel
summary_data1 <- long_data1 %>%
  group_by(Category, facet, Label) %>%
  summarize(
    mean_val = mean(Value),
    se = sd(Value) / sqrt(n()),
    .groups = "drop"
  )

summary_data2 <- long_data2 %>%
  group_by(Category, facet, Label) %>%
  summarize(
    mean_val = mean(Value),
    se = sd(Value) / sqrt(n()),
    .groups = "drop"
  )
# Generate the plot
p1 <- ggplot(summary_data1, aes(x = mean_val, y = Category, color = Label)) +
  geom_point(size = 2, position = position_dodge(width = 0.4)) +
  geom_errorbar(aes(xmin = mean_val - se, xmax = mean_val + se),
                width = 0.2, position = position_dodge(width = 0.4)) +
  scale_color_manual(values = c("black", "red")) +
  facet_wrap(~ facet) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10)
  ) +
  xlab("Atom %") +
  ylab("Fatty acid marker")

p1

p2 <- ggplot(summary_data2, aes(x = mean_val, y = Category, color = Label)) +
  geom_point(size = 2, position = position_dodge(width = 0.4)) +
  geom_errorbar(aes(xmin = mean_val - se, xmax = mean_val + se),
                width = 0.2, position = position_dodge(width = 0.4)) +
  scale_color_manual(values = c("black", "red")) +
  facet_wrap(~ facet) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10)
  ) +
  xlab("Atom %") +
  ylab("Fatty acid marker")

p2

#Plot all together
all<-ggarrange(
  p1, 
  p2,
  common.legend = TRUE, legend="bottom",
  ncol=2,
  nrow=1,
  widths = c(1,1),
  labels = c("a)", "b)")
)
all
ggsave(plot=all,"Fig. S14.png", width = 11, height = 5.5)


