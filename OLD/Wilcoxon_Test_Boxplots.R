Alpha.Wilcoxon <- Data.filtered %>%
  # Group and filter the data with variance == 0
  group_by(Range_Intervals,Metric) %>%
  filter(Var_Value != 0) %>%
  # Group the data
  group_by(Taxa,Metric) %>%
  # Compute the kruskall-test
  wilcox_test(formula = Value ~ Stripe, p.adjust.method = "bonferroni") %>%
  # Transform "group1" and "group2" into numeric
  mutate_at(c("group1", "group2"), as.numeric ) %>%
  # Add a column that is the "stripe distance" between the distribution compared
  mutate(Stripe_Distance = abs(group1 - group2))

# We need a column "y.position" for the plotting of the significance brackets
y.position <- Total.filtered %>%
  # Group the data
  group_by(Taxa,Metric) %>%
  # Add the y.position with an increase of X%. 
  summarise(y.position = max(Value) * 1.05) 
  
# Filter the precedent data_frame to only keep the values of stripe distance == 1
Metric.Wilcoxon.Adj <- Metric.Wilcoxon %>%
  # Filter the data to only keep adjacent stripes eventually 
  # filter(Stripe_Distance == 1) %>%
  # Add the values of y.position 
  left_join(y.position)

# We will use ddboxplot instead of geom boxplot because it is easier to plot the significance brackets on each facet. 

Alpha.Boxplot <- Total.filtered %>%
  # Draw the plot
  ggboxplot(
    x = "Break", y = "Value", fill = "Break",
    facet = c("Metric", "Taxa"),
    scales = "free",
    ggtheme = arrange_theme()
    ) + 
  # Add the significance levels
  # stat_pvalue_manual(Metric.Wilcoxon.Adj, hide.ns = TRUE, step.increase = 0.1, step.group.by = "Metric") +
    # Labels
  xlab("Altitudinal stripes") +
  ylab("Metric values") +
  labs(
    color = "Stripe altitudinal limits",
    title = paste0("BoxPlot of metric values ~ Altitudinal stripes number"),
    subtitle = paste0("Wilcox-tests were realized between adjacent stripes and significant results are displayed.",
                          "\n*: p <= 0.05 / **: p <= 0.01 / ***: p <= 0.001 / ****: p <= 0.0001")
  ) + # Theme
  theme(axis.text.x=element_text(angle = 90, hjust = 0))

