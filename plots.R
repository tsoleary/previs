# Plotting tests ---------------------------------------------------------------

library(ggplot2)

library(plotly)

g <- ggplot(txhousing, aes(x = date, y = sales, group = city)) +
  geom_line(alpha = 0.4)

ggplotly(g, tooltip = c("city"))


