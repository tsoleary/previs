# Plotting tests ---------------------------------------------------------------

library(ggplot2)

library(plotly)

g <- ggplot(txhousing, aes(x = date, y = sales, group = city)) +
  geom_line(alpha = 0.4)

ggplotly(g, tooltip = c("city"))



# From internet
# https://r4ds.had.co.nz/data-visualisation.html
 
library(tidyverse)

mpg

ggplot(data = mpg) +
  geom_point(mapping = aes(x = displ, y = hwy))

