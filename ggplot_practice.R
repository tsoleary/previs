
# standard scatter plot
ggplot(data = mpg) + 
  geom_point(mapping = aes(x = displ, y = hwy, color = class)) 

ggplot(data = mpg) + 
  geom_point(mapping = aes(x = displ, y = hwy, color = manufacturer)) 

# jitter the points to allow for the points that are overlapping to be seen
ggplot(data = mpg) + 
  geom_point(mapping = aes(x = displ, y = hwy, color = class),
             position = "jitter") 

# same as above
ggplot(data = mpg) + 
  geom_jitter(mapping = aes(x = displ, y = hwy, color = class)) 

# make points slightly transparent with alpha
ggplot(data = mpg) + 
  geom_jitter(mapping = aes(x = displ, y = hwy, color = class), alpha = 0.5)

# themes
ggplot(data = mpg) + 
  geom_jitter(mapping = aes(x = displ, y = hwy, color = class), 
              alpha = 0.5) +
  theme_classic()

# size and outline of points
ggplot(data = mpg) + 
  geom_jitter(mapping = aes(x = displ, y = hwy, fill = class), 
              alpha = 0.5, size = 3, pch = 21,  color = "black") +
  theme_classic()

# the one above looks pretty good. I like it


# multiple plots per square
ggplot(data = mpg) + 
  geom_point(mapping = aes(x = displ, y = hwy)) +
  facet_wrap(~ class, nrow = 2)


# linear regression stuff ------------------------------------------------------

# how to get a different linear regression for each class
lin_fit <- function(dat){
  the_fit <- lm(dat$hwy ~ dat$displ, dat)
  p_val <- anova(the_fit)$'Pr(>F)'[1]
  slo_int <- data.frame(t(coef(the_fit)))
  r_sq <- summary(the_fit)$r.squared
  result <- cbind(slo_int, r_sq, p_val)
  colnames(result) <- c("intercept", "slope", "r_squared", "p_value")
  return(result)
}

mpg2 <- mpg %>%
  group_by(class) %>%
  do(lin_fit(.))


mpg3 <- mpg %>%
  group_by(manufacturer, class) %>%
  do(lin_fit(.))


mpg %>% 
  nest(-drv) %>% 
  mutate(model = map(data, ~ lm(hwy~displ, data = .x)),
         adj.r.squared = map_dbl(model, ~ signif(summary(.x)$adj.r.squared, 5)),
         intercept = map_dbl(model, ~ signif(.x$coef[[1]],5)),
         slope = map_dbl(model, ~ signif(.x$coef[[2]], 5)),
         pvalue = map_dbl(model, ~ signif(summary(.x)$coef[2,4], 5)) 
  ) %>% 
  select(-data, -model) %>% 
  left_join(mpg) %>% 
  
  ggplot(aes(displ, hwy)) +
  geom_point() +
  geom_smooth(se = FALSE, method = "lm") +
  facet_wrap(~drv) +
  geom_text(aes(3, 40, label = paste("Adj R2 = ", adj.r.squared, "\n",
                                     "Intercept =",intercept, "\n",
                                     "Slope =", slope, "\n",
                                     "P =", pvalue)))

