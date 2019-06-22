
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
  setNames(data.frame(t(coef(the_fit))), c("intercept", "slope"))
}

setNames(data.frame(t(coef(lm(hwy ~ displ, mpg))), summary(the_fit)$p_value), c("intercept", "slope"))

mpg2 <- mpg %>%
  group_by(manufacturer, class) %>%
  do(lin_fit(.))


