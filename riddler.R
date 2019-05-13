# Riddler ----------------------------------------------------------------------
# Due 11:59 pm EST Sunday



# JOE DIMAGGIO -----------------------------------------------------------------
# Riddler Classic May 5, 2019 JOE DIMAGGIO
# https://fivethirtyeight.com/features/can-the-riddler-bros-beat-joe-dimaggios-
# hitting-streak/
# 
# From Steven Pratt, where have you gone, Joe DiMaggio? To Riddler Nation, 
# perhaps:
#   
#   Five brothers join the Riddler Baseball Independent Society, or RBIs. Each 
# of them enjoys a lengthy career of 20 seasons, with 160 games per season and 
# four plate appearances per game. (To make this simple, assume each plate 
# appearance results in a hit or an out, so there are no sac flies or walks to 
# complicate this math.)
# 
# Given that their batting averages are .200, .250, .300, .350 and .400, what 
# are each brother’s chances of beating DiMaggio’s 56-game hitting streak at 
# some point in his career? (Streaks can span across seasons.)
# 
# By the way, their cousin has a .500 average, but he will get tossed from the 
# league after his 10th season when he tests positive for performance enhancers. 
# What are his chances of beating the streak?


# batting average
avg <- .300

# number of years played
y <- 20

# number of games per year
g <- 160

# length of the streak
str <- 57

# plate appearances per game
pa_g <- 4

# total number of games played
tot_g <- g * y

# total number of plate appearances
tot_pa <- pa_g * g * y

# probability of not getting a hit in a game 
p_no_hit_g <- (1 - avg) ^ pa_g  

# probability of getting at least one hit, is 1 - p_not_hit_g
p_hit_g <- 1 - p_no_hit_g

# probabilities of getting a hit 57 times in a row
# p_str <- (p_hit_g) ^ 57
p_str <- (1 - ((1 - p_hit_g) ^ 57))

# number of opportunities to break the streak over the course of a career
# (can't break streak last 56 games)

op_str <- tot_g - 57 + 1

p_break <- op_str * p_str


# create a function to predict the probabilty 

streak <- function (avg, y = 20, g = 160, pa_g = 4, str = 57) {
  ((1 - ((1 - avg) ^ pa_g)) ^ str) * (y * g - str + 1)
}

streak(.500)
# not right


# End Riddler May 5, 2019 ------------------------------------------------------