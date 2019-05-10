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

# batting average of the five brothers

ba <- .300

# ba_b1 <- .200
# ba_b2 <- .250
# ba_b3 <- .300
# ba_b4 <- .350
# ba_b5 <- .400
# ba_cuz <- .500

# number of years played

y <- 20

# y_b <- 20
# y_cuz <- 10

# number of games per year

g <- 160

# let's look at the totals for fun
#ignoring the cousing for now

pa_g <- 4
tot_g <- g * y
tot_pa <- pa_g * g * y
tot_hit <- tot_pa * ba

# let's find the probability of not getting a hit in a game 

p_not_hit_g <- (1 - ba) ^ pa_g  



# so if you are batting .300 then the probabiliy of not getting a hit in a game
# when you have 4 PA is .2401. 

# so getting at least one hit, NOT(p_not_hit_g)

p__hit_g <- 1 - p_not_hit_g

# then a probability of two games in row
# I assume that it just multiplies but it may be more complicated




# let's define breaking as at least 57 games with hits in a row


# because you cannot break the streak if it starts during the 
# last 56 games of your carreer. Or does that not matter?

tot_g - 56 * 2 






# End Riddler May 5, 2019 ------------------------------------------------------