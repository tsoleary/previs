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

# probability of getting at least one hit in a game, is 1 - p_not_hit_g
p_hit_g <- 1 - p_no_hit_g

# probabilities of getting a hit 57 games in a row


# https://stats.stackexchange.com/questions/21825/probability-of-a-run-of-k-
# successes-in-a-sequence-of-n-bernoulli-trials/23762#23762

# Probability of a run of k successes in a sequence of n Bernoulli trials
# Ask Question
# 13

# I'm trying to find the probability of getting 8 trials in a row correct in 
# a block of 25 trials, you have 8 total blocks (of 25 trials) to get 8 trials 
# correct in a row. The probability of getting any trial correct based on 
# guessing is 1/3, after getting 8 in a row correct the blocks will end (so 
# getting more than 8 in a row correct is technically not possible). How would 
# I go about finding the probability of this occurring? I've been thinking along
# the lines of using (1/3)^8 as the probability of getting 8 in a row correct, 
# there are 17 possible chances to get 8 in a row in a block of 25 trials, if 
# I multiply 17 possibilities * 8 blocks I get 136, would 1-(1-(1/3)^8)^136 give 
# me the likelihood of getting 8 in a row correct in this situation or am I 
# missing something fundamental here?
#   

# I believe the problem with the argument given is that the events considered 
# are not independent. For example, consider a single block. If I tell you 
# that (a) there is no run of eight that starts at position 6, (b) there is a 
# run starting at position 7 and (c) there is no run starting at position 8, 
# what does that tell you about the probability of a run starting at positions, 
# say, 9 through 15? 


# r code to simulate this 
hits8 <- function() {
  x <- rbinom(26, 1, 1/3)                # 25 Binomial trials
  x[1] <- 0                              # ... and a 0 to get started with `diff`
  if(sum(x) >= 8) {                      # Are there at least 8 successes?
    max(diff(cumsum(x), lag=8)) >= 8     # Are there 8 successes in a row anywhere?
  } else {
    FALSE                                # Not enough successes for 8 in a row
  }
}
set.seed(17)                             # this is a way to set the random number
                                         # sequence to be the same for a trial
mean(replicate(10^5, hits8()))

# r code to simulate batting
hit <- function() {
  x <- rbinom(3200, 1, .300)            
  x[1] <- 0                             
  if(sum(x) >= 57) {                      
    # max(diff(cumsum(x), lag=57)) >= 57   https://www.rdocumentation.org/packages/base/versions/3.6.0/topics/diff  
  } else {
    FALSE                                
  }
}

mean(replicate(10^5, hit()))




# number of opportunities to break the streak over the course of a career
# (can't break streak last 56 games)

op_str <- tot_g - 57 + 1

p_break <- op_str * p_str


# create a function to predict the probabilty 

streak <- function (avg, y = 20, g = 160, pa_g = 4, str = 57) {
  ((1 - ((1 - avg) ^ pa_g)) ^ str) * (y * g - str + 1)
}

streak(.500)
# above 1
# not right


# End Riddler May 5, 2019 ------------------------------------------------------