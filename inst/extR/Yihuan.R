
library(TukeyGH77)
remotes::install_github('tingtingzhan/gg.tzh'); library(gg.tzh)

ggplot() + paths_function(dGH, dots = list(
  g = c(0, 0, 0, .1, .2), 
  h = c(0, .1, .4, 0, .2)
), size = 2, hjust = .1) + 
  xlim(-4, 4) + labs(y = NULL, title = 'Yihuan\'s paper')
  
