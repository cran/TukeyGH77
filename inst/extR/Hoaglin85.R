# learn Hoaglin (1985)

library(TukeyGH77)
remotes::install_github('tingtingzhan/gg.tzh'); library(gg.tzh)

theme_bw() |> theme_set()

(gv = (1:5)/5)
(hv = (1:5)/10)

ggplot() + 
  paths_function(\(z,g) expm1(g*z)/g, dots = list(g = gv), hjust = .9, size = 2L) + 
  xlim(-2,2) + 
  labs(y = NULL, title = 'Fig 11-2')

ggplot() + 
  paths_function(\(z,h) dGH(z,h=h)-dnorm(z), dots = list(h = hv), n = 501L, hjust = .2, size = 2L) + 
  xlim(0,8) + 
  labs(y = NULL, title = 'Fig 11-8')

ggplot() + 
  paths_function(\(z,h) z*exp(h*z^2/2), dots = list(h = hv), hjust = .05, size = 2L) + 
  xlim(-2,2) + 
  labs(y = NULL, title = 'Fig 11-9')

ggplot() + 
  paths_function(\(z,h) z*exp(h*z^2/2), dots = list(h = -hv), hjust = .2, size = 2L) + 
  xlim(-6,6) + 
  labs(y = NULL, title = 'Not monotone')


