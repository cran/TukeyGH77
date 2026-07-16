


(1:1e3L) |>
  vapply(FUN = \(i) {
    set.seed(seed = i)
    m = rGH(n = 1e3L, g = -.3, h = .1) |> 
      letterValue_gh()
    # b = x |> letterValue_gh(g_select = 'h.optim')
    # if ((a['h'] == 0) && (b['h'] > 0)) print(i) # salvageable
    # if ((a['h'] == 0) && (b['h'] == 0)) print(i) # hopeless
    ret <- (m['h'] == 0) # hopeless
    on.exit(return(ret))
    # if (ret) print(i) # hopeless
  }, FUN.VALUE = NA) |>
  mean.default()
  