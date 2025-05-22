

system.time(aldex_p1 <- aldex(table[, !colnames(table) %in% "taxonomy"], conditions = condiciones, mc.samples = 128, denom = "all"))


system.time({
  aldex_clr <- aldex.clr(table[, !colnames(table) %in% "taxonomy"],
                         condiciones,
                         mc.samples = 128,
                         denom = "all")
  
  effect_size <- aldex.effect(
    aldex_clr,
    verbose = TRUE,
    include.sample.summary = FALSE,
    useMC = TRUE,
    CI = FALSE)
  
  KW <- aldex.ttest(aldex_clr, verbose = FALSE)
  
  resultado <- cbind(effect_size, KW)
})
