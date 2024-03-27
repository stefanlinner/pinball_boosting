# Sorucing scripts -------------------------------------------------------------

source("libraries.R", encoding = "UTF-8")
source("simulation_setups.R", encoding = "UTF-8")

# Loading data -----------------------------------------------------------------

sim_res_al1brq <- readRDS("simulation_output/simulation_results_rev_al1brq.RDS")
sim_res_l2brq <- readRDS("simulation_output/simulation_results_rev_l2brq.RDS")

quantile_risk <- function(y, f, tau){
  loss <- (tau * (y - f) * ((y - f) > 0) +
             (tau - 1) * (y - f) * ((y - f) <= 0))
  mean(loss)
}

unfill_vec <- function(x) {
  same <- x == dplyr::lag(x)
  ifelse(!is.na(same) & same, NA, x)
}

revision <- TRUE

options(knitr.kable.NA = '')

# Simulation results -----------------------------------------------------------

## Quantile risk ---------------------------------------------------------------

empirical_risk <- 
  function(.contaminated){
    map_dfr(
      unique(sim_res_al1brq$parameter$setup), 
      function(.setup){
        error_out <- map_dfr(
          unique(sim_res_al1brq$parameter$error),
          function(.error){
            method_out <- map_dfr(
              c("QR-QR.Boost", "QR-LS.Boost"),
              function(.method){
                tau_out <- 
                  map_dfc(
                    c(0.1, 0.3, 0.5, 0.7, 0.9),
                    function(.tau){
                      init_set <- 0
                      data.validation <- sim[[.setup]]$data_gen(102, n = 1000, error = .error, tau = .tau, contaminated = .contaminated, revision = revision)
                      y.validation <- data.validation[["y"]]
                      X.validation <- data.validation[, -1] %>% as.matrix() %>% cbind(rep(1, 1000), .)
                      
                      
                      if(.method == "QR-QR.Boost"){
                        
                        results <- 
                          sim_res_al1brq$output %>% 
                          filter(
                            setup == .setup, 
                            error == .error, 
                            tau == .tau, 
                            contaminated == .contaminated, 
                            init == init_set
                          )
                        
                        predictions <- 
                          map(
                            results$value,
                            function(.x){
                              beta <- .x[["coefs"]] %>% map_dbl(~.x)
                              y.prediction <- X.validation %*% beta
                              y.prediction
                            }
                          )
                        
                        risk <- 
                          map_dbl(predictions, function(.preds){
                            quantile_risk(y = y.validation, f = .preds, tau = .tau)
                          }) %>% 
                          mean() %>% 
                          round(3)
                        
                        out <- 
                          list(risk)
                        
                        names(out) <- as.character(.tau)
                        
                        return(out) 
                      }
                      
                      if(.method == "QR-LS.Boost"){
                        
                        results <- 
                          sim_res_l2brq$output %>% 
                          filter(
                            setup == .setup, 
                            error == .error, 
                            tau == .tau, 
                            contaminated == .contaminated, 
                            init == init_set
                          )
                        
                        predictions <- 
                          map(
                            results$value,
                            function(.x){
                              beta <- .x[["coefs"]] %>% map_dbl(~.x)
                              y.prediction <- X.validation %*% beta
                              y.prediction
                            }
                          )
                        
                        risk <- 
                          map_dbl(predictions, function(.preds){
                            quantile_risk(y = y.validation, f = .preds, tau = .tau)
                          }) %>% 
                          mean() %>% 
                          round(3)
                        
                        out <- 
                          list(risk)
                        
                        names(out) <- as.character(.tau)
                        
                        return(out) 
                      }
                      
                    }
                  )
                
                tibble(
                  Method = .method, 
                  tau_out
                )
              }
            )
            tibble(
              `Error distribution` = .error, 
              method_out
            )
          }
        )
        tibble(
          `Parameter setup` = .setup, 
          error_out
        )
      }
    )
  }

risk_not_contaminated <- empirical_risk(.contaminated = FALSE)

risk_not_contaminated %>% 
  group_by(`Parameter setup`, `Error distribution`) %>% 
  mutate(
    across(
      .cols = c(`0.1`, `0.3`, `0.5`, `0.7`, `0.9`),
      .fns = function(.x){
        if_else(.x == min(.x), str_c("\\blue{", .x, "}"), as.character(.x))
      }
    )
  ) %>% 
  ungroup() %>% 
  kbl("latex", escape = FALSE, linesep = "")

risk_contaminated <- empirical_risk(.contaminated = TRUE)

risk_contaminated %>% 
  group_by(`Parameter setup`, `Error distribution`) %>% 
  mutate(
    across(
      .cols = c(`0.1`, `0.3`, `0.5`, `0.7`, `0.9`),
      .fns = function(.x){
        if_else(.x == min(.x), str_c("\\blue{", .x, "}"), as.character(.x))
      }
    )
  ) %>% 
  ungroup() %>% 
  kbl("latex", escape = FALSE, linesep = "")

## MSE -------------------------------------------------------------------------


mse <- 
  function(.contaminated){
    map_dfr(
      unique(sim_res_al1brq$parameter$setup), 
      function(.setup){
        error_out <- map_dfr(
          unique(sim_res_al1brq$parameter$error) %>% setdiff("mixed"),
          function(.error){
            
            if(.setup %in% c("hom", "het")){
              betas <- map_chr(
                0:1, 
                function(.x){
                  paste0("hatbeta_tau",.x)
                }
              )
            } else if(.setup == "highdim"){
              betas <- map_chr(
                0:4, 
                function(.x){
                  paste0("hatbeta_tau",.x)
                }
              )
            } else {
              betas <- map_chr(
                0:6, 
                function(.x){
                  paste0("hatbeta_tau",.x)
                }
              )
            }
            
            beta_out <- 
              map_dfr(
                betas,
                function(.beta){
                  
                  method_out <- map_dfr(
                    c("QR-QR.Boost", "QR-LS.Boost"),
                    function(.method){
                      tau_out <- 
                        map_dfc(
                          c(0.1, 0.3, 0.5, 0.7, 0.9),
                          function(.tau){
                            init_set <- 0
                            
                            if(.setup == "multi2"){
                              
                              tau.vec <- c(0.1, 0.3, 0.5, 0.7, 0.9)
                              
                              betas <- sim[[.setup]]$beta[which(.tau == tau.vec),]
                              
                              if(.error == "norm"){
                                
                                betatau <- map_dbl(seq_along(betas), function(.x){
                                  betas[.x] + sim[[.setup]]$alpha[.x] * qnorm(p = .tau)
                                })
                                
                              } else if(.error == "tdist"){
                                
                                betatau <- map_dbl(seq_along(betas), function(.x){
                                  betas[.x] + sim[[.setup]]$alpha[.x] * qt(p = .tau, df = 2)
                                })
                                
                              } else if(.error == "gamma"){
                                
                                betatau <- map_dbl(seq_along(betas), function(.x){
                                  betas[.x] + sim[[.setup]]$alpha[.x] * qgamma(p = .tau, shape = 2, scale = 1)
                                })
                                
                              } else if(.error == "tdist_1"){
                                
                                betatau <- map_dbl(seq_along(betas), function(.x){
                                  betas[.x] + sim[[.setup]]$alpha[.x] * qt(p = .tau, df = 1)
                                })
                                
                              }
                              
                            } else {
                              if(.error == "norm"){
                                
                                betatau <- map_dbl(seq_along(sim[[.setup]]$beta), function(.x){
                                  sim[[.setup]]$beta[.x] + sim[[.setup]]$alpha[.x] * qnorm(p = .tau)
                                })
                                
                              } else if(.error == "tdist"){
                                
                                betatau <- map_dbl(seq_along(sim[[.setup]]$beta), function(.x){
                                  sim[[.setup]]$beta[.x] + sim[[.setup]]$alpha[.x] * qt(p = .tau, df = 2)
                                })
                                
                              } else if(.error == "gamma"){
                                
                                betatau <- map_dbl(seq_along(sim[[.setup]]$beta), function(.x){
                                  sim[[.setup]]$beta[.x] + sim[[.setup]]$alpha[.x] * qgamma(p = .tau, shape = 2, scale = 1)
                                })
                                
                              } else if(.error == "tdist_1"){
                                
                                betatau <- map_dbl(seq_along(sim[[.setup]]$beta), function(.x){
                                  sim[[.setup]]$beta[.x] + sim[[.setup]]$alpha[.x] * qt(p = .tau, df = 1)
                                })
                                
                              }
                            }
                            
                            if(.setup == "highdim"){
                              betatau <- betatau[1:5]
                            }
                            
                            beta_ind <- .beta %>% stringr::str_extract("[:digit:]") %>% as.numeric()
                            
                            if(.method == "QR-QR.Boost"){
                              
                              
                              results <- 
                                sim_res_al1brq$output %>% 
                                filter(
                                  setup == .setup, 
                                  error == .error, 
                                  tau == .tau, 
                                  contaminated == .contaminated, 
                                  init == init_set
                                )
                              
                              estimates <- 
                                map_dbl(
                                  results$value, 
                                  function(.res){
                                    .res[["coefs"]][[beta_ind + 1]]
                                  }
                                )
                              
                              mse <- 
                                mean((estimates - betatau[beta_ind + 1])^2) %>% 
                                round(3)
                              
                              out <- 
                                list(mse)
                              
                              names(out) <- as.character(.tau)
                              
                              return(out) 
                            }
                            
                            if(.method == "QR-LS.Boost"){
                              
                              results <- 
                                sim_res_l2brq$output %>% 
                                filter(
                                  setup == .setup, 
                                  error == .error, 
                                  tau == .tau, 
                                  contaminated == .contaminated, 
                                  init == init_set
                                )
                              
                              estimates <- 
                                map_dbl(
                                  results$value, 
                                  function(.res){
                                    .res[["coefs"]][[beta_ind + 1]]
                                  }
                                )
                              
                              mse <- 
                                mean((estimates - betatau[beta_ind + 1])^2) %>% 
                                round(3)
                              
                              out <- 
                                list(mse)
                              
                              names(out) <- as.character(.tau)
                              
                              return(out)
                            }
                            
                          }
                        )
                      tibble(
                        Method = .method, 
                        tau_out
                      )
                    }
                  )
                  tibble(
                    MSE = .beta, 
                    method_out
                  )
                }
              )
            tibble(
              `Error distribution` = .error, 
              beta_out
            )
          }
        )
        tibble(
          `Parameter setup` = .setup, 
          error_out
        )
      }
    )
    
  }

mse_not_contaminated <- mse(.contaminated = FALSE)

mse_not_contaminated %>% 
  group_by(`Parameter setup`, `Error distribution`, MSE) %>% 
  mutate(
    across(
      .cols = c(`0.1`, `0.3`, `0.5`, `0.7`, `0.9`),
      .fns = function(.x){
        if_else(.x == min(.x), str_c("\\blue{", .x, "}"), as.character(.x))
      }
    )
  ) %>% 
  ungroup() %>% 
  mutate(`Parameter setup` = as.character(`Parameter setup`)) %>% 
  mutate(
    across(
      .cols = c(`Parameter setup`, `Error distribution`, MSE),
      .fns = function(.x){
        unfill_vec(.x)
      }
    )
  ) %>% 
  kbl("latex", escape = FALSE, linesep = "")

mse_contaminated <- mse(.contaminated = TRUE)

mse_contaminated %>% 
  group_by(`Parameter setup`, `Error distribution`, MSE) %>% 
  mutate(
    across(
      .cols = c(`0.1`, `0.3`, `0.5`, `0.7`, `0.9`),
      .fns = function(.x){
        if_else(.x == min(.x), str_c("\\blue{", .x, "}"), as.character(.x))
      }
    )
  ) %>% 
  ungroup() %>% 
  mutate(`Parameter setup` = as.character(`Parameter setup`)) %>% 
  mutate(
    across(
      .cols = c(`Parameter setup`, `Error distribution`, MSE),
      .fns = function(.x){
        unfill_vec(.x)
      }
    )
  ) %>% 
  kbl("latex", escape = FALSE, linesep = "")

## Bias ------------------------------------------------------------------------

bias <- 
  function(.contaminated){
    map_dfr(
      unique(sim_res_al1brq$parameter$setup), 
      function(.setup){
        error_out <- map_dfr(
          unique(sim_res_al1brq$parameter$error) %>% setdiff("mixed"),
          function(.error){
            
            if(.setup %in% c("hom", "het")){
              betas <- map_chr(
                0:1, 
                function(.x){
                  paste0("hatbeta_tau",.x)
                }
              )
            } else if(.setup == "highdim"){
              betas <- map_chr(
                0:4, 
                function(.x){
                  paste0("hatbeta_tau",.x)
                }
              )
            } else {
              betas <- map_chr(
                0:6, 
                function(.x){
                  paste0("hatbeta_tau",.x)
                }
              )
            }
            
            beta_out <- 
              map_dfr(
                betas,
                function(.beta){
                  method_out <- map_dfr(
                    c("QR-QR.Boost", "QR-LS.Boost"),
                    function(.method){
                      tau_out <- 
                        map_dfc(
                          c(0.1, 0.3, 0.5, 0.7, 0.9),
                          function(.tau){
                            init_set <- 0
                            
                            if(.setup == "multi2"){
                              
                              tau.vec <- c(0.1, 0.3, 0.5, 0.7, 0.9)
                              
                              betas <- sim[[.setup]]$beta[which(.tau == tau.vec),]
                              
                              if(.error == "norm"){
                                
                                betatau <- map_dbl(seq_along(betas), function(.x){
                                  betas[.x] + sim[[.setup]]$alpha[.x] * qnorm(p = .tau)
                                })
                                
                              } else if(.error == "tdist"){
                                
                                betatau <- map_dbl(seq_along(betas), function(.x){
                                  betas[.x] + sim[[.setup]]$alpha[.x] * qt(p = .tau, df = 2)
                                })
                                
                              } else if(.error == "gamma"){
                                
                                betatau <- map_dbl(seq_along(betas), function(.x){
                                  betas[.x] + sim[[.setup]]$alpha[.x] * qgamma(p = .tau, shape = 2, scale = 1)
                                })
                                
                              } else if(.error == "tdist_1"){
                                
                                betatau <- map_dbl(seq_along(betas), function(.x){
                                  betas[.x] + sim[[.setup]]$alpha[.x] * qt(p = .tau, df = 1)
                                })
                                
                              }
                              
                            } else {
                              if(.error == "norm"){
                                
                                betatau <- map_dbl(seq_along(sim[[.setup]]$beta), function(.x){
                                  sim[[.setup]]$beta[.x] + sim[[.setup]]$alpha[.x] * qnorm(p = .tau)
                                })
                                
                              } else if(.error == "tdist"){
                                
                                betatau <- map_dbl(seq_along(sim[[.setup]]$beta), function(.x){
                                  sim[[.setup]]$beta[.x] + sim[[.setup]]$alpha[.x] * qt(p = .tau, df = 2)
                                })
                                
                              } else if(.error == "gamma"){
                                
                                betatau <- map_dbl(seq_along(sim[[.setup]]$beta), function(.x){
                                  sim[[.setup]]$beta[.x] + sim[[.setup]]$alpha[.x] * qgamma(p = .tau, shape = 2, scale = 1)
                                })
                                
                              } else if(.error == "tdist_1"){
                                
                                betatau <- map_dbl(seq_along(sim[[.setup]]$beta), function(.x){
                                  sim[[.setup]]$beta[.x] + sim[[.setup]]$alpha[.x] * qt(p = .tau, df = 1)
                                })
                                
                              }
                            }
                            
                            if(.setup == "highdim"){
                              betatau <- betatau[1:5]
                            }
                            
                            beta_ind <- .beta %>% stringr::str_extract("[:digit:]") %>% as.numeric()
                            
                            if(.method == "QR-QR.Boost"){
                              
                              
                              results <- 
                                sim_res_al1brq$output %>% 
                                filter(
                                  setup == .setup, 
                                  error == .error, 
                                  tau == .tau, 
                                  contaminated == .contaminated, 
                                  init == init_set
                                )
                              
                              estimates <- 
                                map_dbl(
                                  results$value, 
                                  function(.res){
                                    .res[["coefs"]][[beta_ind + 1]]
                                  }
                                )
                              
                              bias <- 
                                mean(estimates - betatau[beta_ind + 1]) %>% 
                                round(3)
                              
                              out <- 
                                list(bias)
                              
                              names(out) <- as.character(.tau)
                              
                              return(out) 
                            }
                            
                            if(.method == "QR-LS.Boost"){
                              
                              results <- 
                                sim_res_l2brq$output %>% 
                                filter(
                                  setup == .setup, 
                                  error == .error, 
                                  tau == .tau, 
                                  contaminated == .contaminated, 
                                  init == init_set
                                )
                              
                              estimates <- 
                                map_dbl(
                                  results$value, 
                                  function(.res){
                                    .res[["coefs"]][[beta_ind + 1]]
                                  }
                                )
                              
                              bias <- 
                                mean(estimates - betatau[beta_ind + 1]) %>% 
                                round(3)
                              
                              out <- 
                                list(bias)
                              
                              names(out) <- as.character(.tau)
                              
                              return(out)
                            }
                            
                          }
                        )
                      tibble(
                        Method = .method, 
                        tau_out
                      )
                    }
                  )
                  tibble(
                    Bias = .beta, 
                    method_out
                  )
                }
              )
            tibble(
              `Error distribution` = .error, 
              beta_out
            )
          }
        )
        tibble(
          `Parameter setup` = .setup, 
          error_out
        )
      }
    )
    
  }

bias_not_contaminated <- bias(.contaminated = FALSE)

bias_not_contaminated %>% 
  filter(`Error distribution` != "tdist_1") %>% 
  group_by(`Parameter setup`, `Error distribution`, Bias) %>% 
  mutate(
    across(
      .cols = c(`0.1`, `0.3`, `0.5`, `0.7`, `0.9`),
      .fns = function(.x){
        if_else(abs(.x) == min(abs(.x)), str_c("\\blue{", .x, "}"), as.character(.x))
      }
    )
  ) %>% 
  ungroup() %>% 
  mutate(`Parameter setup` = as.character(`Parameter setup`)) %>% 
  mutate(
    across(
      .cols = c(`Parameter setup`, `Error distribution`, Bias),
      .fns = function(.x){
        unfill_vec(.x)
      }
    )
  ) %>% 
  kbl("latex", escape = FALSE, linesep = "")

bias_contaminated <- bias(.contaminated = TRUE)

bias_contaminated %>% 
  filter(`Error distribution` != "tdist_1") %>% 
  group_by(`Parameter setup`, `Error distribution`, Bias) %>% 
  mutate(
    across(
      .cols = c(`0.1`, `0.3`, `0.5`, `0.7`, `0.9`),
      .fns = function(.x){
        if_else(abs(.x) == min(abs(.x)), str_c("\\blue{", .x, "}"), as.character(.x))
      }
    )
  ) %>% 
  ungroup() %>% 
  mutate(`Parameter setup` = as.character(`Parameter setup`)) %>% 
  mutate(
    across(
      .cols = c(`Parameter setup`, `Error distribution`, Bias),
      .fns = function(.x){
        unfill_vec(.x)
      }
    )
  ) %>% 
  kbl("latex", escape = FALSE, linesep = "")

# Number of iterations ---------------------------------------------------------

iterations <- 
  function(.contaminated){
    map_dfr(
      unique(sim_res_al1brq$parameter$setup), 
      function(.setup){
        
        error_out <- map_dfr(
          unique(sim_res_al1brq$parameter$error),
          function(.error){
            
            method_out <- map_dfr(
              c("QR-QR.Boost", "QR-LS.Boost"),
              function(.method){
                
                tau_out <- map_dfc(
                  c(0.1, 0.3, 0.5, 0.7, 0.9),
                  function(.tau){
                    
                    init_set <- 0
                    
                    if(.method == "QR-QR.Boost"){
                      
                      zbrq <- sim[[.setup]]$zbrq
                      
                      results <- 
                        sim_res_al1brq$output %>% 
                        filter(
                          setup == .setup, 
                          error == .error, 
                          tau == .tau, 
                          contaminated == .contaminated, 
                          init == init_set
                        )
                      
                      median.mstop <-
                        map_dbl(
                          results$value, 
                          function(.res){
                            mstop <- .res[["emp.risks"]] %>% 
                              which.min() %>% 
                              zbrq[.]
                            mstop
                          }
                        ) %>% 
                        median()
                      
                      out <- 
                        list(median.mstop)
                      
                      names(out) <- as.character(.tau)
                      
                      return(out) 
                    }
                    
                    if(.method == "QR-LS.Boost"){
                      
                      zqb <- sim[[.setup]]$zqb
                      
                      results <- 
                        sim_res_l2brq$output %>% 
                        filter(
                          setup == .setup, 
                          error == .error, 
                          tau == .tau, 
                          contaminated == .contaminated, 
                          init == init_set
                        )
                      
                      median.mstop <-
                        map_dbl(
                          results$value, 
                          function(.res){
                            mstop <- .res[["emp.risks"]] %>% 
                              which.min() %>% 
                              zqb[.]
                            mstop
                          }
                        ) %>% 
                        median()
                      
                      out <- 
                        list(median.mstop)
                      
                      names(out) <- as.character(.tau)
                      
                      return(out) 
                    }
                    
                  }
                )
                
                tibble(
                  Method = .method,
                  tau_out
                )
                
              }
            )
            
            tibble(
              `Error distr.` = .error,
              method_out
            )
            
          }
        )
        tibble(
          `Parameter setup` = .setup,
          error_out
        )
        
      }
      
    )
    
  }

iterations_not_contaminated <- iterations(.contaminated = FALSE)

iterations_not_contaminated %>% 
  filter(`Error distr.` != "tdist_1") %>% 
  group_by(`Parameter setup`, `Error distr.`) %>% 
  mutate(
    across(
      .cols = c(`0.1`, `0.3`, `0.5`, `0.7`, `0.9`),
      .fns = function(.x){
        if_else(.x == min(.x), str_c("\\blue{", .x, "}"), as.character(.x))
      }
    )
  ) %>% 
  ungroup() %>% 
  mutate(
    `Parameter setup` = as.character(`Parameter setup`),
    `Error distr.` = as.character(`Error distr.`)
  ) %>% 
  mutate(
    across(
      .cols = c(`Parameter setup`, `Error distr.`),
      .fns = function(.x){
        unfill_vec(.x)
      }
    )
  ) %>% 
  kbl("latex", escape = FALSE, linesep = "")

iterations_contaminated <- iterations(.contaminated = TRUE)

iterations_contaminated %>% 
  filter(`Error distr.` != "tdist_1") %>% 
  group_by(`Parameter setup`, `Error distr.`) %>% 
  mutate(
    across(
      .cols = c(`0.1`, `0.3`, `0.5`, `0.7`, `0.9`),
      .fns = function(.x){
        if_else(.x == min(.x), str_c("\\blue{", .x, "}"), as.character(.x))
      }
    )
  ) %>% 
  ungroup() %>% 
  mutate(
    `Parameter setup` = as.character(`Parameter setup`),
    `Error distr.` = as.character(`Error distr.`)
  ) %>% 
  mutate(
    across(
      .cols = c(`Parameter setup`, `Error distr.`),
      .fns = function(.x){
        unfill_vec(.x)
      }
    )
  ) %>% 
  kbl("latex", escape = FALSE, linesep = "")

## Sensitivity -----------------------------------------------------------------

sensitivity_highdim <- 
  function(.contaminated){
    map_dfr(
      c("norm", "tdist", "gamma", "mixed"),
      function(.error){
        
        method_out <- map_dfr(
          c("QR-QR.Boost", "QR-LS.Boost"),
          function(.method){
            
            tau_out <- map_dfc(
              c(0.1,0.3,0.5,0.7,0.9),
              function(.tau){
                
                init_set <- 0
                
                if(.method == "QR-QR.Boost"){
                  
                  results <- 
                    sim_res_al1brq$output %>% 
                    filter(
                      setup == "highdim", 
                      error == .error, 
                      tau == .tau, 
                      contaminated == .contaminated, 
                      init == init_set
                    )
                  
                  true.positive.rate <-
                    map_dbl(
                      results$value, 
                      function(.res){
                        mean(.res[["mean.appearance"]][1:4] != 0)
                      }
                    ) %>% 
                    mean() %>% round(3)
                  
                  out <- 
                    list(true.positive.rate)
                  
                  names(out) <- as.character(.tau)
                  
                  return(out) 
                  
                }
                
                
                if(.method == "QR-LS.Boost"){
                  
                  results <- 
                    sim_res_l2brq$output %>% 
                    filter(
                      setup == "highdim", 
                      error == .error, 
                      tau == .tau, 
                      contaminated == .contaminated, 
                      init == init_set
                    )
                  
                  true.positive.rate <-
                    map_dbl(
                      results$value, 
                      function(.res){
                        mean(.res[["mean.appearance"]][1:4] != 0)
                      }
                    ) %>% 
                    mean() %>% round(3)
                  
                  out <- 
                    list(true.positive.rate)
                  
                  names(out) <- as.character(.tau)
                  
                  return(out) 
                  
                }
                
              }
            )
            
            tibble(
              Method = .method, 
              tau_out
            )
            
          }
        )
        
        tibble(
          `Error distribution` = .error, 
          method_out
        )
        
      }
    )
    
  }

sensitivity_highdim_not_contaminated <- sensitivity_highdim(.contaminated = FALSE)

sensitivity_highdim_not_contaminated %>% 
  filter(`Error distribution` != "tdist_1") %>% 
  group_by(`Error distribution`) %>% 
  mutate(
    across(
      .cols = c(`0.1`, `0.3`, `0.5`, `0.7`, `0.9`),
      .fns = function(.x){
        if_else(.x == max(.x), str_c("\\blue{", .x, "}"), as.character(.x))
      }
    )
  ) %>% 
  ungroup() %>% 
  mutate(
    `Error distribution` = as.character(`Error distribution`)
  ) %>% 
  mutate(
    across(
      .cols = c(`Error distribution`),
      .fns = function(.x){
        unfill_vec(.x)
      }
    )
  ) %>% 
  kbl("latex", escape = FALSE, linesep = "")

sensitivity_highdim_contaminated <- sensitivity_highdim(.contaminated = TRUE)

sensitivity_highdim_contaminated %>% 
  filter(`Error distribution` != "tdist_1") %>% 
  group_by(`Error distribution`) %>% 
  mutate(
    across(
      .cols = c(`0.1`, `0.3`, `0.5`, `0.7`, `0.9`),
      .fns = function(.x){
        if_else(.x == max(.x), str_c("\\blue{", .x, "}"), as.character(.x))
      }
    )
  ) %>% 
  ungroup() %>% 
  mutate(
    `Error distribution` = as.character(`Error distribution`)
  ) %>% 
  mutate(
    across(
      .cols = c(`Error distribution`),
      .fns = function(.x){
        unfill_vec(.x)
      }
    )
  ) %>% 
  kbl("latex", escape = FALSE, linesep = "")

## Specificity -----------------------------------------------------------------

specificity_highdim <- 
  function(.contaminated){
    map_dfr(
      c("norm", "tdist", "gamma", "mixed"),
      function(.error){
        
        method_out <- map_dfr(
          c("QR-QR.Boost", "QR-LS.Boost"),
          function(.method){
            
            tau_out <- map_dfc(
              c(0.1,0.3,0.5,0.7,0.9),
              function(.tau){
                
                init_set <- 0
                
                if(.method == "QR-QR.Boost"){
                  
                  results <- 
                    sim_res_al1brq$output %>% 
                    filter(
                      setup == "highdim", 
                      error == .error, 
                      tau == .tau, 
                      contaminated == .contaminated, 
                      init == init_set
                    )
                  
                  spec <-
                    map_dbl(
                      results$value, 
                      function(.res){
                        mean(.res[["mean.appearance"]][5:99] == 0)
                      }
                    ) %>% 
                    mean() %>% round(3)
                  
                  out <- 
                    list(spec)
                  
                  names(out) <- as.character(.tau)
                  
                  return(out) 
                  
                }
                
                
                if(.method == "QR-LS.Boost"){
                  
                  results <- 
                    sim_res_l2brq$output %>% 
                    filter(
                      setup == "highdim", 
                      error == .error, 
                      tau == .tau, 
                      contaminated == .contaminated, 
                      init == init_set
                    )
                  
                  spec <-
                    map_dbl(
                      results$value, 
                      function(.res){
                        mean(.res[["mean.appearance"]][5:99] == 0)
                      }
                    ) %>% 
                    mean() %>% round(3)
                  
                  out <- 
                    list(spec)
                  
                  names(out) <- as.character(.tau)
                  
                  return(out) 
                  
                }
                
              }
            )
            
            tibble(
              Method = .method, 
              tau_out
            )
            
          }
        )
        
        tibble(
          `Error distribution` = .error, 
          method_out
        )
        
      }
    )
    
  }


specificity_highdim_not_contaminated <- specificity_highdim(.contaminated = FALSE)

specificity_highdim_not_contaminated %>% 
  filter(`Error distribution` != "tdist_1") %>% 
  group_by(`Error distribution`) %>% 
  mutate(
    across(
      .cols = c(`0.1`, `0.3`, `0.5`, `0.7`, `0.9`),
      .fns = function(.x){
        if_else(.x == max(.x), str_c("\\blue{", .x, "}"), as.character(.x))
      }
    )
  ) %>% 
  ungroup() %>% 
  mutate(
    `Error distribution` = as.character(`Error distribution`)
  ) %>% 
  mutate(
    across(
      .cols = c(`Error distribution`),
      .fns = function(.x){
        unfill_vec(.x)
      }
    )
  ) %>% 
  kbl("latex", escape = FALSE, linesep = "")

specificity_highdim_contaminated <- specificity_highdim(.contaminated = TRUE)

specificity_highdim_contaminated %>% 
  filter(`Error distribution` != "tdist_1") %>% 
  group_by(`Error distribution`) %>% 
  mutate(
    across(
      .cols = c(`0.1`, `0.3`, `0.5`, `0.7`, `0.9`),
      .fns = function(.x){
        if_else(.x == max(.x), str_c("\\blue{", .x, "}"), as.character(.x))
      }
    )
  ) %>% 
  ungroup() %>% 
  mutate(
    `Error distribution` = as.character(`Error distribution`)
  ) %>% 
  mutate(
    across(
      .cols = c(`Error distribution`),
      .fns = function(.x){
        unfill_vec(.x)
      }
    )
  ) %>% 
  kbl("latex", escape = FALSE, linesep = "")

## Tau-fit ---------------------------------------------------------------------

sim_sample_data_boosting <- readRDS("simulation_output/boosting_sample_data.RDS")

tau_insample <- 
  function(.contaminated){
    map_dfr(
      unique(sim_res_al1brq$parameter$setup), 
      function(.setup){
        error_out <- map_dfr(
          unique(sim_res_al1brq$parameter$error),
          function(.error){
            
            method_out <- map_dfr(
              c("QR-QR.Boost", "QR-LS.Boost"),
              function(.method){
                tau_out <- 
                  map_dfc(
                    c(0.1, 0.3, 0.5, 0.7, 0.9),
                    function(.tau){
                      init_set <- 0
                      
                      if(.method == "QR-QR.Boost"){
                        
                        sample_data <- 
                          sim_sample_data_boosting$output %>% 
                          filter(
                            setup == .setup, 
                            error == .error, 
                            tau == .tau, 
                            contaminated == .contaminated, 
                            init == init_set
                          )
                        
                        results <- 
                          sim_res_al1brq$output %>% 
                          filter(
                            setup == .setup, 
                            error == .error, 
                            tau == .tau, 
                            contaminated == .contaminated, 
                            init == init_set
                          )
                        
                        tau_fit <- 
                          map_dbl(
                            seq_along(results$value),
                            function(.x){
                              
                              X <- cbind(
                                1,
                                sample_data$value[[.x]][, -1]
                              ) %>% 
                                as.matrix()
                              
                              beta <- results$value[[.x]][["coefs"]] %>% map_dbl(~.x)
                              fitted <- (X %*% beta) %>% as.vector()
                              
                              y_quantile <- quantile(sample_data$value[[.x]][["y"]], probs = .tau)
                              
                              nom <- 
                                quantile_risk(y = sample_data$value[[.x]][["y"]], f = fitted, tau = .tau)
                              
                              denom <- 
                                quantile_risk(y = sample_data$value[[.x]][["y"]], f = y_quantile, tau = .tau)
                              
                              1 - (nom/denom)
                              
                            }
                          ) %>% 
                          mean() %>% 
                          round(3)
                        
                        out <- 
                          list(tau_fit)
                        
                        names(out) <- as.character(.tau)
                        
                        return(out) 
                      }
                      
                      if(.method == "QR-LS.Boost"){
                        
                        sample_data <- 
                          sim_sample_data_boosting$output %>% 
                          filter(
                            setup == .setup, 
                            error == .error, 
                            tau == .tau, 
                            contaminated == .contaminated, 
                            init == init_set
                          )
                        
                        results <- 
                          sim_res_l2brq$output %>% 
                          filter(
                            setup == .setup, 
                            error == .error, 
                            tau == .tau, 
                            contaminated == .contaminated, 
                            init == init_set
                          )
                        
                        tau_fit <- 
                          map_dbl(
                            seq_along(results$value),
                            function(.x){
                              
                              X <- cbind(
                                1,
                                sample_data$value[[.x]][, -1]
                              ) %>% 
                                as.matrix()
                              
                              beta <- results$value[[.x]][["coefs"]] %>% map_dbl(~.x)
                              fitted <- (X %*% beta) %>% as.vector()
                              
                              y_quantile <- quantile(sample_data$value[[.x]][["y"]], probs = .tau)
                              
                              nom <- 
                                quantile_risk(y = sample_data$value[[.x]][["y"]], f = fitted, tau = .tau)
                              
                              denom <- 
                                quantile_risk(y = sample_data$value[[.x]][["y"]], f = y_quantile, tau = .tau)
                              
                              1 - (nom/denom)
                              
                            }
                          ) %>% 
                          mean() %>% 
                          round(3)
                        
                        out <- 
                          list(tau_fit)
                        
                        names(out) <- as.character(.tau)
                        
                        return(out) 
                      }
                      
                    }
                  )
                
                tibble(
                  Method = .method, 
                  tau_out
                )
              }
            )
            tibble(
              `Error distribution` = .error, 
              method_out
            )
          }
        )
        tibble(
          `Parameter setup` = .setup, 
          error_out
        )
      }
    )
    
  }

tau_insample_not_contaminated <- tau_insample(.contaminated = FALSE)

tau_insample_not_contaminated %>% 
  filter(`Error distribution` != "tdist_1") %>% 
  group_by(`Parameter setup`, `Error distribution`) %>% 
  mutate(
    across(
      .cols = c(`0.1`, `0.3`, `0.5`, `0.7`, `0.9`),
      .fns = function(.x){
        if_else(abs(.x) == min(abs(.x)), str_c("\\blue{", .x, "}"), as.character(.x))
      }
    )
  ) %>% 
  ungroup() %>% 
  mutate(
    `Parameter setup` = as.character(`Parameter setup`),
    `Error distribution` = as.character(`Error distribution`)
  ) %>% 
  mutate(
    across(
      .cols = c(`Parameter setup`, `Error distribution`),
      .fns = function(.x){
        unfill_vec(.x)
      }
    )
  ) %>% 
  kbl("latex", escape = FALSE, linesep = "")

tau_insample_contaminated <- tau_insample(.contaminated = TRUE)

tau_insample_contaminated %>% 
  filter(`Error distribution` != "tdist_1") %>% 
  group_by(`Parameter setup`, `Error distribution`) %>% 
  mutate(
    across(
      .cols = c(`0.1`, `0.3`, `0.5`, `0.7`, `0.9`),
      .fns = function(.x){
        if_else(abs(.x) == min(abs(.x)), str_c("\\blue{", .x, "}"), as.character(.x))
      }
    )
  ) %>% 
  ungroup() %>% 
  mutate(
    `Parameter setup` = as.character(`Parameter setup`),
    `Error distribution` = as.character(`Error distribution`)
  ) %>% 
  mutate(
    across(
      .cols = c(`Parameter setup`, `Error distribution`),
      .fns = function(.x){
        unfill_vec(.x)
      }
    )
  ) %>% 
  kbl("latex", escape = FALSE, linesep = "")
