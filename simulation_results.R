# Sorucing scripts -------------------------------------------------------------

source("libraries.R", encoding = "UTF-8")
source("simulation_setups.R", encoding = "UTF-8")

revision <- TRUE

# Loading data -----------------------------------------------------------------

if(revision){
  sim_res_al1brq <- readRDS("simulation_output/simulation_results_rev_al1brq.RDS")
  sim_res_l2brq <- readRDS("simulation_output/simulation_results_rev_l2brq.RDS")
} else {
  sim_res_al1brq <- readRDS("simulation_output/simulation_results_al1brq.RDS")
  sim_res_l2brq <- readRDS("simulation_output/simulation_results_l2brq.RDS")
}

sim_res_rq <- readRDS("simulation_output/simulation_results_rq.RDS")

quantile_risk <- function(y, f, tau){
  loss <- (tau * (y - f) * ((y - f) > 0) +
             (tau - 1) * (y - f) * ((y - f) <= 0))
  mean(loss)
}

unfill_vec <- function(x) {
  same <- x == dplyr::lag(x)
  ifelse(!is.na(same) & same, NA, x)
}

options(knitr.kable.NA = '')


# Simulation results -----------------------------------------------------------

## Quantile risk ---------------------------------------------------------------

empirical_risk <- 
  function(.contaminated){
    map_dfr(
      c("hom", "het", "multi", "multi2", "highdim"), 
      function(.setup){
        error_out <- map_dfr(
          unique(sim_res_al1brq$parameter$error),
          function(.error){
            
            if(.setup %in% c("multi", "multi2")){
              methods <- c("QR-QR.Boost", "QR-LS.Boost", "RQ", "RQAic")
            } else {
              methods <- c("QR-QR.Boost", "QR-LS.Boost", "RQ")
            }
            
            method_out <- map_dfr(
              methods,
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
                      
                      if(.method == "RQ"){
                        
                        results <- 
                          sim_res_rq$output %>% 
                          filter(
                            setup == .setup, 
                            error == .error, 
                            tau == .tau, 
                            contaminated == .contaminated
                          )
                        
                        predictions <- 
                          map(
                            results$value,
                            function(.x){
                              beta <- .x[["coefs_rq"]] %>% as.vector(mode = "numeric")
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
                      
                      if(.method == "RQAic"){
                        
                        results <- 
                          sim_res_rq$output %>% 
                          filter(
                            setup == .setup, 
                            error == .error, 
                            tau == .tau, 
                            contaminated == .contaminated
                          )
                        
                        predictions <- 
                          map(
                            results$value,
                            function(.x){
                              beta <- .x[["coefs_rq_bs"]] %>% as.vector(mode = "numeric")
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
  filter(`Error distribution` != "tdist_1") %>% 
  group_by(`Parameter setup`, `Error distribution`) %>% 
  mutate(
    across(
      .cols = c(`0.1`, `0.3`, `0.5`, `0.7`, `0.9`),
      .fns = function(.x){
        if_else(.x == min(.x), str_c("\\blue{", format(.x, digits = 3, nsmall = 3), "}"), as.character(format(.x, digits = 3, nsmall = 3)))
      }
    )
  ) %>% 
  ungroup() %>% 
  mutate(
    across(where(is.factor), as.character)
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

risk_contaminated <- empirical_risk(.contaminated = TRUE)

risk_contaminated %>% 
  filter(`Error distribution` != "tdist_1") %>% 
  group_by(`Parameter setup`, `Error distribution`) %>% 
  mutate(
    across(
      .cols = c(`0.1`, `0.3`, `0.5`, `0.7`, `0.9`),
      .fns = function(.x){
        if_else(.x == min(.x), str_c("\\blue{", format(.x, digits = 3, nsmall = 3), "}"), as.character(format(.x, digits = 3, nsmall = 3)))
      }
    )
  ) %>% 
  ungroup() %>% 
  mutate(
    across(where(is.factor), as.character)
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

## MSE -------------------------------------------------------------------------
# \widehat{\beta}_{\tau 2}

mse <- 
  function(.contaminated){
    map_dfr(
      c("hom", "het", "multi", "multi2", "highdim"), 
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
                  
                  
                  if(.setup %in% c("multi", "multi2")){
                    methods <- c("QR-QR.Boost", "QR-LS.Boost", "RQ", "RQAic")
                  } else {
                    methods <- c("QR-QR.Boost", "QR-LS.Boost", "RQ")
                  }
                  
                  method_out <- map_dfr(
                    methods,
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
                            
                            if(.method == "RQ"){
                              
                              results <- 
                                sim_res_rq$output %>% 
                                filter(
                                  setup == .setup, 
                                  error == .error, 
                                  tau == .tau, 
                                  contaminated == .contaminated
                                )
                              
                              estimates <- 
                                map_dbl(
                                  results$value, 
                                  function(.res){
                                    .res[["coefs_rq"]][beta_ind + 1] %>% as.numeric()
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
                            
                            if(.method == "RQAic"){
                              
                              results <- 
                                sim_res_rq$output %>% 
                                filter(
                                  setup == .setup, 
                                  error == .error, 
                                  tau == .tau, 
                                  contaminated == .contaminated
                                )
                              
                              estimates <- 
                                map_dbl(
                                  results$value, 
                                  function(.res){
                                    .res[["coefs_rq_bs"]][beta_ind + 1] %>% as.numeric()
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
  filter(`Error distribution` != "tdist_1") %>% 
  group_by(`Parameter setup`, `Error distribution`, MSE) %>% 
  mutate(
    across(
      .cols = c(`0.1`, `0.3`, `0.5`, `0.7`, `0.9`),
      .fns = function(.x){
        if_else(.x == min(.x), str_c("\\blue{", format(.x, digits = 3, nsmall = 3), "}"), as.character(format(.x, digits = 3, nsmall = 3)))
      }
    )
  ) %>% 
  ungroup() %>% 
  mutate(
    across(where(is.factor), as.character)
  ) %>% 
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
  filter(`Error distribution` != "tdist_1") %>% 
  group_by(`Parameter setup`, `Error distribution`, MSE) %>% 
  mutate(
    across(
      .cols = c(`0.1`, `0.3`, `0.5`, `0.7`, `0.9`),
      .fns = function(.x){
        if_else(.x == min(.x), str_c("\\blue{", format(.x, digits = 3, nsmall = 3), "}"), as.character(format(.x, digits = 3, nsmall = 3)))
      }
    )
  ) %>% 
  ungroup() %>% 
  mutate(
    across(where(is.factor), as.character)
  ) %>% 
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
      c("hom", "het", "multi", "multi2", "highdim"), 
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
                  
                  
                  if(.setup %in% c("multi", "multi2")){
                    methods <- c("QR-QR.Boost", "QR-LS.Boost", "RQ", "RQAic")
                  } else {
                    methods <- c("QR-QR.Boost", "QR-LS.Boost", "RQ")
                  }
                  
                  method_out <- map_dfr(
                    methods,
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
                            
                            if(.method == "RQ"){
                              
                              results <- 
                                sim_res_rq$output %>% 
                                filter(
                                  setup == .setup, 
                                  error == .error, 
                                  tau == .tau, 
                                  contaminated == .contaminated
                                )
                              
                              estimates <- 
                                map_dbl(
                                  results$value, 
                                  function(.res){
                                    .res[["coefs_rq"]][beta_ind + 1] %>% as.numeric()
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
                            
                            if(.method == "RQAic"){
                              
                              results <- 
                                sim_res_rq$output %>% 
                                filter(
                                  setup == .setup, 
                                  error == .error, 
                                  tau == .tau, 
                                  contaminated == .contaminated
                                )
                              
                              estimates <- 
                                map_dbl(
                                  results$value, 
                                  function(.res){
                                    .res[["coefs_rq_bs"]][beta_ind + 1] %>% as.numeric()
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
        if_else(abs(.x) == min(abs(.x)), str_c("\\blue{", format(.x, digits = 3, nsmall = 3), "}"), as.character(format(.x, digits = 3, nsmall = 3)))
      }
    )
  ) %>% 
  ungroup() %>% 
  mutate(
    across(where(is.factor), as.character)
  ) %>% 
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
        if_else(abs(.x) == min(abs(.x)), str_c("\\blue{", format(.x, digits = 3, nsmall = 3), "}"), as.character(format(.x, digits = 3, nsmall = 3)))
      }
    )
  ) %>% 
  ungroup() %>% 
  mutate(
    across(where(is.factor), as.character)
  ) %>% 
  mutate(
    across(
      .cols = c(`Parameter setup`, `Error distribution`, Bias),
      .fns = function(.x){
        unfill_vec(.x)
      }
    )
  ) %>% 
  kbl("latex", escape = FALSE, linesep = "")

## MPI -------------------------------------------------------------------------

mpi <- 
  function(.contaminated){
    map_dfr(
      c("multi", "multi2"), 
      function(.setup){
        tau_out <- 
          map_dfr(
            c(0.1, 0.3, 0.5, 0.7, 0.9),
            function(.tau){
              
              error_out <- 
                map_dfr(
                  unique(sim_res_al1brq$parameter$error), 
                  function(.error){
                    
                    method_out <- 
                      map_dfr(
                        c("QR-QR.Boost", "QR-LS.Boost"), 
                        function(.method){
                          
                          betas <- map_chr(
                            1:6, 
                            function(.x){
                              paste0("hatbeta_tau",.x)
                            }
                          )
                          
                          beta_out <- 
                            map_dfc(
                              betas,
                              function(.beta){
                                
                                init_set <- 0
                                
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
                                  
                                  mpi <- 
                                    map_dbl(
                                      results$value, 
                                      function(.res){
                                        .res[["mean.appearance"]][[paste0("beta_",beta_ind)]]
                                      }
                                    ) %>% 
                                    mean() %>% 
                                    round(3)
                                  
                                  out <- 
                                    list(mpi)
                                  
                                  names(out) <- .beta
                                  
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
                                  
                                  mpi <- 
                                    map_dbl(
                                      results$value, 
                                      function(.res){
                                        .res[["mean.appearance"]][[paste0("beta_",beta_ind)]]
                                      }
                                    ) %>% 
                                    mean() %>% 
                                    round(3)
                                  
                                  out <- 
                                    list(mpi)
                                  
                                  names(out) <- .beta
                                  
                                  return(out)
                                }
                                
                              }
                            )
                          
                          tibble(
                            Method = .method,
                            beta_out
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
                tau = .tau,
                error_out
              )
            }
          )
        tibble(
          `Parameter setup` = .setup,
          tau_out
        )
      }
    )
    
  }

mpi_not_contaminated <- mpi(.contaminated = FALSE)

mpi_not_contaminated %>% 
  filter(!`Error distribution` %in% c("tdist_1", "mixed")) %>% 
  mutate(
    across(where(is.factor), as.character)
  ) %>% 
  mutate(
    across(
      .cols = c(`Parameter setup`, tau, `Error distribution`),
      .fns = function(.x){
        unfill_vec(.x)
      }
    )
  ) %>% 
  kbl("latex", escape = FALSE, linesep = "")

mpi_contaminated <- mpi(.contaminated = TRUE)

mpi_contaminated %>% 
  filter(!`Error distribution` %in% c("tdist_1", "mixed")) %>% 
  mutate(
    across(where(is.factor), as.character)
  ) %>% 
  mutate(
    across(
      .cols = c(`Parameter setup`, tau, `Error distribution`),
      .fns = function(.x){
        unfill_vec(.x)
      }
    )
  ) %>% 
  kbl("latex", escape = FALSE, linesep = "")

## MFI -------------------------------------------------------------------------

mfi <- 
  function(.contaminated){
    map_dfr(
      c("multi", "multi2"), 
      function(.setup){
        tau_out <- 
          map_dfr(
            c(0.1, 0.3, 0.5, 0.7, 0.9),
            function(.tau){
              
              error_out <- 
                map_dfr(
                  unique(sim_res_al1brq$parameter$error), 
                  function(.error){
                    
                    method_out <- 
                      map_dfr(
                        c("QR-QR.Boost", "QR-LS.Boost"), 
                        function(.method){
                          
                          betas <- map_chr(
                            1:6, 
                            function(.x){
                              paste0("hatbeta_tau",.x)
                            }
                          )
                          
                          beta_out <- 
                            map_dfc(
                              betas,
                              function(.beta){
                                
                                init_set <- 0
                                
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
                                  
                                  mfi <- 
                                    map_dbl(
                                      results$value, 
                                      function(.res){
                                        mfi_out <- .res[["first.appearance"]][[paste0("beta_",beta_ind)]]
                                        if(is.infinite(mfi_out)){
                                          mfi_out <- NA_real_
                                        }
                                        mfi_out
                                      }
                                    ) %>% 
                                    mean(na.rm = TRUE) %>% 
                                    round(3)
                                  
                                  out <- 
                                    list(mfi)
                                  
                                  names(out) <- .beta
                                  
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
                                  
                                  mfi <- 
                                    map_dbl(
                                      results$value, 
                                      function(.res){
                                        mfi_out <- .res[["first.appearance"]][[paste0("beta_",beta_ind)]]
                                        if(is.infinite(mfi_out)){
                                          mfi_out <- NA_real_
                                        }
                                        mfi_out
                                      }
                                    ) %>% 
                                    mean(na.rm = TRUE) %>% 
                                    round(3)
                                  
                                  out <- 
                                    list(mfi)
                                  
                                  names(out) <- .beta
                                  
                                  return(out)
                                }
                                
                              }
                            )
                          
                          tibble(
                            Method = .method,
                            beta_out
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
                tau = .tau,
                error_out
              )
            }
          )
        tibble(
          `Parameter setup` = .setup,
          tau_out
        )
      }
    )
    
  }

mfi_not_contaminated <- mfi(.contaminated = FALSE)

mfi_not_contaminated %>% 
  filter(!`Error distribution` %in% c("tdist_1", "mixed")) %>% 
  mutate(
    across(where(is.factor), as.character)
  ) %>% 
  mutate(
    across(
      .cols = c(`Parameter setup`, tau, `Error distribution`),
      .fns = function(.x){
        unfill_vec(.x)
      }
    )
  ) %>% 
  kbl("latex", escape = FALSE, linesep = "")

mfi_contaminated <- mfi(.contaminated = TRUE)

mfi_contaminated %>% 
  filter(!`Error distribution` %in% c("tdist_1", "mixed")) %>% 
  mutate(
    across(where(is.factor), as.character)
  ) %>% 
  mutate(
    across(
      .cols = c(`Parameter setup`, tau, `Error distribution`),
      .fns = function(.x){
        unfill_vec(.x)
      }
    )
  ) %>% 
  kbl("latex", escape = FALSE, linesep = "")

## PER -------------------------------------------------------------------------

per <- 
  function(.contaminated){
    map_dfr(
      c("multi", "multi2"), 
      function(.setup){
        tau_out <- 
          map_dfr(
            c(0.1, 0.3, 0.5, 0.7, 0.9),
            function(.tau){
              
              error_out <- 
                map_dfr(
                  unique(sim_res_al1brq$parameter$error), 
                  function(.error){
                    
                    method_out <- 
                      map_dfr(
                        c("QR-QR.Boost", "QR-LS.Boost", "RQAic"), 
                        function(.method){
                          
                          betas <- map_chr(
                            1:6, 
                            function(.x){
                              paste0("hatbeta_tau",.x)
                            }
                          )
                          
                          beta_out <- 
                            map_dfc(
                              betas,
                              function(.beta){
                                
                                init_set <- 0
                                
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
                                  
                                  mean.appearance <-
                                    map_dbl(
                                      results$value, 
                                      function(.res){
                                        .res[["mean.appearance"]][[beta_ind]]
                                      }
                                    )
                                  
                                  excluded_runs <- 
                                    (sum(mean.appearance == 0)/100) %>% 
                                    round(3)
                                  
                                  out <- 
                                    list(excluded_runs)
                                  
                                  names(out) <- .beta
                                  
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
                                  
                                  mean.appearance <-
                                    map_dbl(
                                      results$value, 
                                      function(.res){
                                        .res[["mean.appearance"]][[beta_ind]]
                                      }
                                    )
                                  
                                  excluded_runs <- 
                                    (sum(mean.appearance == 0)/100) %>% 
                                    round(3)
                                  
                                  out <- 
                                    list(excluded_runs)
                                  
                                  names(out) <- .beta
                                  
                                  return(out)
                                }
                                
                                if(.method == "RQAic"){
                                  
                                  results <- 
                                    sim_res_rq$output %>% 
                                    filter(
                                      setup == .setup, 
                                      error == .error, 
                                      tau == .tau, 
                                      contaminated == .contaminated
                                    )
                                  
                                  estimates <- 
                                    map_dbl(
                                      results$value, 
                                      function(.res){
                                        .res[["coefs_rq_bs"]][beta_ind + 1] %>% as.numeric()
                                      }
                                    )
                                  
                                  excluded_runs <- 
                                    (sum(estimates == 0)/100) %>% 
                                    round(3)
                                  
                                  out <- 
                                    list(excluded_runs)
                                  
                                  names(out) <- .beta
                                  
                                  return(out)
                                  
                                }
                                
                              }
                            )
                          
                          tibble(
                            Method = .method,
                            beta_out
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
                tau = .tau,
                error_out
              )
            }
          )
        tibble(
          `Parameter setup` = .setup,
          tau_out
        )
      }
    )
    
  }

per_not_contaminated <- per(.contaminated = FALSE)

per_not_contaminated %>% 
  filter(`Error distribution` != "tdist_1") %>% 
  filter(`Parameter setup` == "multi") %>% 
  group_by(`Parameter setup`, tau, `Error distribution`) %>% 
  mutate(
    across(
      .cols = c(hatbeta_tau1, hatbeta_tau2, hatbeta_tau3, hatbeta_tau4),
      .fns = function(.x){
        if_else(.x == min(.x), str_c("\\blue{", format(.x, digits = 3, nsmall = 3), "}"), as.character(format(.x, digits = 3, nsmall = 3)))
      }
    )
  ) %>% 
  mutate(
    across(
      .cols = c(hatbeta_tau5, hatbeta_tau6),
      .fns = function(.x){
        if_else(.x == max(.x), str_c("\\blue{", format(.x, digits = 3, nsmall = 3), "}"), as.character(format(.x, digits = 3, nsmall = 3)))
      }
    )
  ) %>% 
  ungroup() %>% 
  mutate(
    across(where(is.factor), as.character)
  ) %>% 
  mutate(
    across(
      .cols = c(`Parameter setup`, tau, `Error distribution`),
      .fns = function(.x){
        unfill_vec(.x)
      }
    )
  ) %>% 
  kbl("latex", escape = FALSE, linesep = "")

per_not_contaminated %>% 
  filter(`Error distribution` != "tdist_1") %>% 
  filter(`Parameter setup` == "multi2") %>% 
  group_by(`Parameter setup`, tau, `Error distribution`) %>% 
  mutate(
    across(
      .cols = c(hatbeta_tau1),
      .fns = function(.x){
        if_else(
          tau >= 0.5,
          if_else(.x == max(.x), str_c("\\blue{", format(.x, digits = 3, nsmall = 3), "}"), as.character(format(.x, digits = 3, nsmall = 3))),
          if_else(.x == min(.x), str_c("\\blue{", format(.x, digits = 3, nsmall = 3), "}"), as.character(format(.x, digits = 3, nsmall = 3)))
        )
      }
    )
  ) %>% 
  mutate(
    across(
      .cols = c(hatbeta_tau3),
      .fns = function(.x){
        if_else(
          tau <= 0.5,
          if_else(.x == max(.x), str_c("\\blue{", format(.x, digits = 3, nsmall = 3), "}"), as.character(format(.x, digits = 3, nsmall = 3))),
          if_else(.x == min(.x), str_c("\\blue{", format(.x, digits = 3, nsmall = 3), "}"), as.character(format(.x, digits = 3, nsmall = 3)))
        )
      }
    )
  ) %>% 
  mutate(
    across(
      .cols = c(hatbeta_tau2, hatbeta_tau4),
      .fns = function(.x){
        if_else(.x == min(.x), str_c("\\blue{", format(.x, digits = 3, nsmall = 3), "}"), as.character(format(.x, digits = 3, nsmall = 3)))
      }
    )
  ) %>% 
  mutate(
    across(
      .cols = c(hatbeta_tau5, hatbeta_tau6),
      .fns = function(.x){
        if_else(.x == max(.x), str_c("\\blue{", format(.x, digits = 3, nsmall = 3), "}"), as.character(format(.x, digits = 3, nsmall = 3)))
      }
    )
  ) %>% 
  ungroup() %>% 
  mutate(
    across(where(is.factor), as.character)
  ) %>% 
  mutate(
    across(
      .cols = c(`Parameter setup`, tau, `Error distribution`),
      .fns = function(.x){
        unfill_vec(.x)
      }
    )
  ) %>% 
  kbl("latex", escape = FALSE, linesep = "")

per_contaminated <- per(.contaminated = TRUE)

per_contaminated %>% 
  filter(`Error distribution` != "tdist_1") %>% 
  filter(`Parameter setup` == "multi") %>% 
  group_by(`Parameter setup`, tau, `Error distribution`) %>% 
  mutate(
    across(
      .cols = c(hatbeta_tau1, hatbeta_tau2, hatbeta_tau3, hatbeta_tau4),
      .fns = function(.x){
        if_else(.x == min(.x), str_c("\\blue{", format(.x, digits = 3, nsmall = 3), "}"), as.character(format(.x, digits = 3, nsmall = 3)))
      }
    )
  ) %>% 
  mutate(
    across(
      .cols = c(hatbeta_tau5, hatbeta_tau6),
      .fns = function(.x){
        if_else(.x == max(.x), str_c("\\blue{", format(.x, digits = 3, nsmall = 3), "}"), as.character(format(.x, digits = 3, nsmall = 3)))
      }
    )
  ) %>% 
  ungroup() %>% 
  mutate(
    across(where(is.factor), as.character)
  ) %>% 
  mutate(
    across(
      .cols = c(`Parameter setup`, tau, `Error distribution`),
      .fns = function(.x){
        unfill_vec(.x)
      }
    )
  ) %>% 
  kbl("latex", escape = FALSE, linesep = "")

per_contaminated %>% 
  filter(`Error distribution` != "tdist_1") %>% 
  filter(`Parameter setup` == "multi2") %>% 
  group_by(`Parameter setup`, tau, `Error distribution`) %>% 
  mutate(
    across(
      .cols = c(hatbeta_tau1),
      .fns = function(.x){
        if_else(
          tau >= 0.5,
          if_else(.x == max(.x), str_c("\\blue{", format(.x, digits = 3, nsmall = 3), "}"), as.character(format(.x, digits = 3, nsmall = 3))),
          if_else(.x == min(.x), str_c("\\blue{", format(.x, digits = 3, nsmall = 3), "}"), as.character(format(.x, digits = 3, nsmall = 3)))
        )
      }
    )
  ) %>% 
  mutate(
    across(
      .cols = c(hatbeta_tau3),
      .fns = function(.x){
        if_else(
          tau <= 0.5,
          if_else(.x == max(.x), str_c("\\blue{", format(.x, digits = 3, nsmall = 3), "}"), as.character(format(.x, digits = 3, nsmall = 3))),
          if_else(.x == min(.x), str_c("\\blue{", format(.x, digits = 3, nsmall = 3), "}"), as.character(format(.x, digits = 3, nsmall = 3)))
        )
      }
    )
  ) %>% 
  mutate(
    across(
      .cols = c(hatbeta_tau2, hatbeta_tau4),
      .fns = function(.x){
        if_else(.x == min(.x), str_c("\\blue{", format(.x, digits = 3, nsmall = 3), "}"), as.character(format(.x, digits = 3, nsmall = 3)))
      }
    )
  ) %>% 
  mutate(
    across(
      .cols = c(hatbeta_tau5, hatbeta_tau6),
      .fns = function(.x){
        if_else(.x == max(.x), str_c("\\blue{", format(.x, digits = 3, nsmall = 3), "}"), as.character(format(.x, digits = 3, nsmall = 3)))
      }
    )
  ) %>% 
  ungroup() %>% 
  mutate(
    across(where(is.factor), as.character)
  ) %>% 
  mutate(
    across(
      .cols = c(`Parameter setup`, tau, `Error distribution`),
      .fns = function(.x){
        unfill_vec(.x)
      }
    )
  ) %>% 
  kbl("latex", escape = FALSE, linesep = "")

## True population parameters --------------------------------------------------

true_betas <- 
  map_dfr(
    c("multi", "multi2"), 
    function(.setup){
      
      tau_out <- 
        map_dfr(
          c(0.1, 0.3, 0.5, 0.7, 0.9),
          function(.tau){
            
            error_out <- map_dfr(
              unique(sim_res_al1brq$parameter$error) %>% setdiff(c("mixed", "tdist_1")),
              function(.error){
                
                if(.setup == "multi2"){
                  
                  tau.vec <- c(0.1, 0.3, 0.5, 0.7, 0.9)
                  
                  betas <- sim[[.setup]]$beta[which(.tau == tau.vec),]
                  
                  if(.error == "norm"){
                    
                    betatau <- map(seq_along(betas), function(.x){
                      betas[.x] + sim[[.setup]]$alpha[.x] * qnorm(p = .tau)
                    })
                    
                  } else if(.error == "tdist"){
                    
                    betatau <- map(seq_along(betas), function(.x){
                      betas[.x] + sim[[.setup]]$alpha[.x] * qt(p = .tau, df = 2)
                    })
                    
                  } else if(.error == "gamma"){
                    
                    betatau <- map(seq_along(betas), function(.x){
                      betas[.x] + sim[[.setup]]$alpha[.x] * qgamma(p = .tau, shape = 2, scale = 1)
                    })
                    
                  } else if(.error == "tdist_1"){
                    
                    betatau <- map(seq_along(betas), function(.x){
                      betas[.x] + sim[[.setup]]$alpha[.x] * qt(p = .tau, df = 1)
                    })
                    
                  }
                  
                } else {
                  if(.error == "norm"){
                    
                    betatau <- map(seq_along(sim[[.setup]]$beta), function(.x){
                      sim[[.setup]]$beta[.x] + sim[[.setup]]$alpha[.x] * qnorm(p = .tau)
                    })
                    
                  } else if(.error == "tdist"){
                    
                    betatau <- map(seq_along(sim[[.setup]]$beta), function(.x){
                      sim[[.setup]]$beta[.x] + sim[[.setup]]$alpha[.x] * qt(p = .tau, df = 2)
                    })
                    
                  } else if(.error == "gamma"){
                    
                    betatau <- map(seq_along(sim[[.setup]]$beta), function(.x){
                      sim[[.setup]]$beta[.x] + sim[[.setup]]$alpha[.x] * qgamma(p = .tau, shape = 2, scale = 1)
                    })
                    
                  } else if(.error == "tdist_1"){
                    
                    betatau <- map(seq_along(sim[[.setup]]$beta), function(.x){
                      sim[[.setup]]$beta[.x] + sim[[.setup]]$alpha[.x] * qt(p = .tau, df = 1)
                    })
                    
                  }
                }
                
                names(betatau) <- paste0("betatau",0:6)
                
                tibble(
                  `Error distribution` = .error, 
                  map_dfc(betatau, ~.x)
                )
              }
            )
            
            tibble(
              Tau = .tau,
              error_out
            )
            
          }
        )
      
      tibble(
        `Parameter setup` = .setup, 
        tau_out
      )
    }
  )

true_betas


## Computational times ---------------------------------------------------------

iterations <- 
  function(.contaminated){
    map_dfr(
      c("hom", "het", "multi", "multi2", "highdim"), 
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
  mutate(
    across(where(is.factor), as.character)
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
  mutate(
    across(where(is.factor), as.character)
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
    across(where(is.factor), as.character)
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
    across(where(is.factor), as.character)
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
    across(where(is.factor), as.character)
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
    across(where(is.factor), as.character)
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
sim_sample_data_rq <- readRDS("simulation_output/rq_sample_data.RDS")

tau_insample <- 
  function(.contaminated){
    map_dfr(
      c("hom", "het", "multi", "multi2", "highdim"), 
      function(.setup){
        error_out <- map_dfr(
          unique(sim_res_al1brq$parameter$error),
          function(.error){
            
            if(.setup %in% c("multi", "multi2")){
              methods <- c("QR-QR.Boost", "QR-LS.Boost", "RQ", "RQAic")
            } else {
              methods <- c("QR-QR.Boost", "QR-LS.Boost", "RQ")
            }
            
            method_out <- map_dfr(
              methods,
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
                      
                      if(.method == "RQ"){
                        
                        sample_data <- 
                          sim_sample_data_rq$output %>% 
                          filter(
                            setup == .setup, 
                            error == .error, 
                            tau == .tau, 
                            contaminated == .contaminated
                          )
                        
                        results <- 
                          sim_res_rq$output %>% 
                          filter(
                            setup == .setup, 
                            error == .error, 
                            tau == .tau, 
                            contaminated == .contaminated
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
                              
                              beta <- results$value[[.x]][["coefs_rq"]] %>% as.numeric()
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
                      
                      if(.method == "RQAic"){
                        
                        sample_data <- 
                          sim_sample_data_rq$output %>% 
                          filter(
                            setup == .setup, 
                            error == .error, 
                            tau == .tau, 
                            contaminated == .contaminated
                          )
                        
                        results <- 
                          sim_res_rq$output %>% 
                          filter(
                            setup == .setup, 
                            error == .error, 
                            tau == .tau, 
                            contaminated == .contaminated
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
                              
                              beta <- results$value[[.x]][["coefs_rq_bs"]] %>% as.numeric()
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
        if_else(.x == max(.x), str_c("\\blue{", format(.x, digits = 3, nsmall = 3), "}"), as.character(format(.x, digits = 3, nsmall = 3)))
      }
    )
  ) %>% 
  ungroup() %>% 
  mutate(
    across(where(is.factor), as.character)
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
        if_else(.x == max(.x), str_c("\\blue{", format(.x, digits = 3, nsmall = 3), "}"), as.character(format(.x, digits = 3, nsmall = 3)))
      }
    )
  ) %>% 
  ungroup() %>% 
  mutate(
    across(where(is.factor), as.character)
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
