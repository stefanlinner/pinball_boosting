# Sourcing necessary scripts ---------------------------------------------------

source("libraries.R", encoding = "UTF-8")
source("simulation_setups.R", encoding = "UTF-8")


revision <- TRUE # FALSE for original results

# simulation parameter combinations --------------------------------------------

parameter_combinations <- 
  list(
    setup = c("hom", "het", "multi", "highdim", "multi2"),
    error = c("norm", "tdist", "gamma", "mixed", "tdist_1"),
    tau = c(0.1, 0.3, 0.5, 0.7, 0.9),
    contaminated = c(TRUE, FALSE)
  )

parameter_combinations_df <- 
  expand.grid(parameter_combinations)

parameter_combinations_boosting <- 
  parameter_combinations_df

parameter_combinations_boosting$init <- 0

parameter_combinations_boosting <- 
  rbind(
    parameter_combinations_boosting,
    expand.grid(
      list(
        setup = "hom",
        error = "norm",
        tau = c(0.1, 0.3, 0.5, 0.7, 0.9), 
        contaminated = c(TRUE, FALSE),
        init = c(0.1, 0.3, 0.5, 0.7, 0.9)
      )
    )
  )



# test data --------------------------------------------------------------------

if(!revision){
  test_data <- 
    map(
      seq_len(nrow(parameter_combinations_df)),
      function(row){
        sim[[parameter_combinations_df[row, "setup"]]]$data_gen(
          seed = 101, 
          n = 1000, 
          error = parameter_combinations_df[row, "error"], 
          tau = parameter_combinations_df[row, "tau"], 
          contaminated = parameter_combinations_df[row, "contaminated"],
          revision = revision
        )
      }
    ) %>% 
    set_names(
      paste(
        parameter_combinations_df[["setup"]],
        parameter_combinations_df[["error"]],
        parameter_combinations_df[["tau"]],
        parameter_combinations_df[["contaminated"]], 
        sep = ","
      )
    )
} else {
  test_data <- 
    map(
      seq_len(nrow(parameter_combinations_df)),
      function(row){
        sim[[parameter_combinations_df[row, "setup"]]]$data_gen(
          seed = 101, 
          n = 100, 
          error = parameter_combinations_df[row, "error"], 
          tau = parameter_combinations_df[row, "tau"], 
          contaminated = parameter_combinations_df[row, "contaminated"],
          revision = revision
        )
      }
    ) %>% 
    set_names(
      paste(
        parameter_combinations_df[["setup"]],
        parameter_combinations_df[["error"]],
        parameter_combinations_df[["tau"]],
        parameter_combinations_df[["contaminated"]], 
        sep = ","
      )
    )
}



# simulation functions ---------------------------------------------------------

## AL1BRQ simulation function --------------------------------------------------

sim_fun_al1brq <- function(
    setup,
    error,
    tau,
    init,
    nu, 
    exact.fit, 
    contaminated,
    sim_param, 
    test_data, 
    revision
){
  
  if(init == 0){
    init <- tau
  }
  
  test.data <- test_data[[paste(setup, error, tau, contaminated, sep = ",")]]
  
  z.brq <- sim_param[[setup]]$zbrq
  n <-  sim_param[[setup]]$n
  p <- sim_param[[setup]]$p
  covars <- all.vars(sim_param[[setup]]$formulabrq()[[3]])
  
  beta_names <- map_chr(seq_len(p), function(.x){
    paste("beta_", .x - 1, sep = "")
  })
  
  brq_names <- map_chr(covars, function(.x){
    paste("brq(", .x, ")", sep = "")
  })
  
  emp.risks <- vector(mode = "numeric", length = length(z.brq))
  
  
  train.data <- sim_param[[setup]]$data_gen(
    n = n, 
    error = error, 
    tau = tau, 
    contaminated = contaminated,
    revision = revision
  )
  
  brq.model <- 
    boostrq(
      formula = sim_param[[setup]]$formulabrq(), 
      data = train.data,
      mstop = z.brq[1], 
      nu = nu,
      tau = tau, 
      digits = 10, 
      offset = rep(sim_param[[setup]]$offset(
        y = train.data$y, 
        exact.fit = exact.fit, 
        init = init
      ), n), 
      exact.fit = exact.fit
    )
  
  emp.risks <- 
    map_dbl(seq_len(length(z.brq)), function(.j){
      
      brq.model$subset(z.brq[.j])
      
      y.test.pred <- predict(brq.model, test.data)
      
      round(
        sum(
          tau * (test.data$y - y.test.pred) * ((test.data$y - y.test.pred) > 0) + 
            (tau - 1) * (test.data$y - y.test.pred) * ((test.data$y - y.test.pred) <= 0)
        ), 
        3)
      
    })
  
  brq.mstop <- z.brq[which.min(emp.risks)]
  
  brq.model$subset(brq.mstop)
  
  num.slope.iters <- brq.mstop - sum(brq.model$xselect() == 0)
  
  coefs <- 
    map(
      .x = seq_len(p),
      function(.x){
        if(beta_names[.x] == "beta_0"){
          brq.model$offset + map_dbl(brq_names, function(.y){
            coef(brq.model)[[.y]][1]
          }) %>% sum()
        } else {
          coef(brq.model)[[brq_names[.x - 1]]][2]
        }
      }
    ) %>% 
    set_names(beta_names)
  
  mean.appearance <- map(
    seq_len(p-1), 
    function(.x){
      sum(brq.model$xselect() == .x) / num.slope.iters
    }
  ) %>% 
    set_names(beta_names[-1])
  
  first.appearance <- map(
    seq_len(p-1), 
    function(.x){
      min(which(brq.model$xselect()[brq.model$xselect() != 0] == .x)) / num.slope.iters
    }
  ) %>% 
    set_names(beta_names[-1])
  
  list(
    coefs = coefs,
    emp.risks = emp.risks,
    mean.appearance = mean.appearance,
    first.appearance = first.appearance
  )
  
}


# L2BRQ simulation function ----------------------------------------------------

sim_fun_l2brq <- 
  function(
    setup,
    error,
    tau,
    init,
    nu, 
    exact.fit, 
    contaminated,
    sim_param, 
    test_data, 
    revision
  ){
    
    if(init == 0 & !revision){
      init <- 0.5
    } 
    
    if(init == 0 & revision){
      init <- tau
    }
    
    if(revision){
      nu <- -0.05 * abs(tau - 0.5) + 0.11
    } else {
      nu <- 0.1
    }
    
    test.data <- test_data[[paste(setup, error, tau, contaminated, sep = ",")]]
    
    n <-  sim_param[[setup]]$n
    p <- sim_param[[setup]]$p
    covars <- all.vars(sim_param[[setup]]$formulabrq()[[3]])
    
    beta_names <- map_chr(seq_len(p), function(.x){
      paste("beta_", .x - 1, sep = "")
    })
    
    z.qb <- sim_param[[setup]]$zqb
    
    qb_names <- map_chr(covars, function(.x){
      paste("bols(", .x, ")", sep = "")
    })
    
    emp.risks <- vector(mode = "numeric", length = length(z.qb))
    
    
    train.data <- 
      sim_param[[setup]]$data_gen(
        n = n, 
        error = error, 
        tau = tau, 
        contaminated = contaminated,
        revision = revision
      )
    
    qb.model <- 
      mboost(
        formula = sim_param[[setup]]$formulaqb(),
        family = QuantReg(tau = c(tau)), 
        control = boost_control(mstop = z.qb[1], nu = nu), 
        offset = sim_param[[setup]]$offset(
          y = train.data$y, 
          exact.fit = exact.fit, 
          init = init
        ), 
        data = train.data
      )
    
    emp.risks <- 
      map_dbl(seq_len(length(z.qb)), function(.j){
        
        qb.model$subset(z.qb[.j])
        
        y.test.pred <- predict(qb.model, test.data)
        
        round(
          sum(
            tau * (test.data$y - y.test.pred) * ((test.data$y - y.test.pred) > 0) + 
              (tau - 1) * (test.data$y - y.test.pred) * ((test.data$y - y.test.pred) <= 0)
          ), 
          3)
        
      })
    
    qb.mstop <- z.qb[which.min(emp.risks)]
    
    qb.model$subset(qb.mstop)
    
    coefs <- 
      map(
        .x = seq_len(p),
        function(.x){
          if(beta_names[.x] == "beta_0"){
            qb.model$offset + map_dbl(qb_names, function(.y){
              coef(qb.model, which = "")[[.y]][1]
            }) %>% sum()
          } else{
            coef(qb.model, which = "")[[qb_names[.x - 1]]][2]
          }
        }
      ) %>% 
      set_names(beta_names)
    
    mean.appearance <- map(
      seq_len(p-1), 
      function(.x){
        sum(qb.model$xselect() == .x) / qb.mstop
      }
    ) %>% 
      set_names(beta_names[-1])
    
    first.appearance <- map(
      seq_len(p-1), 
      function(.x){
        min(which(qb.model$xselect() == .x)) / qb.mstop
      }
    ) %>% 
      set_names(beta_names[-1])
    
    list(
      coefs = coefs,
      emp.risks = emp.risks,
      mean.appearance = mean.appearance,
      first.appearance = first.appearance
    )
    
  }


## RQ & RQAic simulation function ----------------------------------------------


sim_fun_rq <- 
  function(
    setup,
    error,
    tau,
    contaminated,
    sim_param,
    revision
  ){
    
    n <-  sim_param[[setup]]$n
    p <- sim_param[[setup]]$p
    covars <- all.vars(sim_param[[setup]]$formulabrq()[[3]])
    
    beta_names <- map_chr(seq_len(p), function(.x){
      paste("beta_", .x - 1, sep = "")
    })
    
    train.data <- sim_param[[setup]]$data_gen(
      n = n,
      error = error, 
      tau = tau, 
      contaminated = contaminated,
      revision = revision 
    )
    
    rq.model <- 
      rq(y ~ ., tau = tau, data = train.data)
    
    coefs <- coef(rq.model)
    names(coefs) <- NULL
    
    
    coefs_rq <- data.frame(
      map_dbl(
        .x = seq_len(p), 
        function(.x){
          coefs[.x]
        }) %>% 
        map(.f = function(.z) .z)
    ) %>% 
      set_names(c(beta_names))
    
    coefs_rq_bs <- NULL
    
    if(setup %in% c("multi", "multi2")){
      
      best.subset <- glmulti(
        y = "y",
        xr = covars,
        data = train.data,
        fitfunction = "rq", 
        crit = "aic", 
        method = "h",
        level = 1, 
        tau = tau, 
        plotty = FALSE,
        report = FALSE
      )
      
      rq.model_bs <- 
        rq(
          as.formula(summary(best.subset)$bestmodel), 
          tau = tau,
          data = train.data
        )
      
      coefs_bs <- coef(rq.model_bs)
      names(coefs_bs) <- NULL
      
      coefs_rq_bs <- 
        data.frame( 
          map_dbl(
            .x = seq_len(p), 
            function(.x){
              coefs_bs[.x]
            }) %>% 
            map(.f = function(.z).z)
        ) %>% 
        set_names(c(beta_names))
      
      coefs_rq_bs[is.na(coefs_rq_bs)] <- 0
      
    }
    
    list(
      coefs_rq = coefs_rq,
      coefs_rq_bs = coefs_rq_bs
    )
    
  }


get_sample_data <- 
  function(
    setup,
    error,
    tau,
    init,
    contaminated,
    sim_param,
    revision
  ) {
    
    n <-  sim_param[[setup]]$n
    
    train.data <- 
      sim_param[[setup]]$data_gen(
        n = n, 
        error = error, 
        tau = tau, 
        contaminated = contaminated,
        revision = revision
      )
    
    train.data
  }


# simulation workers -----------------------------------------------------------

saveRDS("test", here("simulation_output/test_save.RDS"))
message(future::availableCores())

message("Al1brq")

set.seed(123)
sim_al1brq <-
  future_mc(
    fun = sim_fun_al1brq,
    repetitions = 100,
    param_table = parameter_combinations_boosting,
    check = FALSE,
    nu = NULL,
    exact.fit = FALSE,
    sim_param = sim,
    test_data = test_data, 
    revision = revision
  )

saveRDS(sim_al1brq, here("simulation_output/simulation_results_rev_al1brq.RDS"))
rm(sim_al1brq)

message("l2brq")

set.seed(123)
sim_l2brq <-
  future_mc(
    fun = sim_fun_l2brq,
    repetitions = 100,
    param_table = parameter_combinations_boosting,
    check = FALSE,
    nu = NULL,
    exact.fit = FALSE,
    sim_param = sim,
    test_data = test_data, 
    revision = revision
  )

saveRDS(sim_l2brq, here("simulation_output/simulation_results_rev_l2brq.RDS"))
rm(sim_l2brq)

# set.seed(123)
# sim_rq <-
#   future_mc(
#     fun = sim_fun_rq,
#     repetitions = 100,
#     param_list =  parameter_combinations,
#     check = FALSE,
#     sim_param = sim,
#     parallel = FALSE, 
#     revision = revision
#   )
# 
# saveRDS(sim_rq, "simulation_output/simulation_results_rq.RDS")
# rm(sim_rq)
# 

message("data")

set.seed(123)
sim_sample_data_boosting <-
  future_mc(
    fun = get_sample_data,
    repetitions = 100,
    param_table = parameter_combinations_boosting,
    check = FALSE,
    sim_param = sim,
    revision = revision
  )
saveRDS(sim_sample_data_boosting, here("simulation_output/boosting_sample_data_rev.RDS"))
rm(sim_sample_data_boosting)

message("done")
# 
# set.seed(123)
# sim_sample_data_rq <-
#   future_mc(
#     fun = get_sample_data,
#     repetitions = 100,
#     param_list = parameter_combinations,
#     check = FALSE,
#     sim_param = sim,
#     init = 0,
#     parallel = FALSE, 
#     revision = revision
#   )
# saveRDS(sim_sample_data_rq, "simulation_output/rq_sample_data.RDS")
# rm(sim_sample_data_rq)
