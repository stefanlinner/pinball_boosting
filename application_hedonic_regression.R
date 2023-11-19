# Sourcing libraries -----------------------------------------------------------

source("libraries.R", encoding = "UTF-8")


# Reading data -----------------------------------------------------------------

dat <- read.table(
  "ag-data.fil",
  col.names = c(
    "sell", "lot", "bdms", "fb", "sty", "drv", "rec", "ffin", "ghw", "ca", "gar", "reg"
  )
)

dat_interactions <- dat

dat_interactions$bdms <- 
  factor(
    dat_interactions$bdms, 
    levels = 1:6
  )

dat_interactions$fb <- 
  factor(
    dat_interactions$fb, 
    levels = 1:4
  )

dat_interactions$sty <- 
  factor(
    dat_interactions$sty, 
    levels = 1:4
  )


# Model matrices ---------------------------------------------------------------


model_fb <- model.matrix(
  ~ - 1 + fb, 
  data = dat_interactions
) %>% 
  as.data.frame()

model_sty <- model.matrix(
  ~ -1 + sty, 
  data = dat_interactions
) %>% 
  as.data.frame()

model_bdms <- model.matrix(
  ~ -1 + bdms, 
  data = dat_interactions
) %>% 
  as.data.frame()

model_lot <- 
  model.matrix(
    ~ -1 + log(lot) + I(log(lot)^2), 
    data = dat_interactions
  ) %>% 
  as.data.frame()

names(model_lot) <- c("llot", "lot2")

model_fb_sty <- 
  map_dfc(
    names(model_fb),
    function(fb){
      map_dfc(
        names(model_sty), 
        function(sty){
          out <- list(model_fb[[fb]] * model_sty[[sty]])
          names(out) <- paste(fb, sty, sep =".")
          out
        }
      )
    }
  )

model_fb_bdms <- 
  map_dfc(
    names(model_fb),
    function(fb){
      map_dfc(
        names(model_bdms), 
        function(bdms){
          out <- list(model_fb[[fb]] * model_bdms[[bdms]])
          names(out) <- paste(fb, bdms, sep =".")
          out
        }
      )
    }
  )

model_sty_bdms <- 
  map_dfc(
    names(model_sty),
    function(sty){
      map_dfc(
        names(model_bdms), 
        function(bdms){
          out <- list(model_sty[[sty]] * model_bdms[[bdms]])
          names(out) <- paste(sty, bdms, sep =".")
          out
        }
      )
    }
  )

model_cat_interactions <- 
  cbind(model_fb_sty, model_fb_bdms, model_sty_bdms)

model_lot_interactions <- 
  map_dfc(
    names(model_cat_interactions),
    function(cat){
      map_dfc(
        names(model_lot), 
        function(lot){
          out <- list(model_cat_interactions[[cat]] * model_lot[[lot]])
          names(out) <- paste(cat, lot, sep =".")
          out
        }
      )
    }
  )

model_full_interactions <- 
  cbind(
    model_fb,
    model_bdms, 
    model_sty, 
    model_lot,
    model_cat_interactions, 
    model_lot_interactions, 
    dat_interactions[, c("sell", "drv", "rec", "ffin", "ghw", "ca", "gar", "reg")]
  )


# Model formula ----------------------------------------------------------------

formula_brq <- 
  map_chr(
    names(model_full_interactions)[!(names(model_full_interactions) %in% c("sell"))], 
    function(var){
      paste("brq(", var, ")", sep = "")
    }
  ) %>% 
  paste(collapse = " + ") %>% 
  paste("log(sell) ~", .) %>% 
  as.formula()

highdim <- 
  function(taus){
    map_dfr(
      taus, 
      function(tau){
        
        al1brq_small <- boostrq(
          log(sell) ~ brq(drv) + brq(rec) + brq(ffin) + brq(ghw) + brq(ca) + brq(gar) + 
            brq(reg) + brq(log(lot)) + brq(bdms) + brq(fb) + brq(sty),
          data = dat, 
          tau = tau, 
          nu = NULL, 
          mstop = 1
        )
        
        set.seed(101)
        cvkbrq_small <-
          cvrisk(
            al1brq_small,
            grid = seq(10, 500, 10),
            folds = mboost::cv(al1brq_small$weights, type = "kfold", B = 5)
          )
        
        if(tau == 0.1){
          risk_al1brq_small <- mean(cvkbrq_small[,names(which.min(round(colMeans(cvkbrq_small),8)))]) %>% round(3)
        } else {
          risk_al1brq_small <- mean(cvkbrq_small[,as.character(mstop(cvkbrq_small))]) %>% round(3)
        }
        
        al1brq_large <- 
          boostrq(
            formula_brq,
            data = model_full_interactions, 
            tau = tau, 
            nu = NULL, 
            mstop = 1
          )
        
        set.seed(101)
        cvkbrq_large <-
          cvrisk(
            al1brq_large,
            grid = seq(10, 500, 10),
            folds = mboost::cv(al1brq_large$weights, type = "kfold", B = 5)
          )
        
        mstop(al1brq_large) <- mstop(cvkbrq_large)
        
        risk_al1brq_large <- mean(cvkbrq_large[,as.character(mstop(cvkbrq_large))]) %>% round(3)
        
        not_selected <- 
          setdiff(1:215, al1brq_large$xselect() %>% unique()) %>% 
          length()
        
        data.frame(
          tau = tau, 
          AL1BRQ_small = risk_al1brq_small, 
          AL1BRQ_large = risk_al1brq_large,
          Not_selected = not_selected
        )
        
        
      }
    ) %>% 
      kable(
        format = "latex", 
        caption = "Risk AL1BRQ"
      )
  }


aggr_coef <- 
  function(model){
    
    coefs <- 
      lapply(coef(model), 
             function(coefs){
               coefs[2]
             }) %>%
      unlist()
    
    names(coefs) <- names(coef(model))
    
    intercept <- 
      lapply(coef(model), 
             function(coefs){
               coefs[1]
             }) %>%
      unlist() %>% 
      sum() %>% 
      sum(attributes(coef(model))$offset)
    
    
    coefs <- c("(Intercept)" = intercept, coefs)
    coefs
  }

k_fold_cv <- 
  function(tau, method = "rq"){
    
    folds <- mboost::cv(rep(1, nrow(dat)), type = "kfold", B = 5)
    
    lapply(
      seq_len(ncol(folds)), 
      function(fold){
        
        if(method == "rq"){
          model <- rq(
            log(sell) ~ drv + rec+ ffin + ghw + ca + gar + reg + log(lot) + bdms + fb + sty, 
            tau = tau,
            data = dat, 
            weights = folds[, fold]
          )
        }
        
        if(method == "ols"){
          model <- lm(
            log(sell) ~ drv + rec+ ffin + ghw + ca + gar + reg + log(lot) + bdms + fb + sty, 
            data = dat, 
            weights = folds[, fold]
          )
        }
        
        preds <- predict(model, newdata = dat[(folds[, fold] - 1) * (-1), ])
        
        boostrq:::quantile.risk(
          y = log(dat[(folds[, fold] - 1) * (-1), "sell"]), 
          f = preds, 
          weights = rep(1, length(preds)), 
          tau = tau
        ) / sum((folds[, fold] - 1) * (-1))
      }
    ) %>% 
      unlist() %>% 
      mean()
  }

coef_table <- 
  function(tau){
    
    rq <- rq(
      log(sell) ~ drv + rec+ ffin + ghw + ca + gar + reg + log(lot) + bdms + fb + sty, 
      tau = tau,
      data = dat
    )
    
    al1brq <- boostrq(
      log(sell) ~ brq(drv) + brq(rec) + brq(ffin) + brq(ghw) + brq(ca) + brq(gar) + 
        brq(reg) + brq(log(lot)) + brq(bdms) + brq(fb) + brq(sty),
      data = dat, 
      tau = tau, 
      nu = NULL, 
      mstop = 1
    )
    
    set.seed(101)
    cvkbrq <-
      cvrisk(
        al1brq,
        grid = seq(10, 500, 10),
        folds = mboost::cv(al1brq$weights, type = "kfold", B = 5)
      )
    
    if(tau == 0.1){
      mstop(al1brq) <- as.numeric(names(which.min(round(colMeans(cvkbrq),8))))
    } else {
      mstop(al1brq) <- mstop(cvkbrq)
    }
    
    l2brq <- mboost(
      log(sell) ~ drv + rec + ffin + ghw + ca + gar + 
        reg + log(lot) + bdms + fb + sty,
      data = dat, 
      baselearner = "bols",
      family = QuantReg(tau = tau, qoffset = 0.5), 
      control = boost_control(mstop = 200, nu = 0.1)
    )
    
    set.seed(101)
    cvkl2 <-
      cvrisk(
        l2brq,
        grid = seq(10, 500, 10),
        folds = mboost::cv(model.weights(l2brq), type = "kfold", B = 5)
      )
    
    mstop(l2brq) <- mstop(cvkl2)
    
    if(tau == 0.5){
      
      ols <- lm(
        log(sell) ~ drv + rec+ ffin + ghw + ca + gar + reg + log(lot) + bdms + fb + sty, 
        data = dat
      )
      
      l2brm <-  mboost(
        log(sell) ~ drv + rec + ffin + ghw + ca + gar + 
          reg + log(lot) + bdms + fb + sty,
        data = dat, 
        baselearner = "bols",
        family = Gaussian(), 
        control = boost_control(mstop = 200, nu = 0.1)
      )
      
      set.seed(101)
      cvkbrm <-
        cvrisk(
          l2brm,
          grid = seq(10, 500, 10),
          folds = mboost::cv(model.weights(l2brm), type = "kfold", B = 5)
        )
      
      mstop(l2brm) <- mstop(cvkbrm)
      
      out <- 
        data.frame(
          RQ = coef(rq) %>% round(3), 
          AL1BRQ = coef(al1brq, aggregate = "sum_aggr") %>% round(3), 
          L2BRQ = aggr_coef(l2brq) %>% round(3), 
          OLS = coef(ols) %>% round(3),
          L2BRM = aggr_coef(l2brm) %>% round(3)
        ) %>% 
        kable(
          format = "latex", 
          row.names = TRUE, 
          caption = paste("tau =", tau)
        )
      
      return(out)
      
    } else{
      out <- data.frame(
        RQ = coef(rq) %>% round(3), 
        AL1BRQ = coef(al1brq, aggregate = "sum_aggr") %>% round(3), 
        L2BRQ = aggr_coef(l2brq) %>% round(3)
      ) %>% 
        kable(
          format = "latex", 
          row.names = TRUE, 
          caption = paste("tau =", tau)
        )
      return(out)
    }
    
    data.frame(
      RQ = coef(rq) %>% round(3), 
      AL1BRQ = coef(al1brq, aggregate = "sum_aggr") %>% round(3), 
      L2BRQ = aggr_coef(l2brq) %>% round(3)
    ) %>% 
      kable(
        format = "latex", 
        row.names = TRUE, 
        caption = paste("tau =", tau)
      )
  }

risk_table <- 
  function(taus){
    
    map_dfr(
      taus,
      function(tau){
        
        
        set.seed(101)
        rq_risk <- k_fold_cv(tau) %>% round(3)
        
        al1brq <- boostrq(
          log(sell) ~ brq(drv) + brq(rec) + brq(ffin) + brq(ghw) + brq(ca) + brq(gar) + 
            brq(reg) + brq(log(lot)) + brq(bdms) + brq(fb) + brq(sty),
          data = dat, 
          tau = tau, 
          nu = NULL, 
          mstop = 1
        )
        
        set.seed(101)
        cvkbrq <-
          cvrisk(
            al1brq,
            grid = seq(10, 500, 10),
            folds = mboost::cv(al1brq$weights, type = "kfold", B = 5)
          )
        
        if(tau == 0.1){
          risk_al1brq <- mean(cvkbrq[,names(which.min(round(colMeans(cvkbrq),8)))]) %>% round(3)
        } else {
          risk_al1brq <- mean(cvkbrq[,as.character(mstop(cvkbrq))]) %>% round(3)
        }
        
        l2brq <- mboost(
          log(sell) ~ drv + rec + ffin + ghw + ca + gar + 
            reg + log(lot) + bdms + fb + sty,
          data = dat, 
          baselearner = "bols",
          family = QuantReg(tau = tau, qoffset = 0.5), 
          control = boost_control(mstop = 200, nu = 0.1)
        )
        
        set.seed(101)
        cvkl2 <-
          cvrisk(
            l2brq,
            grid = seq(10, 500, 10),
            folds = mboost::cv(model.weights(l2brq), type = "kfold", B = 5)
          )
        
        risk_l2brq <- mean(cvkl2[,as.character(mstop(cvkl2))]) %>% round(3)
        
        risk_ols <- NA
        risk_l2brm <- NA
        
        if(tau == 0.5){
          
          set.seed(101)
          risk_ols <- k_fold_cv(tau, method = "ols") %>% round(3)
          
          l2brm <-  mboost(
            log(sell) ~ drv + rec + ffin + ghw + ca + gar + 
              reg + log(lot) + bdms + fb + sty,
            data = dat, 
            baselearner = "bols",
            family = Gaussian(), 
            control = boost_control(mstop = 200, nu = 0.1)
          )
          
          set.seed(101)
          cvkbrm <-
            cvrisk(
              l2brm,
              grid = seq(10, 500, 10),
              folds = mboost::cv(model.weights(l2brm), type = "kfold", B = 5)
            )
          
          risk_l2brm <- mean(cvkbrm[,as.character(mstop(cvkbrm))]) %>% round(3)
          
        }
        
        data.frame(
          tau = tau, 
          RQ = rq_risk, 
          AL1BRQ = risk_al1brq,
          L2BRQ = risk_l2brq, 
          OLS = risk_ols, 
          L2BRM = risk_l2brm
        )
      }) %>% 
      kable(
        format = "latex", 
        caption = "Risk estimated by 5-fold CV"
      )
  }

# Application results ----------------------------------------------------------

coef_table(0.1)
coef_table(0.3)
coef_table(0.5)
coef_table(0.7)
coef_table(0.9)

coef_table(0.75)
coef_table(0.25)

risk_table(seq(0.1, 0.9, 0.2))
highdim(c(0.25, 0.5, 0.75))

