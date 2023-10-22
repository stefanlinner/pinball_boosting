# Parameters for different cases -----------------------------------------------

sim <- 
  list(
    # homoskedastic
    hom = list(
      zbrq = 1:200*20,
      zqb = 1:250*100, 
      beta = c(3, 1), 
      alpha = c(4, 0),
      n = 200,
      p = 2,
      data_gen = function(seed = NULL, n, error, tau, contaminated) {
        beta <- c(3, 1)
        alpha <- c(4, 0)
        
        if(!is.null(seed)){
          set.seed(seed)
        }
        x1 <- runif(n = n, min = 0, max = 10)
        
        if(error == "norm"){
          u  <- rnorm(n)
        } else if(error == "tdist"){
          u  <- rt(n, df = 2)
        } else if(error == "gamma"){
          u  <- rgamma(n, shape = 2, scale = 1)
        } else if(error == "mixed"){
          u <- numeric(n)
          u[x1<median(x1)] <- rt(sum(x1<median(x1)), df = 2)
          u[x1>=median(x1)] <- 
            rt(sum(x1>=median(x1)), df = 2) + rgamma(sum(x1>=median(x1)), shape = 2, scale = 1)
        } else if(error == "tdist_1"){
          u  <- rt(n, df = 1)
        } else {
          stop("invalid error")
        }
        
        X  <- matrix(c(rep(1, n), x1), nrow = n, dimnames = list(NULL, c("Int", "x1")))
        y  <- X %*% beta + as.vector((X %*% alpha)) * u
        
        if(contaminated){
          y[which.max(y)] <- max(y) * 10
        }
        
        dat <- 
          cbind(y = y, X[, -1]) %>% 
          as.data.frame() %>%
          set_names(c("y", colnames(X)[-1]))
        
        dat
      },
      offset = function(exact.fit, y, init){
        ifelse(
          exact.fit, 
          quantile(y, probs = init, type = 1), 
          quantile(y, probs = init)
        )
      },
      formulabrq = function() as.formula(y ~ brq(x1)),
      formulaqb = function() as.formula(y ~ bols(x1))
    ), 
    
    # heteroskedastic
    het = list(
      zbrq = 1:250*20, 
      zqb = 1:400*100, 
      beta = c(4, 2), 
      alpha = c(4, 1),
      n = 200,
      p = 2,
      data_gen = function(seed = NULL, n, error, tau, contaminated) {
        beta <- c(4, 2)
        alpha <- c(4, 1)
        
        if(!is.null(seed)){
          set.seed(seed)
        }
        x1 <- runif(n = n, min = 0, max = 10)
        
        if(error == "norm"){
          u  <- rnorm(n)
        } else if(error == "tdist"){
          u  <- rt(n, df = 2)
        } else if(error == "gamma"){
          u  <- rgamma(n, shape = 2, scale = 1)
        } else if(error == "mixed"){
          u <- numeric(n)
          u[x1<median(x1)] <- rt(sum(x1<median(x1)), df = 2)
          u[x1>=median(x1)] <- 
            rt(sum(x1>=median(x1)), df = 2) + rgamma(sum(x1>=median(x1)), shape = 2, scale = 1)
        } else if(error == "tdist_1") {
          u  <- rt(n, df = 1)
        } else {
          stop("invalid error")
        }
        X  <- matrix(c(rep(1, n), x1), nrow = n, dimnames = list(NULL, c("Int", "x1")))
        y  <- X %*% beta + as.vector((X %*% alpha)) * u
        
        if(contaminated){
          y[which.max(y)] <- max(y) * 10
        }
        
        dat <- 
          cbind(y = y, X[, -1]) %>% 
          as.data.frame() %>% 
          set_names(c("y", colnames(X)[-1]))
        
        dat
      },
      offset = function(exact.fit, y, init){
        ifelse(
          exact.fit, 
          quantile(y, probs = init, type = 1), 
          quantile(y, probs = init)
        )
      },
      formulabrq = function() as.formula(y ~ brq(x1)), 
      formulaqb = function() as.formula(y ~ bols(x1))
    ),
    
    # multivariate
    multi = list(
      zbrq = 1:200*100, 
      zqb = 1:250*250, 
      beta = c(5, 8, -5, 2, -2, 0, 0), 
      alpha = c(1, 0, 2, 0, 1, 0, 0),
      n = 500,
      p = 7,
      data_gen = function(seed = NULL, n, error, tau, contaminated) {
        
        
        beta <- c(5, 8, -5, 2, -2, 0, 0)
        alpha <- c(1, 0, 2, 0, 1, 0, 0)
        
        nam <- map_chr(seq_len(6), function(.x){
          paste("x", .x, sep = "")
        })
        
        if(!is.null(seed)){
          set.seed(seed)
        }
        X <- cbind(
          Int = rep(1, n),
          map_dfc(seq_len(6), function(.x){
            x.dat <- data.frame(runif(n = n, min = 0, max = 10))
            names(x.dat) <- paste("x_", .x, sep = "")
            x.dat
          }) %>% 
            as.data.frame() %>% 
            set_names(nam) %>% 
            as.matrix()
        )
        
        if(error == "norm"){
          u  <- rnorm(n)
        } else if(error == "tdist"){
          u  <- rt(n, df = 2)
        } else if(error == "gamma"){
          u  <- rgamma(n, shape = 2, scale = 1)
        } else if(error == "mixed"){
          u <- numeric(n)
          u[X[,"x1"]<median(X[,"x1"])] <- rt(sum(X[,"x1"]<median(X[,"x1"])), df = 2)
          u[X[,"x1"]>=median(X[,"x1"])] <- 
            rt(sum(X[,"x1"]>=median(X[,"x1"])), df = 2) + rgamma(sum(X[,"x1"]>=median(X[,"x1"])), shape = 2, scale = 1)
        } else if(error == "tdist_1") {
          u  <- rt(n, df = 1)
        } else {
          stop("invalid error")
        }
        
        y  <- X %*% beta + as.vector((X %*% alpha)) * u
        
        if(contaminated){
          y[which.max(y)] <- max(y) * 10
        }
        
        dat <- 
          cbind(y = y, X[, -1]) %>%
          as.data.frame() %>%
          set_names(c("y", colnames(X)[-1]))
        
        dat
      },
      offset = function(exact.fit, y, init){
        ifelse(
          exact.fit, 
          quantile(y, probs = init, type = 1), 
          quantile(y, probs = init)
        )
      },
      formulabrq = function(){
        
        covars <- map_chr(seq_len(6), function(.x){
          paste("x", .x, sep = "")
        })
        
        map_chr(covars, function(.x){
          paste("brq(", .x, ")", sep = "")
        }) %>% 
          paste(collapse = " + ") %>% 
          paste("y ~ ", ., sep = "") %>% 
          as.formula()
      }, 
      formulaqb = function(){
        
        covars <- map_chr(seq_len(6), function(.x){
          paste("x", .x, sep = "")
        })
        
        map_chr(covars, function(.x){
          paste("bols(", .x, ")", sep = "")
        }) %>% 
          paste(collapse = " + ") %>% 
          paste("y ~ ", ., sep = "") %>% 
          as.formula()
      }
    ), 
    
    # highdimensional 
    highdim = list(
      zbrq = 1:150*150, 
      zqb = 1:200*250, 
      beta = c(5, 8, -5, 2, -2, rep(0, 95)), 
      alpha = c(1, 0, 2, 0, 1, rep(0, 95)),
      n = 200,
      p = 100,
      data_gen = function(seed = NULL, n, error, tau, contaminated){
        
        beta <- c(5, 8, -5, 2, -2, rep(0, 95))
        alpha <- c(1, 0, 2, 0, 1, rep(0, 95))
        
        nam <- map_chr(seq_len(99), function(.x){
          paste("x", .x, sep = "")
        })
        
        if(!is.null(seed)){
          set.seed(seed)
        }
        X <- cbind(
          Int = rep(1, n),
          map_dfc(seq_len(99), function(.x){
            x.dat <- data.frame(runif(n = n, min = 0, max = 10))
            names(x.dat) <- paste("x_", .x, sep = "")
            x.dat
          }) %>% 
            as.data.frame() %>% 
            set_names(nam) %>% 
            as.matrix()
        )
        
        if(error == "norm"){
          u  <- rnorm(n)
        } else if(error == "tdist"){
          u  <- rt(n, df = 2)
        } else if(error == "gamma"){
          u  <- rgamma(n, shape = 2, scale = 1)
        } else if(error == "mixed"){
          u <- numeric(n)
          u[X[,"x1"]<median(X[,"x1"])] <- rt(sum(X[,"x1"]<median(X[,"x1"])), df = 2)
          u[X[,"x1"]>=median(X[,"x1"])] <- 
            rt(sum(X[,"x1"]>=median(X[,"x1"])), df = 2) + rgamma(sum(X[,"x1"]>=median(X[,"x1"])), shape = 2, scale = 1)
        } else if(error == "tdist_1") {
          u  <- rt(n, df = 1)
        } else {
          stop("invalid error")
        }
        
        y  <- X %*% beta + as.vector((X %*% alpha)) * u
        
        if(contaminated){
          y[which.max(y)] <- max(y) * 10
        }
        
        dat <- 
          cbind(y = y, X[, -1]) %>% 
          as.data.frame() %>% 
          set_names(c("y", colnames(X)[-1]))
        dat
      },
      offset = function(exact.fit, y, init){
        ifelse(
          exact.fit, 
          quantile(y, probs = init, type = 1), 
          quantile(y, probs = init)
        )
      },
      formulabrq = function(){
        
        covars <- map_chr(seq_len(99), function(.x){
          paste("x", .x, sep = "")
        })
        
        map_chr(covars, function(.x){
          paste("brq(", .x, ")", sep = "")
        }) %>% 
          paste(collapse = " + ") %>% 
          paste("y ~ ", ., sep = "") %>% 
          as.formula()
      }, 
      formulaqb = function(){
        
        covars <- map_chr(seq_len(99), function(.x){
          paste("x", .x, sep = "")
        })
        
        map_chr(covars, function(.x){
          paste("bols(", .x, ")", sep = "")
        }) %>% 
          paste(collapse = " + ") %>% 
          paste("y ~ ", ., sep = "") %>% 
          as.formula()
      }
    ),
    
    # multivariate2
    multi2 = list(
      zbrq = 1:150*100, 
      zqb = 1:200*250, 
      beta = matrix(
        c(
          rep(c(5, 8, -5, 0, -2, 0, 0), 2), 
          c(5, 0, -5, 0, -2, 0, 0), 
          rep(c(5, 0, -5, 2, -2, 0, 0),2)
        ), 
        byrow = TRUE, nrow = 5
      ), 
      alpha = c(1, 0, 2, 0, 1, 0, 0),
      n = 500,
      p = 7,
      data_gen = function(seed = NULL, n, error, tau, contaminated) {
        
        tau.vec <- c(0.1, 0.3, 0.5, 0.7, 0.9)
        beta.mat <-  matrix(
          c(
            rep(c(5, 8, -5, 0, -2, 0, 0), 2), 
            c(5, 0, -5, 0, -2, 0, 0), 
            rep(c(5, 0, -5, 2, -2, 0, 0),2)
          ), 
          byrow = TRUE, nrow = 5
        )
        beta <- beta.mat[which(tau == tau.vec),]
        alpha <- c(1, 0, 2, 0, 1, 0, 0)
        
        nam <- map_chr(seq_len(6), function(.x){
          paste("x", .x, sep = "")
        })
        
        if(!is.null(seed)){
          set.seed(seed)
        }
        X <- cbind(
          Int = rep(1, n),
          map_dfc(seq_len(6), function(.x){
            x.dat <- data.frame(runif(n = n, min = 0, max = 10))
            names(x.dat) <- paste("x_", .x, sep = "")
            x.dat
          }) %>% 
            as.data.frame() %>% 
            set_names(nam) %>% 
            as.matrix()
        )
        
        if(error == "norm"){
          u  <- rnorm(n)
        } else if(error == "tdist"){
          u  <- rt(n, df = 2)
        } else if(error == "gamma"){
          u  <- rgamma(n, shape = 2, scale = 1)
        } else if(error == "mixed"){
          u <- numeric(n)
          u[X[,"x1"]<median(X[,"x1"])] <- rt(sum(X[,"x1"]<median(X[,"x1"])), df = 2)
          u[X[,"x1"]>=median(X[,"x1"])] <- 
            rt(sum(X[,"x1"]>=median(X[,"x1"])), df = 2) + rgamma(sum(X[,"x1"]>=median(X[,"x1"])), shape = 2, scale = 1)
        } else if(error == "tdist_1") {
          u  <- rt(n, df = 1)
        } else {
          stop("invalid error")
        }
        
        y  <- X %*% beta + as.vector((X %*% alpha)) * u
        
        if(contaminated){
          y[which.max(y)] <- max(y) * 10
        }
        
        dat <- 
          cbind(y = y, X[, -1]) %>% 
          as.data.frame() %>% 
          set_names(c("y", colnames(X)[-1]))
        
        dat
      },
      offset = function(exact.fit, y, init){
        ifelse(
          exact.fit, 
          quantile(y, probs = init, type = 1), 
          quantile(y, probs = init)
        )
      },
      formulabrq = function(){
        
        covars <- map_chr(seq_len(6), function(.x){
          paste("x", .x, sep = "")
        })
        
        map_chr(covars, function(.x){
          paste("brq(", .x, ")", sep = "")
        }) %>% 
          paste(collapse = " + ") %>% 
          paste("y ~ ", ., sep = "") %>% 
          as.formula()
      }, 
      formulaqb = function(){
        
        covars <- map_chr(seq_len(6), function(.x){
          paste("x", .x, sep = "")
        })
        
        map_chr(covars, function(.x){
          paste("bols(", .x, ")", sep = "")
        }) %>% 
          paste(collapse = " + ") %>% 
          paste("y ~ ", ., sep = "") %>% 
          as.formula()
      }
    )
  )
