bgmmm <- function(y, K, priors = NULL, options = NULL, parallel = FALSE){
  
  y       <- as.matrix(y)
  priors  <- .get_priors(priors, y)
  options <- .get_options(options)
  
  object          <- NULL
  object$call     <- match.call()
  object$settings <- list(
    K       = K,
    priors  = priors,
    options = options
  )
  
  ### fit the models
  object$models <- .prepare_models(K)
  
  if(parallel == FALSE){
    for(m in 1:length(object$models)){
      object$models[[m]] <- .fit_jags(y, object$models[[m]][["K"]], object$models[[m]][["model_type"]], priors, options)
    }
  }else{
    cl <- parallel::makePSOCKcluster(min(c(parallel::detectCores(), length(object$models))))
    #parallel::clusterEvalQ(cl, {library("bglmmm")})
    parallel::clusterExport(cl, c("y", "object",".fit_jags_parallel", ".fit_jags", ".transform_post", ".get_par_bounds", ".get_marglik_post", ".marglik_func", ".prepare_inits", ".prepare_data", ".get_model", ".marglik_fail"), envir = environment())
    object$models <- parallel::clusterApplyLB(cl, length(object$models):1, .fit_jags_parallel, y = y, object = object, priors = priors, options = options)[length(object$models):1]
    parallel::stopCluster(cl)
  }

  
  ### posterior probabilities and inclusion Bayes factors
  object$inference <- .bglmmm_inference(object$models)

  object$models <- .name_models(object$models)
  
  class(object) <- "bgmmm"
  return(object)
}

#### fitting functions ####
.prepare_models    <- function(K){
  
  models <- list(
    list(K = 1, model_type = "mc")
  )
  
  i <- 1
  for(k in 2:K){
    for(model_type in c("Mc", "mC", "MC")){
      i           <- i + 1
      models[[i]] <- list(K = k, model_type = model_type)
    }
  }
    
  return(models)
}
.name_models       <- function(models){
  names(models) <- sapply(models, function(model)paste0(model$model_type, "_", model$K))
  return(models)
}
.prepare_data      <- function(y, K, priors, model_type){
  
  # get the prior values
  prior_k     <- priors[["prior_k"]]
  prior_r     <- priors[["prior_r"]]
  prior_alpha <- priors[["prior_alpha"]]
  prior_mu    <- priors[["prior_mu"]]
  prior_S     <- priors[["prior_S"]]
  
  # data structure information
  N <- nrow(y)
  R <- ncol(y)
  
  # prior for precision of the Wishard distribution
  prior_R <- diag(prior_r, R)
  
  # prior for the means
  prior_mu <- rep(prior_mu, R)
  prior_S  <- diag(prior_S, R)
  
  # prior for the mixing proportions
  prior_alpha <- rep(prior_alpha, K)
  
  jags_data <- list(
    R = R,
    y = y,
    N = nrow(y),
    
    prior_R     = prior_R,
    prior_k     = prior_k,
    prior_mu    = prior_mu,
    prior_S     = prior_S
  )
  
  if(model_type != "mc"){
    jags_data$K           <- K
    jags_data$prior_alpha <- prior_alpha
  }
  
  return(jags_data)
}
.prepare_inits     <- function(jags_data, model_type, chains, seed = NULL){
  
  K <- if(model_type == "mc") 1 else jags_data$K
  R <- jags_data$R
  y <- jags_data$y
  N <- jags_data$N
  prior_R     <- jags_data$prior_R
  prior_k     <- jags_data$prior_k
  prior_mu    <- jags_data$prior_mu
  prior_S     <- jags_data$prior_S
  prior_alpha <- jags_data$prior_alpha
  
  
  if(is.null(seed)){
    seed <- sample(666666, 1)
  }
  inits <- vector(mode = "list", chains)
  set.seed(seed)
  
  for(i in 1:chains){
    
    temp_init <- list()
    
    # starting values
    mu_init <- array(numeric(K*R), c(K,R))
    for(k in 1:K){
      mu_init[k,] <- MASS::mvrnorm(1, prior_mu, solve(prior_S))
    }
    
    tau_init <- array(numeric(K*R*R), c(K,R,R))
    for(k in 1:K){
      tau_init[k,,] = rWishart(1, prior_k, prior_R)[,,1]
    }
    
    if(model_type != "mc"){
      c_init <- apply(LaplacesDemon::rdirichlet(nrow(y), prior_alpha), 1, function(x)order(order(x))[1])  
    }else{
      c_init <- NULL
    }
    
    # remove starting values for multiple components if they are not present
    if(grepl("m", model_type, ignore.case = FALSE, fixed = TRUE)){
      mu_init  <- mu_init[1,]
    }
    if(grepl("c", model_type, ignore.case = FALSE, fixed = TRUE)){
      tau_init <- tau_init[1,,]
    }
    
    temp_init <- list(
      mu  = mu_init,
      tau = tau_init
    )
    
    if(model_type != "mc"){
      temp_init$c <- c_init
    }
    
    temp_init[[".RNG.seed"]] <- seed + i
    temp_init[[".RNG.name"]] <- "base::Super-Duper"
    
    inits[[i]] <- temp_init
  }
  
  return(inits)
}
.marglik_func      <- function(samples.row, data, model_type){
  
  # extract data
  K <- if(model_type == "mc") 1 else data$K
  N <- data$N
  R <- data$R
  y <- data$y
  
  # extract priors
  prior_R     <- data$prior_R
  prior_k     <- data$prior_k
  prior_mu    <- data$prior_mu
  prior_S     <- data$prior_S
  prior_alpha <- data$prior_alpha
  
  
  ### get parameter estimates depending on the model type
  pi    <- NULL
  mu    <- NULL
  tau   <- NULL
  chol  <- NULL
  
  # mixing proportion parameters & classes
  if(model_type != "mc"){
    pi <- samples.row[ paste0("pi[",1:(K-1),"]") ]
    pi <- c(pi, 1 - sum(pi))
  }
  
  # mean parameters
  if(grepl("M", model_type, ignore.case = FALSE, fixed = TRUE)){
    for(k in 1:K){
      mu[[k]] <- samples.row[ paste0("mu[",k,",",1:R,"]") ]
    }
  }else{
    mu <- samples.row[ paste0("mu[",1:R,"]") ]
  }
  
  # precision matrices - construct back from Cholesky + log transformation
  if(grepl("C", model_type, ignore.case = FALSE, fixed = TRUE)){
    for(k in 1:K){
      chol_names <- sapply(1:R, function(r)paste0("chol[",k,",",1:R,",",r,"]"))
      chol_names <- chol_names[upper.tri(chol_names, diag = TRUE)]
      
      chol[[k]]  <- matrix(0, ncol = R, nrow = R)
      chol[[k]][upper.tri(chol[[k]], diag = TRUE)] <- samples.row[ chol_names ]
      diag(chol[[k]]) <- exp(diag(chol[[k]]))
      tau[[k]]  <- t(chol[[k]]) %*% chol[[k]]
    }
  }else{
    chol_names <- sapply(1:R, function(r)paste0("chol[",1:R,",",r,"]"))
    chol_names <- chol_names[upper.tri(chol_names, diag = TRUE)]
    
    chol  <- matrix(0, ncol = R, nrow = R)
    chol[upper.tri(chol, diag = TRUE)] <- samples.row[ chol_names ]
    diag(chol) <- exp(diag(chol))
    tau  <- t(chol) %*% chol
  }
  
  
  ### compute marginal likelihood
  log_lik <- 0
  
  # adding Jacobean due to the Cholesky transformation used in bridge sampling
  # http://rstudio-pubs-static.s3.amazonaws.com/486816_440106f76c944734a7d4c84761e37388.html
  # https://mc-stan.org/docs/2_21/stan-users-guide/changes-of-variables.html
  # https://mc-stan.org/docs/2_18/reference-manual/covariance-matrices-1.html
  log_Jacobean <- 0
  
  # priors
  if(model_type != "mc"){
    log_lik <- log_lik + LaplacesDemon::ddirichlet(pi, prior_alpha, log = TRUE)
  }
  
  if(grepl("M", model_type, ignore.case = FALSE, fixed = TRUE)){
    for(k in 1:K){
      log_lik <- log_lik + LaplacesDemon::dmvn(mu[[k]], mu = prior_mu, Sigma = LaplacesDemon::Prec2Cov(prior_S), log = TRUE)
    }
  }else{
    log_lik <- log_lik + LaplacesDemon::dmvn(mu, mu = prior_mu, Sigma = LaplacesDemon::Prec2Cov(prior_S), log = TRUE)
  }
  
  if(grepl("C", model_type, ignore.case = FALSE, fixed = TRUE)){
    for(k in 1:K){
      log_lik      <- log_lik      + LaplacesDemon::dwishart(tau[[k]], prior_k, prior_R, log = TRUE)
      log_Jacobean <- log_Jacobean + log(2^R) + sum( sapply(1:R,function(x)log( chol[[k]][x,x]^(R-x+2)) ) )
    }
  }else{
    log_lik      <- log_lik      + LaplacesDemon::dwishart(tau, prior_k, prior_R, log = TRUE)
    log_Jacobean <- log_Jacobean + log(2^R) + sum( sapply(1:R,function(x)log( chol[x,x]^(R-x+2)) ) )
  }
  
  # data
  if(model_type == "mc"){
    log_lik  <- log_lik + sum(LaplacesDemon::dmvn(y, mu = mu, Sigma = LaplacesDemon::Prec2Cov(tau), log = TRUE))
  }else{
    # marginalizing over the indicies https://mc-stan.org/users/documentation/case-studies/identifying_mixture_models.html
    dat_lik  <- sapply(1:K, function(k){
      
      if(grepl("M", model_type, ignore.case = FALSE, fixed = TRUE)){
        temp_mu <- mu[[k]]
      }else{
        temp_mu <- mu
      }
      
      if(grepl("C", model_type, ignore.case = FALSE, fixed = TRUE)){
        temp_tau <- tau[[k]]
      }else{
        temp_tau <- tau
      }
      
      LaplacesDemon::dmvn(y, mu = temp_mu, Sigma = LaplacesDemon::Prec2Cov(temp_tau), log = TRUE) + log(pi[k])
    })
    log_lik <- log_lik + sum(log(rowSums(exp(dat_lik))))
  }
  
  return(log_lik + log_Jacobean)
}
.marglik_fail      <- function(message){
  marg_lik           <- NULL
  marg_lik$logml     <- -Inf
  marg_lik$message   <- message
  class(marg_lik)    <- "bridge"
  return(marg_lik)
}
.get_marglik_post  <- function(posterior, jags_data, model_type){
  
  K <- if(model_type == "mc") 1 else jags_data$K
  R <- jags_data$R
  
  ### remove in indicators
  posterior <- posterior[, !grepl("c",  colnames(posterior))]
  
  ### prepare the precision matrices for bridge-sampling
  # transform precision tau into covariance
  for(k in 1:K){
    
    if(grepl("C", model_type, ignore.case = FALSE, fixed = TRUE)){
      tau_names    <- sapply(1:R, function(r)paste0("tau[",k,",",1:R,",",r,"]"))
    }else{
      tau_names    <- sapply(1:R, function(r)paste0("tau[",1:R,",",r,"]"))
    }
    
    tau_samples  <- posterior[,as.vector(tau_names)]
    
    # Cholesky + log transform
    chol_samples <- t(sapply(1:nrow(tau_samples), function(i){
      temp_chol       <- chol(matrix(tau_samples[i,], R, R))
      diag(temp_chol) <- log(diag(temp_chol))
      temp_chol[upper.tri(temp_chol, diag = TRUE)]
    }))
    colnames(chol_samples) <- gsub("tau", "chol", as.vector(tau_names[upper.tri(tau_names, diag = TRUE)]))
    
    # replace posterior samples by cholesky transform
    posterior <- posterior[,!colnames(posterior) %in% as.vector(tau_names)]
    posterior <- cbind(posterior,chol_samples)
    
    if(grepl("c", model_type, ignore.case = FALSE, fixed = TRUE))break
  }
  
  return(posterior)
}
.get_par_bounds    <- function(jags_data, model_type){
  
  K <- jags_data$K
  N <- jags_data$N
  R <- jags_data$R
  
  ### get parameters based on the model type
  pars        <- NULL
  lb          <- NULL
  ub          <- NULL
  param_types <- NULL
  
  # mixing proportion parameters
  if(model_type != "mc"){
    pars        <- c(pars,        paste0("pi[",1:K,"]"))
    lb          <- c(lb,          rep(0, K))
    ub          <- c(ub,          rep(1, K))
    param_types <- c(param_types, rep("simplex", K))
  }
  
  # mean parameters
  if(grepl("M", model_type, ignore.case = FALSE, fixed = TRUE)){
    for(k in 1:K){
      pars        <- c(pars,        paste0("mu[",k,",",1:R,"]"))
      lb          <- c(lb,          rep(-Inf, R))
      ub          <- c(ub,          rep( Inf, R))
      param_types <- c(param_types, rep("real", R))
    }
  }else{
    pars        <- c(pars,        paste0("mu[",1:R,"]"))
    lb          <- c(lb,          rep(-Inf, R))
    ub          <- c(ub,          rep( Inf, R))
    param_types <- c(param_types, rep("real", R))
  }
  
  # precision matrices
  if(grepl("C", model_type, ignore.case = FALSE, fixed = TRUE)){
    for(k in 1:K){
      
      tau_names   <- sapply(1:R, function(r)paste0("tau[",k,",",1:R,",",r,"]"))
      chol_names  <- gsub("tau", "chol", as.vector(tau_names[upper.tri(tau_names, diag = TRUE)]))
      
      pars        <- c(pars,        chol_names)
      lb          <- c(lb,          rep(-Inf,   length(chol_names)))
      ub          <- c(ub,          rep( Inf,   length(chol_names))) 
      param_types <- c(param_types, rep("real", length(chol_names)))
    }
  }else{
    
    tau_names   <- sapply(1:R, function(r)paste0("tau[",1:R,",",r,"]"))
    chol_names  <- gsub("tau", "chol", as.vector(tau_names[upper.tri(tau_names, diag = TRUE)]))
    
    pars        <- c(pars,        chol_names)
    lb          <- c(lb,          rep(-Inf,   length(chol_names)))
    ub          <- c(ub,          rep( Inf,   length(chol_names))) 
    param_types <- c(param_types, rep("real", length(chol_names)))
  }
  
  names(lb) <- names(ub) <- names(param_types) <- pars
  
  return(list(
    pars        = pars,
    lb          = lb,
    ub          = ub,
    param_types = param_types
  ))
}
.transform_post    <- function(posterior, jags_data, model_type){
  
  # extract data
  K <- if(model_type == "mc") 1 else jags_data$K
  N <- jags_data$N
  R <- jags_data$R
  y <- jags_data$y
  
  output <- list(
    mu    = NULL,
    Sigma = NULL
  )
  
  # mixing proportion parameters & classes
  if(model_type != "mc"){
    output[["pi"]] <- posterior[,paste0("pi[",1:K,"]") ]
    output[["c" ]] <- posterior[,paste0( "c[",1:N,"]") ]
  }
  
  # mean parameters
  if(grepl("M", model_type, ignore.case = FALSE, fixed = TRUE)){
    for(k in 1:K){
      output[["mu"]][[k]] <- posterior[,paste0("mu[",k,",",1:R,"]")]
    }
  }else{
    output[["mu"]] <- posterior[,paste0("mu[",1:R,"]")]
  }
  
  # covariance matrices - construct from precision
  if(grepl("C", model_type, ignore.case = FALSE, fixed = TRUE)){
    for(k in 1:K){
      tau_names              <- sapply(1:R, function(r)paste0("tau[",k,",",1:R,",",r,"]"))
      output[["Sigma"]][[k]] <- array(unlist(sapply(1:nrow(posterior), function(i){
        LaplacesDemon::Prec2Cov(matrix(posterior[i, as.vector(tau_names)], R, R))
      }, simplify = FALSE)), dim = c(R, R, nrow(posterior)))
    }
  }else{
    tau_names         <- sapply(1:R, function(r)paste0("tau[",1:R,",",r,"]"))
    output[["Sigma"]] <- array(unlist(sapply(1:nrow(posterior), function(i){
      LaplacesDemon::Prec2Cov(matrix(posterior[i, as.vector(tau_names)], R, R))
    }, simplify = FALSE)), dim = c(R, R, nrow(posterior)))
  }
  
  return(output)
}
.get_options       <- function(options){
  if(is.null(options)){
    options[["chains"]]      <- 1
    options[["adapt"]]       <- 500
    options[["burnin"]]      <- 2000
    options[["iter"]]        <- 5000
    options[["bridge_iter"]] <- 100000
  }else{
    if(is.null(options[["chains"]])){
      options[["chains"]]      <- 1
    }else{
      options[["chains"]]      <- options[["chains"]]
    }
    if(is.null(options[["adapt"]])){
      options[["adapt"]]       <- 500
    }else{
      options[["adapt"]]       <- options[["adapt"]]
    } 
    if(is.null(options[["burnin"]])){
      options[["burnin"]]      <- 2000
    }else{
      options[["burnin"]]      <- options[["burnin"]]
    }
    if(is.null(options[["iter"]])){
      options[["iter"]]        <- 5000
    }else{
      options[["iter"]]        <- options[["iter"]]
    }
    if(is.null(options[["bridge_iter"]])){
      options[["bridge_iter"]] <- 100000
    }else{
      options[["bridge_iter"]] <- options[["bridge_iter"]]
    }
  }
  
  return(options)
}
.get_priors        <- function(priors, y){
  # set default priors if they are unspecified
  if(is.null(priors)){
    priors[["prior_k"]]     <- 2*ncol(y)
    priors[["prior_r"]]     <- 1/ncol(y)
    priors[["prior_alpha"]] <- nrow(y) * 0.05
    priors[["prior_mu"]]    <- 0
    priors[["prior_S"]]     <- 1
  }else{
    if(is.null(priors[["prior_k"]])){
      priors[["prior_k"]]     <- 2*ncol(y)
    }else{
      priors[["prior_k"]]     <- priors[["prior_k"]] 
    }
    if(is.null(priors[["prior_r"]])){
      priors[["prior_r"]]     <- 1/ncol(y)
    }else{
      priors[["prior_r"]]     <- priors[["prior_r"]] 
    }
    if(is.null(priors[["prior_alpha"]])){
      priors[["prior_alpha"]] <- nrow(y) * 0.05
    }else{
      priors[["prior_alpha"]] <- priors[["prior_alpha"]] 
    }
    if(is.null(priors[["prior_mu"]])){
      priors[["prior_mu"]]    <- 0
    }else{
      priors[["prior_mu"]]    <- priors[["prior_mu"]] 
    }
    if(is.null(priors[["prior_S"]])){
      priors[["prior_S"]]     <- 1
    }else{
      priors[["prior_S"]]     <- priors[["prior_S"]] 
    }
  }
  
  return(priors)
}
.fit_jags          <- function(y, K, model_type, priors, options){
  
  fit_data  <- .prepare_data(y, K, priors, model_type)
  fit_inits <- .prepare_inits(fit_data, model_type, options[["chains"]])
  
  
  # fit the model
  # fit the model
  fit_jags  <- rjags::jags.model(
    file     = textConnection(.get_model(model_type)), 
    data     = fit_data,
    inits    = fit_inits,
    n.chains = options[["chains"]],
    n.adapt  = options[["adapt"]])
  update(fit_jags, options[["burnin"]])
  posterior <- rjags::coda.samples(fit_jags, c("mu", "tau", if(model_type != "mc")c("pi", "c")), options[["iter"]])
  posterior <- do.call(rbind, posterior)
  posterior <- coda::as.mcmc(posterior)
  
  # fit_jags <- runjags::run.jags(
  #   model           = .get_model(model_type),
  #   data            = fit_data,
  #   inits           = fit_inits,
  #   adapt           = options[["adapt"]],
  #   burnin          = options[["burnin"]], 
  #   sample          = options[["iter"]],
  #   n.chains        = options[["chains"]],
  #   monitor         = c("mu", "tau", if(model_type != "mc")c("pi", "c")),
  #   summarise       = FALSE
  # )
  # posterior <- coda::as.mcmc(fit_jags)
  
  
  # compute marginal likelihood
  posterior_bounds  <- .get_par_bounds(fit_data, model_type)
  marglik_post      <- .get_marglik_post(posterior, fit_data, model_type)
  marglik_post      <- marglik_post[, posterior_bounds$pars]
  
  marg_lik        <- tryCatch(suppressWarnings(bridgesampling::bridge_sampler(
    samples          = marglik_post,
    data             = fit_data,
    model_type       = model_type,
    log_posterior    = .marglik_func,
    lb               = posterior_bounds$lb,
    ub               = posterior_bounds$ub,
    param_types      = posterior_bounds$param_types,
    maxiter          = options[["bridge_iter"]],
    silent           = TRUE
  )),error = function(e)return(e))
  
  if(any(class(marg_lik) %in% c("simpleError", "error"))){
    marg_lik <- .marglik_fail(marg_lik$message)
  }
  
  if(is.na(marg_lik$logml)){
    marg_lik <- .marglik_fail("marg_lik NA")
  }
  
  
  # create the model object
  model <- list(
    K          = K,
    model_type = model_type,
    posterior  = .transform_post(posterior, fit_data, model_type),
    marg_lik   = marg_lik
  )
  class(model) <- "bgmmm.model"
  
  return(model)
}
.fit_jags_parallel <- function(x, y, object, priors = NULL, options){
  .fit_jags(y, object$models[[x]][["K"]], object$models[[x]][["model_type"]], priors, options)
}

#### Bayes factor functions ####
.bglmmm_inference <- function(models){
  
  K <- length(models)
  
  marg_liks   <- sapply(models, function(model)model$marg_lik$logml)
  model_Ks    <- sapply(models, function(model)model$K)
  model_types <- sapply(models, function(model)model$model_type)
  names(marg_liks) <- paste0(model_types, "_", model_Ks)
  
  # number of clusters
  prior_odds_Ks <- sapply(model_Ks, function(model_K)1/sum(model_Ks == model_K))
  prior_prob_Ks <- prior_odds_Ks/sum(prior_odds_Ks)
  post_prob_Ks  <- bridgesampling::post_prob(marg_liks, prior_prob = prior_prob_Ks)
  inclusion_K   <- sapply(unique(model_Ks), function(model_K).inclusion_BF(prior_prob_Ks, post_prob_Ks, model_Ks == model_K))
  post_prob_K   <- sapply(unique(model_Ks), function(model_K)sum(post_prob_Ks[model_Ks == model_K]))
  names(inclusion_K) <- names(post_prob_K) <- unique(model_Ks)
  
  # model types
  prior_odds_types <- sapply(model_types, function(model_type)1/sum(model_types == model_type))
  prior_prob_types <- prior_odds_types/sum(prior_odds_types)
  post_prob_types  <- bridgesampling::post_prob(marg_liks, prior_prob = prior_prob_types)
  inclusion_type   <- sapply(unique(model_types), function(model_type).inclusion_BF(prior_prob_types, post_prob_types, model_types == model_type))
  post_prob_type   <- sapply(unique(model_types), function(model_type)sum(post_prob_types[model_types == model_type]))
  names(inclusion_type) <- names(post_prob_type) <- unique(model_types)
  
  # conditional 
  prior_prob_conditional <- rep(1/(length(marg_liks[-1])), length(marg_liks[-1]))
  post_prob_conditional   <- bridgesampling::post_prob(marg_liks[-1], prior_prob = prior_prob_conditional)
  inclusion_conditional   <- sapply(1:length(marg_liks[-1]), function(i).inclusion_BF(prior_prob_conditional, post_prob_conditional, 1:length(marg_liks[-1]) == i))
  names(post_prob_conditional)   <- names(inclusion_conditional) <- names(models)[-1]
  
  prior_odds_structures <- c(1, prior_prob_conditional)
  prior_prob_structures <- prior_odds_structures/sum(prior_odds_structures)
  post_prob_structures  <- bridgesampling::post_prob(marg_liks, prior_prob = prior_prob_structures)
  inclusion_structure   <- c(
    .inclusion_BF(prior_prob_structures, post_prob_structures, 1 == 1:length(marg_liks)),
    .inclusion_BF(prior_prob_structures, post_prob_structures, 1 != 1:length(marg_liks)))
  post_prob_structure   <- c(post_prob_structures[1], sum(post_prob_structures[-1]))
  names(inclusion_structure) <- names(post_prob_structure) <- c("No", "Yes")
  
  return(list(
    models_marglik = marg_liks,
    post_prob      = list(
      K               = post_prob_K,
      model_types     = post_prob_type,
      structure       = post_prob_structure,
      conditional     = post_prob_conditional
    ),
    BF             = list(
      K               = inclusion_K,
      model_types     = inclusion_type,
      structure       = inclusion_structure,
      conditional     = inclusion_conditional
    )
  ))
}
.inclusion_BF     <- function(prior_weights, posterior_weights, conditional_models){
  (sum(posterior_weights[conditional_models])/sum(posterior_weights[!conditional_models]))  /
    (sum(prior_weights[conditional_models])/sum(prior_weights[!conditional_models]))
}

#### printing functions ####
print.bgmmm       <- function(x, ...){
  
  # components
  overview_K <- cbind(
      names(x$inference$BF$K),
      format(round(1/length(x$inference$BF$K), 3), nsmall = 3),
      format(round(x$inference$post_prob$K, 3),    nsmall = 3),
      format(round(x$inference$BF$K, 3),           nsmall = 3)
  )
  colnames(overview_K) <- c("Components", "Prior Prob. ", "Post. Prob.", "Incl. BF")
  rownames(overview_K) <- rep("", nrow(overview_K))
  
  # structure
  overview_types <- cbind(
      ifelse(grepl("m", names(x$inference$BF$model_types), ignore.case = FALSE), "equivalent", "different"),
      ifelse(grepl("c", names(x$inference$BF$model_types), ignore.case = FALSE), "equivalent", "different"),
      format(round(1/length(x$inference$BF$model_types), 3), nsmall = 3),
      format(round(x$inference$post_prob$model_types, 3),    nsmall = 3),
      format(round(x$inference$BF$model_types, 3),           nsmall = 3)
  )
  colnames(overview_types) <- c("Means", "Covariances", "Prior Prob. ", "Post. Prob.", "Incl. BF")
  rownames(overview_types) <- rep("", nrow(overview_types))
  
  # conditional
  overview_conditional <- matrix(format(round(x$inference$post_prob$conditional, 3), nsmall = 3), nrow = 3)
  colnames(overview_conditional) <- 2:x$settings$K
  rownames(overview_conditional) <- c(" different | equivalent", "equivalent | different", " different | different")
  
  
  cat("Call:\n")
  print(x$call)
  cat("\n")
  
  cat("Bayesian Gaussian Multivariate Mixture Models\n")
  cat("\n")
  
  cat("Test for number of components:\n")
  print(overview_K, quote = FALSE, right = TRUE)
  cat("\n")

  cat("Test for model structure:\n")  
  print(overview_types, quote = FALSE, right = TRUE)
  cat("\n")
  
  cat("Condtional posterior model probabilities (means|covariances \\ components):\n")  
  print(overview_conditional, quote = FALSE, right = TRUE)

}
print.bgmmm.model <- function(x, ...){
  
  cat("Bayesian Gaussian Multivariate Mixture Model\n")
  cat(paste0(
    if(x$K == 1) "1 component\n" else paste0(x$K, " components, ", ifelse(grepl("m", x$model_type, ignore.case = FALSE), "equivalent", "different"), " means, ", ifelse(grepl("c", x$model_type, ignore.case = FALSE), "equivalent", "different"), " covariances\n")
  ))
  cat("\n")
  
  if(x$K != 1){
    cat("mixing proportions (pi):\n")
    print(format(round(apply(x$posterior$pi, 2, mean), 3), nsmall = 3), quote = FALSE)
  }
  cat("\n")
  
  cat("means (mu):\n")
  if(grepl("m", x$model_type, ignore.case = FALSE)){
    print(format(round(apply(x$posterior$mu, 2, mean), 3), nsmall = 3), quote = FALSE)
  }else{
    mu <- sapply(x$posterior$mu, function(mu)apply(mu, 2, mean))
    colnames(mu) <- 1:x$K
    rownames(mu) <- paste0("mu[",1:nrow(mu),"]")
    print(format(round(mu, 3), nsmall = 3), quote = FALSE, right = TRUE)
  }
  cat("\n")
  
  cat("correlations\\covariances (Sigma):\n")
  if(grepl("c", x$model_type, ignore.case = FALSE)){
    Sigma <- .CorCovMat(x$posterior$Sigma)
    Sigma <- apply(Sigma, c(1,2), mean)
    temp_Sigma <- format(round(Sigma, 3), nsmall = 3)
    if(all(temp_Sigma >= 0))temp_Sigma <- matrix(paste0(" ", temp_Sigma), ncol = ncol(temp_Sigma), nrow = nrow(temp_Sigma))
    print(temp_Sigma, quote = FALSE, right = TRUE)
  }else{
    Sigma <- sapply(x$posterior$Sigma, .CorCovMat, simplify = FALSE)
    Sigma <- sapply(Sigma, function(S)apply(S, c(1,2), mean), simplify = F)
    for(k in 1:x$K){
      cat(paste0("Sigma[",k,"]\n"))
      temp_Sigma <- format(round(Sigma[[k]], 3), nsmall = 3)
      if(all(temp_Sigma >= 0))temp_Sigma <- matrix(paste0(" ", temp_Sigma), ncol = ncol(temp_Sigma), nrow = nrow(temp_Sigma))
      print(temp_Sigma, quote = FALSE, right = TRUE)
      cat("\n")
    }
  }
}

.CorCovMat <- function(x){
  for(i in 1:dim(x)[3]){
    x[,,i][lower.tri(x[,,i])] <- cov2cor(x[,,i])[lower.tri(x[,,i])]
  }
  return(x)
} 

#### models ####
.get_model <- function(model_type){
  # based on: https://sourceforge.net/p/mcmc-jags/discussion/610037/thread/54347058/
  
  if(model_type == "MC"){
    return('
      model{
  
      # Each data sample picks a cluster
      pi ~ ddirch(prior_alpha)
  
      for(i in 1:N) {
        c[i]  ~ dcat(pi)
      }
  
      # Each cluster has a gaussian
      for(k in 1:K) {
        mu[k,1:R]        ~ dmnorm(prior_mu[],prior_S[,])
        tau[k,1:R,1:R]   ~ dwish(prior_R[,], prior_k)
      }
  
      # Each data sample is distributed from one of the clusters
      for(i in 1:N) {
        y[i,1:R] ~ dmnorm(mu[c[i],], tau[c[i],,])
      }
   
      }')
  }else if(model_type == "Mc"){
    return(
      'model{

      # Each data sample picks a cluster
      pi ~ ddirch(prior_alpha)
  
      for(i in 1:N) {
        c[i]  ~ dcat(pi)
      }
  
      # Common covariance matrix
      tau   ~ dwish(prior_R[,], prior_k)
  
      # Each cluster has a gaussian
      for(k in 1:K) {
        mu[k,1:R]        ~ dmnorm(prior_mu[],prior_S[,])
      }
  
      # Each data sample is distributed from one of the clusters
      for(i in 1:N) {
        y[i,1:R] ~ dmnorm(mu[c[i],],tau[,])
      }
  
      }')
  }else if(model_type == "mC"){
    return(
      'model{

      # Each data sample picks a cluster
      pi ~ ddirch(prior_alpha)
  
      for(i in 1:N) {
        c[i]  ~ dcat(pi)
      }
  
      # Common mean
      mu ~ dmnorm(prior_mu, prior_S[,])
  
      # Each cluster has a gaussian
      for(k in 1:K) {
        tau[k,1:R,1:R]   ~ dwish(prior_R[,], prior_k)
      }
  
      # Each data sample is distributed from one of the clusters
      for(i in 1:N) {
        y[i,1:R] ~ dmnorm(mu[],tau[c[i],,])
      }
  
      }')
  }else if(model_type == "mc"){
    return(
      'model{

      # Common mean and covariance
      mu ~ dmnorm(prior_mu, prior_S[,])
      tau[1:R,1:R]   ~ dwish(prior_R[,], prior_k)
  
      # Each data sample is distributed from one of the clusters
      for(i in 1:N) {
        y[i,1:R] ~ dmnorm(mu[],tau[,])
      }
  
      }')
  }
}
