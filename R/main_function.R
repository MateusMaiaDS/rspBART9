# source("R/debugging_rspBART.R")
# rm(list=ls())
# source("R/other_functions.R")
# source("R/sim_functions.R")
# source("R/tree_functions.R")
# source("inst/debugging_rspBART.R")
# devtools::load_all()
set.seed(42)

# Creating the main function from the rspBART
rspBART <- function(x_train,
                    y_train,
                    x_test,
                    n_tree = 10,
                    node_min_size = 15,
                    n_mcmc = 2000,
                    n_burn = 500,
                    alpha = 0.95,
                    beta = 2,
                    df = 3,
                    sigquant = 0.9,
                    kappa = 2,
                    # Splines parameters
                    nIknots = 3,
                    dif_order = 1,
                    tau = 100,
                    scale_bool = TRUE,
                    stump = FALSE,
                    numcut = 100L, # Defining the grid of split rules
                    usequants = FALSE,
                    motrbart_bool = FALSE,
                    use_bs = FALSE,
                    plot_preview = FALSE,
                    all_var = FALSE,
                    scale_init = FALSE,
                    update_tau_beta = FALSE
) {

  # Verifying if x_train and x_test are matrices
  if(!is.data.frame(x_train) || !is.data.frame(x_test)){
    stop("Insert valid data.frame for both data and xnew.")
  }


  if(all_var){
    stop("Wrong version. To use all variables call rspBART8")
  }

  # Getting the valid
  dummy_x <- base_dummyVars(x_train)

  # Create a list
  if(length(dummy_x$facVars)!=0){
    for(i in 1:length(dummy_x$facVars)){
      # See if the levels of the test and train matches
      if(!all(levels(x_train[[dummy_x$facVars[i]]])==levels(x_test[[dummy_x$facVars[i]]]))){
        levels(x_test[[dummy_x$facVars[[i]]]]) <- levels(x_train[[dummy_x$facVars[[i]]]])
      }
      df_aux <- data.frame( x = x_train[,dummy_x$facVars[i]],y)
      formula_aux <- stats::aggregate(y~x,df_aux,mean)
      formula_aux$y <- rank(formula_aux$y)
      x_train[[dummy_x$facVars[i]]] <- as.numeric(factor(x_train[[dummy_x$facVars[[i]]]], labels = c(formula_aux$y)))-1

      # Doing the same for the test set
      x_test[[dummy_x$facVars[i]]] <- as.numeric(factor(x_test[[dummy_x$facVars[[i]]]], labels = c(formula_aux$y)))-1

    }
  }

  # Getting the train and test set
  x_train_scale <- as.matrix(x_train)
  x_test_scale <- as.matrix(x_test)

  # Scaling x
  x_min <- apply(as.matrix(x_train_scale),2,min)
  x_max <- apply(as.matrix(x_train_scale),2,max)

  # Storing the original
  x_train_original <- x_train
  x_test_original <- x_test


  # Normalising all the columns
  for(i in 1:ncol(x_train)){
    x_train_scale[,i] <- normalize_covariates_bart(y = x_train_scale[,i],a = x_min[i], b = x_max[i])
    x_test_scale[,i] <- normalize_covariates_bart(y = x_test_scale[,i],a = x_min[i], b = x_max[i])
  }



  # Creating the numcuts matrix of splitting rules
  xcut_m <- matrix(NA,nrow = numcut,ncol = ncol(x_train_scale))
  for(i in 1:ncol(x_train_scale)){

    if(nrow(x_train_scale)<numcut){
      xcut_m[,i] <- sort(x_train_scale[,i])
    } else {
      xcut_m[,i] <- seq(min(x_train_scale[,i]),
                        max(x_train_scale[,i]),
                        length.out = numcut+2)[-c(1,numcut+2)]
    }
  }




  # =========================================================================================================
  # Getting the Splines Basis functions
  # =========================================================================================================

  # Setting new parameters for the spline
  ndx <- nIknots+1
  ord_ <- 4
  degree_ <- 3
  x_min_sp <- apply(x_train_scale,2,min)
  x_max_sp <- apply(x_train_scale,2,max)
  dx <- (x_max_sp-x_min_sp)/ndx

  # New_knots
  new_knots <- matrix()
  new_knots <- matrix(mapply(x_min_sp,x_max_sp,dx, FUN = function(MIN,MAX,DX){seq(from = MIN-(ord_-1)*DX, to = MAX+(ord_-1)*DX, by = DX)}), ncol = length(dummy_x$continuousVars)) # MIN and MAX are 0 and 1 respectively, because of the scale
  colnames(new_knots) <- dummy_x$continuousVars

  # Selecting which one gonna be used bs or splines.des()
  if(use_bs){
    D_train <- matrix(NA,
                      nrow = nrow(x_train_scale),
                      ncol = (nIknots+degree_)*length(dummy_x$continuousVars))
    D_test <- matrix(NA,
                     nrow = nrow(x_test_scale),
                     ncol = (nIknots+degree_)*length(dummy_x$continuousVars))

  } else {
    D_train <- matrix(NA,
                      nrow = nrow(x_train_scale),
                      ncol = (nrow(new_knots)-ord_)*length(dummy_x$continuousVars))

    D_test <- matrix(NA,
                     nrow = nrow(x_test_scale),
                     ncol = (nrow(new_knots)-ord_)*length(dummy_x$continuousVars))
  }

  # Selecting the basis size.
  if(use_bs){
    basis_size <- (nIknots+degree_)
  } else {
    basis_size <- (nrow(new_knots)-ord_) # Change this value to the desired size of each sublist
  }

  D_seq <- 1:ncol(D_train)  # Replace this with the columns of D

  # Creating a vector
  basis_subindex <- split(D_seq, rep(1:((length(D_seq) %/% basis_size)), each = basis_size, length.out = length(D_seq)))

  # Creating the natural B-spline for each predictor
  for(i in 1:length(basis_subindex)){

    if(use_bs){
      B_train_obj <- splines::bs(x = x_train_scale[,dummy_x$continuousVars[i], drop = FALSE],
                                 df = nIknots+3,
                                 degree = 3,intercept = FALSE,
                                 Boundary.knots = c(-2,2)*range(x_train_scale[,dummy_x$continuousVars[i]]),warn.outside = TRUE)
    } else {
      B_train_obj <- splines::spline.des(x = x_train_scale[,dummy_x$continuousVars[i], drop = FALSE],
                                         knots = new_knots[,dummy_x$continuousVars[i]],
                                         ord = ord_,
                                         derivs = 0*x_train_scale[,dummy_x$continuousVars[i], drop = FALSE],outer.ok = TRUE)$design
    }

    # Returning to MOTR-BART
    if(length(basis_subindex[[i]])!= ncol(B_train_obj)){
      stop("Error on the basis generation")
    }

    D_train[,basis_subindex[[i]]] <- as.matrix(B_train_obj)


    # For the test setting
    if(use_bs){
      B_test_obj <- predict(object = B_train_obj,newx = x_test_scale[,dummy_x$continuousVars[i], drop = FALSE])
    } else {
      B_test_obj <- splines::spline.des(x = x_test_scale[,dummy_x$continuousVars[i], drop = FALSE],
                                        knots = new_knots[,dummy_x$continuousVars[i]],
                                        ord = ord_,
                                        derivs = 0*x_test_scale[,dummy_x$continuousVars[i], drop = FALSE],outer.ok = TRUE)$design
    }

    # Returning to MOTR-BART
    if(length(basis_subindex[[i]])!= ncol(B_test_obj)){
      stop("Error on the basis generation")
    }

    D_test[,basis_subindex[[i]]] <- as.matrix(B_test_obj)

  }

  # Visualizing the basis

  # ==== COMMENTTED FUNCTIONS BELOW NOT RUN WITH IF NOT INSIDE THE FUNCTION
  # selected_var <- 1
  # D_subset <- D_train[,basis_subindex[[selected_var]]]
  # plot(NULL,ylim = range(D_subset), xlim = range(x_train_scale[,selected_var]), main = "Use BS: FALSE")
  # for(i in 1:ncol(D_subset)){
  #   points(x_train_scale[,selected_var], D_subset[,i], pch = 20, col = ggplot2::alpha(i,0.5))
  # }


  if(motrbart_bool){
    D_train <- x_train_scale
    D_test <- x_test_scale

    basis_size <- 1 # Change this value to the desired size of each sublist
    D_seq <- 1:ncol(D_train)  # Replace this with the columns of D

    # Creating a vector
    basis_subindex <- split(D_seq, rep(1:(length(D_seq) %/% basis_size), each = basis_size, length.out = length(D_seq)))
  }

  # R-th difference order matrix
  # if(dif_order!=0){
  #   D <- D_gen(p = ncol(D_train),n_dif = dif_order)
  # } else {
  #   D <- diag(nrow = ncol(D_train))
  # }

  # Scaling the y
  min_y <- min(y_train)
  max_y <- max(y_train)

  # Getting the min and max for each column
  min_x <- apply(x_train_scale,2,min)
  max_x <- apply(x_train_scale, 2, max)

  # Scaling "y"
  if(scale_bool){
    y_scale <- normalize_bart(y = y_train,a = min_y,b = max_y)

    # New update
    m_tilda <- mean(diag(tcrossprod(D_train)))
    # Maybe need to change that in the future

    # tau_mu <- 4*n_tree*(kappa^2)
    tau_mu <- 4*n_tree*(kappa^2)*(m_tilda)



  } else {
    y_scale <- y_train

    # New parameter update
    m_tilda <- mean(diag(tcrossprod(D_train)))
    tau_mu <- (4*n_tree*(kappa^2))/((max_y-min_y)^2)
    # tau_mu <- (4*n_tree*(kappa^2)*(m_tilda)*(nIknots-1))/((max_y-min_y)^2)
    tau_mu <- (4*n_tree*(kappa^2)*(m_tilda))/((max_y-min_y)^2)


  }


  # Getting the naive sigma value
  nsigma <- naive_sigma(x = x_train_scale,y = y_scale)

  # Calculating tau hyperparam
  a_tau <- df/2

  # Calculating lambda
  qchi <- stats::qchisq(p = 1-sigquant,df = df,lower.tail = 1,ncp = 0)
  lambda <- (nsigma*nsigma*qchi)/df
  d_tau <- (lambda*df)/2


  # Getting hyperparameters for \tau_beta_j
  # df <- 10
  a_tau_beta_j <- df/2
  sigquant_beta <- 0.99
  nsigma_beta <- tau_mu^(-1/2)

  # Calculating lambda
  qchi_beta <- stats::qchisq(p = 1-sigquant_beta,df = df,lower.tail = 1,ncp = 0)
  lambda_beta <- (nsigma_beta*nsigma_beta*qchi_beta)/df
  d_tau_beta_j <- (lambda_beta*df)/2



  # Visualising the prior
  rgamma(n = 1000,shape = a_tau_beta_j,rate = d_tau_beta_j) |> density() |> plot(main = "density prior \tau_beta_j")
  mean(rgamma(n = 1000,shape = a_tau_beta_j,rate = d_tau_beta_j))

  # Call the bart function
  tau_init <- nsigma^(-2)

  mu_init <- mean(y_scale)

  # Creating the vector that stores all trees
  all_tree_post <- vector("list",length = round(n_mcmc-n_burn))


  # =====================================================================
  # ========= From here I gonna initialise the BART function itself =====
  # =====================================================================
  n_post <- (n_mcmc-n_burn)
  all_trees <- vector("list", n_mcmc)
  all_betas <- vector("list",n_mcmc)
  tau_beta <- rep(tau_mu,NCOL(x_train_scale)) # In this first scenario we are going to work with a single value of \tau
  all_tau_beta <- matrix(NA, nrow = (n_mcmc), ncol = NCOL(x_train_scale))
  # all_delta <- numeric(n_mcmc)
  all_tau <- numeric(n_mcmc)

  all_y_hat <- matrix(NA,nrow = n_mcmc,ncol = nrow(x_train_scale))
  all_y_hat_test <- matrix(NA, nrow = n_mcmc, ncol = nrow(x_test_scale))
  all_trees_fit <- vector("list",n_mcmc)
  all_trees <- vector("list",n_mcmc)
  forest <- vector("list",n_tree)

  # Partial component pieces
  partial_train_fits <- vector("list", n_tree)

  proposal_outcomes <- setNames(data.frame(matrix(nrow = 0, ncol = 6)),
                                c("tree_number" , "proposal", "status","mcmc_iter", "new_tree_loglike", "old_tree_loglike"))
  all_train_indexes <- data.frame(matrix(data = NA,nrow = nrow(xcut_m),ncol = ncol(xcut_m)))

  # Gonna create a list of lists to store all the indexes for all split rules and cutpoints
  all_var_splits <- vector("list",ncol(x_train_scale))
  names(all_var_splits) <- colnames(x_train_scale)

  # Iterating over all possible x.columns
  for(i in 1:length(all_var_splits)){

    # Creating the dummy for a list of index to store all numeric split values
    all_cut_points <- vector("list", nrow(xcut_m))


    for(j in 1:length(all_cut_points)){

      # Getting the node indexes object
      left_train_list <- vector("list",length = 1L)
      names(left_train_list) <- "left_train"
      right_train_list <- vector("list",length = 1L)
      names(right_train_list) <- "right_train"
      left_test_list <- vector("list",length = 1L)
      names(left_test_list) <- "left_test"
      right_test_list <- vector("list",length = 1L)
      names(right_test_list) <- "right_test"

      node_index <- append(left_train_list, right_train_list) |>
        append(left_test_list) |> append(right_test_list)

      all_cut_points[[j]]$left_train <-  which(x_train_scale[,i] < xcut_m[j,i])
      all_cut_points[[j]]$right_train <-  which(x_train_scale[,i] >= xcut_m[j,i])
      all_cut_points[[j]]$left_test <-  which(x_test_scale[,i] < xcut_m[j,i])
      all_cut_points[[j]]$right_test <-  which(x_test_scale[,i] >= xcut_m[j,i])

    }

    all_var_splits[[i]] <- all_cut_points

  }

  # Creating the penalty matrix

  all_P <- replicate(NCOL(x_train_scale),
                     P_gen(D_train_ = B_train_obj,dif_order_ = dif_order,tau_mu_ = 1),
                     simplify = FALSE)

  P_train <- as.matrix(Matrix::bdiag(all_P))

  # Adjusting D_train
  # D_train <- do.call(cbind,replicate(NCOL(x_train_scale),D_train,simplify = FALSE))
  # D_test <- do.call(cbind,replicate(NCOL(x_train_scale),D_test,simplify = FALSE))


  #most of the functions
  data <- list(x_train = x_train_scale,
               x_test = x_test_scale,
               y_train = y_scale,
               xcut_m = xcut_m,
               D_train = D_train,
               D_test = D_test,
               dif_order = dif_order,
               alpha = alpha,
               beta = beta,
               basis_subindex = basis_subindex,
               all_var_splits = all_var_splits,
               n_tree = n_tree,
               tau_mu = tau_mu,
               tau = tau,
               a_tau = a_tau,
               d_tau = d_tau,
               tau_beta = tau_beta,
               # delta = delta,
               P = P_train,
               node_min_size = node_min_size,
               all_var = all_var,
               stump = stump,
               a_tau_beta_j = a_tau_beta_j,
               d_tau_beta_j = d_tau_beta_j)

  #   So to simply interepret the element all_var_splits each element correspond
  #to each variable. Afterwards each element corresponds to a cutpoint; Finally,
  #inside that level we would have the index for the the left and right nodes;

  # Initialing for storing post samples
  post <- 0

  # Initialising all the stumps
  for(k in 1:data$n_tree){
    forest[[k]] <- stump(data = data,init_var = k)
  }

  # and tree predictions
  trees_fit <- matrix(0,nrow = n_tree,ncol = nrow(x_train_scale))
  trees_fit_test <- matrix(0,nrow = n_tree, ncol  = nrow(x_test_scale))

  # For cases where the tree is greater than one;
  if(scale_init){
      if(n_tree>1){
        # Initial prediction
        for(i in 1:n_tree){
          trees_fit[i,] <- y_scale/n_tree
        }
      }
  }

  # Initialsing the loop
  for(i in 1:n_mcmc){

    # Initialising the partial train tree fits
    partial_train_fits[[i]] <- vector("list",data$n_tree)
    names(partial_train_fits[[i]]) <- paste0("tree",1:data$n_tree)

    # Initialising orogress bar
    progress <- i / n_mcmc * 100

    x1_pred <- numeric(nrow(x_train))

    for(t in 1:data$n_tree){

      # Calculating the partial residuals
      if(n_tree>1){
        partial_residuals <- y_scale-colSums(trees_fit[-t,,drop = FALSE])
      } else {
        partial_residuals <- y_scale
      }


      # Sample a verb
      verb <- sample(c("grow","prune", "change"), prob = c(0.3,0.3,0.4),size = 1)

      # Forcing to grow when only have a stump
      if(length(forest[[t]])==1){
        if(!data$all_var){
          verb <- sample(c("grow","change"),size = 1)
        } else {
          verb <- "grow"
        }
      }


      # Sampling a verb
      if(verb == "grow"){
        if(runif(n = 1)<0.5){
        forest[[t]] <- grow_predictors(tree = forest[[t]],
                            curr_part_res = partial_residuals,
                            data = data)
        } else {
          forest[[t]] <- grow(tree = forest[[t]],
                          curr_part_res = partial_residuals,
                          data = data)
        }
      } else if (verb == "prune"){
        if(runif(n = 1)<0.5){
          forest[[t]] <- prune_predictors(tree = forest[[t]],
                               curr_part_res = partial_residuals,
                               data = data)
        } else {
          forest[[t]] <- prune(tree = forest[[t]],
                                          curr_part_res = partial_residuals,
                                          data = data)
        }
      } else if (verb == "change"){

        if(runif(n = 1)<0.5){
          forest[[t]] <- change_predictors(tree = forest[[t]],
                                curr_part_res = partial_residuals,
                                data = data)
        } else {
          forest[[t]] <- change(tree = forest[[t]],
                                           curr_part_res = partial_residuals,
                                           data = data)
        }

      }

      # cat("Forest size:", (length(forest)),"\n")
      # cat("Forest verb:", verb,"\n")

      # Updating the betas
      forest[[t]] <- updateBetas(tree = forest[[t]],
                                 curr_part_res = partial_residuals,
                                 data = data)

      # Getting the predictions
      tree_predictions <- getPredictions(tree = forest[[t]],
                                         data = data)

      trees_fit[t,] <- rowSums(tree_predictions$y_train_hat)
      trees_fit_test[t,] <- rowSums(tree_predictions$y_hat_test)
      partial_train_fits[[t]] <- tree_predictions$y_train_hat

      # selected_var_ <- 9
      # plot(data$x_train[,selected_var_],tree_predictions$y_train_hat[,selected_var_])

      if(plot_preview){
          choose_dimension <- 1
          if(t==1){
            plot(x_train_scale[,choose_dimension],tree_predictions$y_train_hat[,choose_dimension], pch = 20, main = paste0("X",choose_dimension," partial pred"),ylim = range(y_scale),
                 col = ggplot2::alpha("black",0.2))
          } else {
            points(x_train_scale[,choose_dimension],tree_predictions$y_train_hat[,choose_dimension], pch=20, col = ggplot2::alpha(t,0.2))
          }
          trees_fit[t,] <- rowSums(tree_predictions$y_train_hat)
          trees_fit_test[t,] <- rowSums(tree_predictions$y_hat_test)
          partial_train_fits[[t]] <- tree_predictions$y_train_hat

          x1_pred <- x1_pred + tree_predictions$y_train_hat[,choose_dimension]
      }
    }

    if(plot_preview){
      points(x_train_scale[,choose_dimension],x1_pred, pch=20, col = "blue")
      x1_pred <- numeric(nrow(x_train))
    }

    # Getting final predcition
    y_hat <- colSums(trees_fit)
    y_hat_test <- colSums(trees_fit_test)


    # Seeing the results for the unidimensional cases.
    if(plot_preview){
        plot(x_train_scale,y_scale)
        for(plot_i in 1:n_tree){
          points(x_train_scale,trees_fit[plot_i,],pch=20,col = ggplot2::alpha(plot_i,0.2))
        }
        points(x_train_scale,y_hat,col = "blue",pch=20)
    }


    # Updating all other parameters
    if(update_tau_beta & data$all_var){
      data$tau_beta <- update_tau_betas_j(forest = forest,data = data)
    } else if (update_tau_beta){
      data$tau_beta <- update_tau_betas(forest = forest,data = data)
    }


    # Updating delta
    # data$delta <- update_delta(data = data)


    # Getting tau
    data$tau <- update_tau(y_train_hat = y_hat,
                           data = data)


    # Storing all predictions
    all_trees[[i]] <- forest
    all_tau[[i]] <- data$tau
    all_trees_fit[[i]] <- partial_train_fits
    all_y_hat[i,] <- y_hat
    all_y_hat_test[i,] <- y_hat_test
    all_tau_beta[i,] <- data$tau_beta
    # all_delta[i] <- data$delta


    # Print progress bar
    cat("\rProgress: [", paste(rep("=", floor(progress / 5)), collapse = ""),
        paste(rep(" ", floor((100 - progress) / 5)), collapse = ""),
        "] ", sprintf("%.1f%%", progress))

    # Flush the output
    flush.console()

    # Simulate some work
    Sys.sleep(0.1)

  }


  # Normalising elements
  all_tau_norm <- numeric(n_mcmc)
  all_tau_beta_norm <- matrix(NA,nrow = n_mcmc, NCOL(x_train_scale))
  all_trees_fit_norm <- vector("list",n_mcmc)
  all_y_hat_norm <- matrix(NA,nrow = nrow(all_y_hat),ncol = ncol(all_y_hat))
  all_y_hat_test_norm <- matrix(NA,nrow = nrow(all_y_hat_test),ncol = ncol(all_y_hat_test))

  # Returning to the original scale
  if(scale_bool){

    all_tau_norm <- all_tau/((max_y-min_y)^2)
    all_tau_beta_norm <- all_tau_beta/((max_y-min_y)^2)

    for(post_iter in 1:n_mcmc){

      # ==============
      # FIX THIS LATER
      # ==============

      for(tree_number in 1:n_tree){
        all_trees_fit_norm[[post_iter]][[tree_number]] <- unnormalize_bart(z = all_trees_fit[[post_iter]][[tree_number]],a = min_y,b = max_y)
      }
      all_y_hat_norm[post_iter,] <- unnormalize_bart(z = all_y_hat[post_iter,],a = min_y,b = max_y)
      all_y_hat_test_norm[post_iter, ] <- unnormalize_bart(z = all_y_hat_test[post_iter,],a = min_y,b = max_y)
    }
  } else {
    all_tau_norm <- all_tau
    all_tau_beta_norm <- all_tau_beta
    all_trees_fit_norm <-all_trees_fit

    all_y_hat_norm <- all_y_hat
    all_y_hat_test_norm <- all_y_hat_test
  }
  # plot(colMeans(all_y_hat_norm),y_train)
  # sqrt(crossprod((colMeans(all_y_hat_norm)-y_train))/n_)

  # ====== Few analyses from the results ======
  #           (Uncomment to run those)
  # ===========================================
  #

  # plot(all_tau,type = "l")
  # plot(all_tau,type = "l")
  #
  # y1_hat <- matrix(0,nrow = n_post,ncol = nrow(data$x_train))
  #
  # for(i in 1:nrow(y1_hat)){
  #   for(t in 1:n_tree){
  #     y1_hat[i,] <- all_trees_fit[[i]][[t]][,1]
  #   }
  # }
  #
  # plot(x_train_scale[,1],colMeans(y1_hat[1801:3000,,drop = FALSE]))
  # plot(colMeans(y1_hat[1801:3000,,drop = FALSE]),y_train)
  # # plot(x_train[,1],10 * sin(pi * x_train[, 1])) # For x1
  # # plot(x_train[,2],20 * (x_train[, 2] - 0.5)^2) # For x1
  #
  # plot(all_tau_beta, type = "l")
  # plot(all_delta, type = "l")
  # plot(x_train_scale,y_scale)
  # points(x_train_scale,colMeans(all_y_hat[1501:2000,]),pch= 20, col = "blue")
  # # ============================================
  #
  # curr <- 0
  # tree_lengths <- numeric()
  #
  # for(i in 1:length(all_trees)){
  #
  #     curr <- curr + 1
  #
  #     for(j in 1:n_tree){
  #         tree_lengths[curr] <- length(all_trees[[i]][[j]])
  #     }
  #
  # }




  # Return the list with all objects and parameters
  return(list(y_train_hat = all_y_hat_norm,
              y_test_hat = all_y_hat_test_norm,
              all_tau = all_tau_norm,
              all_tau_beta = all_tau_beta_norm,
              prior = list(n_tree = n_tree,
                           alpha = alpha,
                           beta = beta,
                           tau_mu = tau_mu,
                           a_tau = a_tau,
                           d_tau = d_tau),
              mcmc = list(n_mcmc = n_mcmc,
                          n_burn = n_burn,
                          all_trees = all_trees),
              data = list(x_train = x_train,
                          y_train = y_train,
                          D_train = D_train,
                          x_test = x_test,
                          D_test = D_test,
                          basis_subindex = basis_subindex)))

}




