# Creating a stump for a tree
stump <- function(data, init_var){

  if(init_var > NCOL(data$x_train)){
    init_var <- sample(1:NCOL(data$x_train),size = 1)
  }

  # Creating the base node
  node <- list()

  if(!data$all_var){
      node[["node0"]] <- list(
        # Creating the node number
        node_number = 0,
        isRoot = TRUE,
        # Creating a vector with the tranining index
        train_index = 1:nrow(data$x_train),
        test_index = 1:nrow(data$x_test),
        depth_node = 0,
        node_var = NA,
        node_cutpoint_index = NA,
        left = NA,
        right = NA,
        parent_node = NA,
        ancestors = init_var,
        terminal = TRUE,
        betas_vec = NULL
      )
  } else {
    node[["node0"]] <- list(
      # Creating the node number
      node_number = 0,
      isRoot = TRUE,
      # Creating a vector with the tranining index
      train_index = 1:nrow(data$x_train),
      test_index = 1:nrow(data$x_test),
      depth_node = 0,
      node_var = NA,
      node_cutpoint_index = NA,
      left = NA,
      right = NA,
      parent_node = NA,
      ancestors = 1:ncol(data$x_train),
      terminal = TRUE,
      betas_vec = NULL
    )
  }
  # Returning the node
  return(node)

}

# Get all the terminal nodes
get_terminals <- function(tree){

  # Return the name of the termianl nodes
  return(names(tree)[unlist(lapply(tree, function(x){x$terminal}),use.names =  TRUE)])
}

# Get nog terminal nodes
get_nogs <- function(tree){

  # Return the name of the termianl nodes
  non_terminal <- names(tree)[!unlist(lapply(tree, function(x){x$terminal}),use.names =  TRUE)]

  # In case there are non nonterminal nondes
  if(length(non_terminal)==0){
    return(non_terminal)
  }

  bool_nog <- vector("logical",length = length(non_terminal))
  for(i in 1:length(bool_nog)){
    # Checking if both children are terminal
    if( tree[[tree[[non_terminal[i]]]$left]]$terminal & tree[[tree[[non_terminal[i]]]$right]]$terminal) {
      bool_nog[i] <- TRUE
    }
  }

  return(  non_terminal[bool_nog])
}

# Getting the maximum node index number
get_max_node <- function(tree){

  # Return the name of the termianl nodes
  return(max(unlist(lapply(tree, function(x){x$node_number}),use.names =  TRUE)))
}




# A function to calculate the loglikelihood
nodeLogLike <- function(curr_part_res,
                        ancestors,
                        index_node,
                        data){

  # Subsetting the residuals
  curr_part_res_leaf <- curr_part_res[index_node]

  # Getting the number of observationsin the terminal node
  n_leaf <- length(index_node)
  d_basis <- length(ancestors)
  ones <- matrix(1,nrow = n_leaf)
  # Getting the index of all covariates for each basis from each ancestors
  D_subset_index <- unlist(data$basis_subindex[ancestors])


  # Getting the p_{tell} i.e: number of betas of the current terminal node
  d_basis <- length(D_subset_index)
  D_leaf <- data$D_train[index_node,D_subset_index, drop = FALSE]

  if(NCOL(D_leaf)==0){
    stop(" Node Log-likelihood: No ancestors")
  }

  # Using the Andrew's approach I would have
  mean_aux <- rep(0,length(curr_part_res_leaf))

  if(data$all_var){
    diag_tau_beta_inv <- diag(x = rep(data$tau_beta, each = NCOL(D_leaf)/NCOL(data$x_train)))
  } else {
    diag_tau_beta_inv <- diag(x = unique(data$tau_beta), nrow = NCOL(D_leaf))
  }

  if(data$dif_order==0){
      cov_aux <- diag(x = (data$tau^(-1)),nrow = n_leaf) + D_leaf%*%tcrossprod(diag_tau_beta_inv,D_leaf)
  } else {
      cov_aux <- diag(x = (data$tau^(-1)),nrow = n_leaf) + D_leaf%*%tcrossprod(data$P[D_subset_index,D_subset_index],D_leaf)

  }

  result <- mvnfast::dmvn(X = curr_part_res_leaf,mu = mean_aux,sigma = cov_aux,log = TRUE)


  return(c(result))

}



# Grow a tree
grow <- function(tree,
                 curr_part_res,
                 data){

  # Getting the maximum index number
  max_index <- get_max_node(tree)

  # Sampling a terminal node
  terminal_nodes <- get_terminals(tree)
  n_t_nodes <- length(terminal_nodes)
  nog_nodes <- get_nogs(tree)
  n_nog_nodes <- length(nog_nodes)
  g_node_name <- sample(terminal_nodes,size = 1)
  g_node <- tree[[g_node_name]]


  valid_terminal_node <- TRUE
  valid_count <- 0
  while(valid_terminal_node){
    # Convinience while to avoid terminal nodes of 2
    # Sample a split var
    p_var <- sample(1:ncol(data$x_train),size = 1)

    # Selecting an available cutpoint from this terminal node
    valid_range_grow <- range(data$x_train[g_node$train_index,p_var])

    # Case of invalid range
    if(length(valid_range_grow)==0){
      return(tree)
    }

    # Subsetting the indexes of
    valid_cutpoint <- which(data$xcut_m[,p_var]>valid_range_grow[1] & data$xcut_m[,p_var]<valid_range_grow[2])

    # When there's no valid cutpoint on the sampled terminal node
    if(length(valid_cutpoint)==0){
      return(tree)
    }

    # Getting which cutpoints are valid and sample onde index
    sample_cutpoint <- sample(valid_cutpoint,
                              size = 1)

    # Getting the left & right index
    left_index  <- data$all_var_splits[[p_var]][[sample_cutpoint]]$left_train[data$all_var_splits[[p_var]][[sample_cutpoint]]$left_train %in% g_node$train_index]
    right_index  <- data$all_var_splits[[p_var]][[sample_cutpoint]]$right_train[data$all_var_splits[[p_var]][[sample_cutpoint]]$right_train %in% g_node$train_index]

    left_test_index  <- data$all_var_splits[[p_var]][[sample_cutpoint]]$left_test[data$all_var_splits[[p_var]][[sample_cutpoint]]$left_test %in% g_node$test_index]
    right_test_index  <- data$all_var_splits[[p_var]][[sample_cutpoint]]$right_test[data$all_var_splits[[p_var]][[sample_cutpoint]]$right_test %in% g_node$test_index]



    # Verifying that the correct number was used
    if((length(left_index)+length(right_index))!=length(g_node$train_index)){
      stop("Something went wrong here --- train grown index doest match")
    }

    if((length(left_test_index)+length(right_test_index))!=length(g_node$test_index)){
      stop("Something went wrong here --- test grown index doest match")
    }

    if( (length(left_index) > data$node_min_size) & (length(right_index)>data$node_min_size)){
      # Getting out of the while
      break
    } else {

      # Adding one to the counter
      valid_count = valid_count + 1

      # Stop trying to search for a valid cutpoint
      if(valid_count > 2) {
        valid_terminal_node = FALSE
        return(tree)
      }
    }
  }

  # For convinience we are going to avoid terminal nodes less than 2
  if( (length(left_index)<2) || (length(right_index) < 2)) {
    stop("Error of invalid terminal node")
  }

  # Calculating loglikelihood for the grown node, the left and the right node

  g_loglike <- nodeLogLike(curr_part_res = curr_part_res,
                           ancestors = unique(g_node$ancestors),
                           index_node = g_node$train_index,
                           data = data)


  left_loglike <-  nodeLogLike(curr_part_res = curr_part_res,
                               index_node = left_index,
                               ancestors = unique(c(g_node$ancestors,p_var)),
                               data = data)

  right_loglike <-  nodeLogLike(curr_part_res = curr_part_res,
                                index_node = right_index,
                                ancestors = unique(c(g_node$ancestors,p_var)),
                                data = data)

  # Calculating the prior
  prior_loglike <- log(data$alpha*(1+g_node$depth_node)^(-data$beta)) + # Prior of the grown node becoming nonterminal
    2*log(1-data$alpha*(1+g_node$depth_node+1)^(-data$beta)) - # plus the prior of the two following nodes being terminal
    log(1-data$alpha*(1+g_node$depth_node)^(-data$beta)) # minus the probability of the grown node being terminal

  # Transition prob
  log_trasition_prob  = log(0.3/(n_nog_nodes+1))-log(0.3/n_t_nodes)

  # Calculating the acceptance probability
  acceptance <- exp(-g_loglike+left_loglike+right_loglike+prior_loglike+log_trasition_prob)

  if(data$stump) {
    acceptance <- acceptance*(-1)
  }

  # Getting the training the left and the right index for the the grown node
  if(stats::runif(n = 1)<acceptance){


    # Verifying if uses all variables or not
    if(!data$all_var){
        left_node <- list(node_number = max_index+1,
                          isRoot = FALSE,
                          train_index = left_index,
                          test_index = left_test_index,
                          depth_node = g_node$depth_node+1,
                          node_var = p_var,
                          node_cutpoint_index = sample_cutpoint,
                          left = NA,
                          right = NA,
                          parent_node = g_node_name,
                          ancestors = c(g_node$ancestors,p_var),
                          terminal = TRUE,
                          betas_vec = rep(0,ncol(data$D_train)))

        right_node <- list(node_number = max_index+2,
                           isRoot = FALSE,
                           train_index = right_index,
                           test_index = right_test_index,
                           depth_node = g_node$depth_node+1,
                           node_var = p_var,
                           node_cutpoint_index = sample_cutpoint,
                           left = NA,
                           right = NA,
                           parent_node = g_node_name,
                           ancestors = c(g_node$ancestors,p_var),
                           terminal = TRUE,
                           betas_vec = rep(0,ncol(data$D_train)))
    } else {

      if(!identical(g_node$ancestors,1:NCOL(data$x_train))){
        stop("No match in ancestors")
      }

      left_node <- list(node_number = max_index+1,
                        isRoot = FALSE,
                        train_index = left_index,
                        test_index = left_test_index,
                        depth_node = g_node$depth_node+1,
                        node_var = p_var,
                        node_cutpoint_index = sample_cutpoint,
                        left = NA,
                        right = NA,
                        parent_node = g_node_name,
                        ancestors = g_node$ancestors,
                        terminal = TRUE,
                        betas_vec = rep(0,ncol(data$D_train)))

      right_node <- list(node_number = max_index+2,
                         isRoot = FALSE,
                         train_index = right_index,
                         test_index = right_test_index,
                         depth_node = g_node$depth_node+1,
                         node_var = p_var,
                         node_cutpoint_index = sample_cutpoint,
                         left = NA,
                         right = NA,
                         parent_node = g_node_name,
                         ancestors = g_node$ancestors,
                         terminal = TRUE,
                         betas_vec = rep(0,ncol(data$D_train)))
    }

    # Modifying the current node
    tree[[g_node_name]]$left = paste0("node",max_index+1)
    tree[[g_node_name]]$right = paste0("node",max_index+2)
    tree[[g_node_name]]$terminal = FALSE

    tree[[paste0("node",max_index+1)]] <- left_node
    tree[[paste0("node",max_index+2)]] <- right_node


  } else {

    # Do nothing

  }

  # Return the new tree
  return(tree)
}

# Grow a tree
grow_predictors <- function(tree,
                 curr_part_res,
                 data){

  if(data$all_var){
    stop("This VERB doesn't allow use all vars")
  }

  # Getting the maximum index number
  max_index <- get_max_node(tree)

  # Sampling a terminal node
  terminal_nodes <- get_terminals(tree)
  n_t_nodes <- length(terminal_nodes)
  nog_nodes <- get_nogs(tree)
  n_nog_nodes <- length(nog_nodes)

  if(n_nog_nodes==0){
    g_node_name <- "node0"
  } else {
    g_node_name <- sample(nog_nodes,size = 1) # VERIFY THE CASE WHEN NOG == 0
  }

  g_node <- tree[[g_node_name]]


  valid_terminal_node <- TRUE
  valid_count <- 0

  candidates <- (1:NCOL(data$x_train))[!(1:NCOL(data$x_train) %in% g_node$ancestors)]

  if(length(candidates)==0){
    return(tree)
  } else {
    p_var <- sample(x = candidates,size = 1)
  }


  # ==========================
  # For the case for the STUMP
  # ==========================

  if(length(tree)==1){

          # Calculating loglikelihood for the grown node, the left and the right node
          g_loglike <- nodeLogLike(curr_part_res = curr_part_res,
                                   ancestors = unique(g_node$ancestors),
                                   index_node = g_node$train_index,
                                   data = data)


          new_g_loglike <- nodeLogLike(curr_part_res = curr_part_res,
                                       ancestors = unique(c(g_node$ancestors,p_var)),
                                       index_node = g_node$train_index,
                                       data = data)


          # Transition prob (MAY NEED TO CHANGE IN THE FUTURE)
          log_trasition_prob  = log(0.3/(n_nog_nodes+1))-log(0.3/n_t_nodes)

          # Calculating the acceptance probability
          acceptance <- exp(-g_loglike+new_g_loglike)

          if(data$stump) {
            acceptance <- acceptance*(-1)
          }

          # Getting the training the left and the right index for the the grown node
          if(stats::runif(n = 1)<acceptance){
            tree[[g_node_name]]$ancestors <- unique(c(g_node$ancestors,p_var))
          }

      } else { # FOR THE CASE THAT ISN'T STUMP

        # Getting the new ancestors
        new_ancestor <- c(unique(tree[[g_node$left]]$ancestors),p_var)

        # Calculating loglikelihood for the grown node, the left and the right node
        g_loglike_left <- nodeLogLike(curr_part_res = curr_part_res,
                                 ancestors = unique(tree[[g_node$left]]$ancestors),
                                 index_node = tree[[g_node$left]]$train_index,
                                 data = data)

        g_loglike_right <- nodeLogLike(curr_part_res = curr_part_res,
                                      ancestors = unique(tree[[g_node$right]]$ancestors),
                                      index_node = tree[[g_node$right]]$train_index,
                                      data = data)

        # Calculating loglikelihood for the grown node, the left and the right node
        new_g_loglike_left <- nodeLogLike(curr_part_res = curr_part_res,
                                      ancestors = new_ancestor,
                                      index_node = tree[[g_node$left]]$train_index,
                                      data = data)

        new_g_loglike_right <- nodeLogLike(curr_part_res = curr_part_res,
                                       ancestors = new_ancestor,
                                       index_node = tree[[g_node$right]]$train_index,
                                       data = data)



        # Calculating the acceptance probability
        acceptance <- exp(-g_loglike_left-g_loglike_right+new_g_loglike_left+new_g_loglike_right)

        if(data$stump) {
          acceptance <- acceptance*(-1)
        }

        # Getting the training the left and the right index for the the grown node
        if(stats::runif(n = 1)<acceptance){
          tree[[g_node$left]]$ancestors <- new_ancestor
          tree[[g_node$right]]$ancestors <- new_ancestor
        }

      }



  # Return the new tree
  return(tree)
}


# Pruning a tree
prune <- function(tree,
                  curr_part_res,
                  data){


  # Getting the maximum index number
  max_index <- get_max_node(tree)

  # Sampling a terminal node
  terminal_nodes <- get_terminals(tree)
  n_t_nodes <- length(terminal_nodes)
  nog_nodes <- get_nogs(tree)
  n_nog_nodes <- length(nog_nodes)

  # Just in case to avoid errors
  if(n_nog_nodes==0){
    return(tree)
  }

  # Selecting a node to be pruned
  p_node_name <- sample(nog_nodes,size = 1)
  p_node <- tree[[p_node_name]]

  # Getting the indexes from the left and right children from the pruned node
  children_left_index <- tree[[p_node$left]]$train_index
  children_right_index <- tree[[p_node$right]]$train_index
  children_left_ancestors <- tree[[p_node$left]]$ancestors
  children_right_ancestors <- tree[[p_node$right]]$ancestors

  # Calculating loglikelihood for the grown node, the left and the right node

  p_loglike <- nodeLogLike(curr_part_res = curr_part_res,
                           index_node = p_node$train_index,
                           data = data,
                           ancestors = unique(p_node$ancestors))


  p_left_loglike <-  nodeLogLike(curr_part_res = curr_part_res,
                                 index_node = children_left_index,
                                 ancestors = unique(children_left_ancestors),
                                 data = data)

  p_right_loglike <-  nodeLogLike(curr_part_res = curr_part_res,
                                  index_node = children_right_index,
                                  ancestors = unique(children_right_ancestors),
                                  data = data)

  # Calculating the prior
  prior_loglike <- log(1-data$alpha*(1+p_node$depth_node)^(-data$beta)) - # Prior of the new terminal node
    log(data$alpha*(1+p_node$depth_node)^(-data$beta)) - # Prior of the grown node becoming nonterminal
    2*log(1-data$alpha*(1+p_node$depth_node+1)^(-data$beta))  # plus the prior of the two following nodes being terminal
  # minus the probability of the grown node being terminal

  # Transition prob
  log_trasition_prob  = log(0.3/(n_t_nodes))-log(0.3/n_nog_nodes)

  # Calculating the acceptance probability
  acceptance <- exp(p_loglike-p_left_loglike-p_right_loglike+prior_loglike+log_trasition_prob)

  # Getting the training the left and the right index for the the grown node
  if(stats::runif(n = 1)<acceptance){

    # Erasing the terminal nodes
    tree[[p_node$left]] <- NULL
    tree[[p_node$right]] <- NULL

    # Modifying back the pruned node
    tree[[p_node_name]]$left <- NA
    tree[[p_node_name]]$right <- NA
    tree[[p_node_name]]$terminal <- TRUE

  } else {
    # Do nothing
  }

  # Return the new tree
  return(tree)

}

# Pruning a tree
prune_predictors <- function(tree,
                  curr_part_res,
                  data){


  # Getting the maximum index number
  max_index <- get_max_node(tree)

  # Sampling a terminal node
  terminal_nodes <- get_terminals(tree)
  n_t_nodes <- length(terminal_nodes)
  nog_nodes <- get_nogs(tree)
  n_nog_nodes <- length(nog_nodes)

  # Just in case to avoid errors
  if(n_nog_nodes==0){
    return(tree)
  }

  # Selecting a node to be pruned
  p_node_name <- sample(nog_nodes,size = 1)
  p_node <- tree[[p_node_name]]

  # Getting the indexes from the left and right children from the pruned node
  if(length(candidates)==1){
    return(tree)
  } else {
    p_var_prune <- sample(tree[[p_node$left]]$ancestors,size = 1)
  }

  new_ancestors <- tree[[p_node$left]]$ancestors[!(tree[[p_node$left]]$ancestors %in% p_var_prune)]

  # Calculating loglikelihood for the grown node, the left and the right node
  p_left_loglike <- nodeLogLike(curr_part_res = curr_part_res,
                           index_node = tree[[p_node$left]]$train_index,
                           data = data,
                           ancestors = tree[[p_node$left]]$ancestors)

  p_right_loglike <- nodeLogLike(curr_part_res = curr_part_res,
                                 index_node = tree[[p_node$right]]$train_index,
                                 data = data,
                                 ancestors = tree[[p_node$right]]$ancestors)

  new_p_left_loglike <-  nodeLogLike(curr_part_res = curr_part_res,
                                 index_node =  tree[[p_node$left]]$train_index,
                                 data = data,
                                 ancestors = new_ancestors)

  new_p_right_loglike <- nodeLogLike(curr_part_res = curr_part_res,
                                 index_node = tree[[p_node$right]]$train_index,
                                 data = data,
                                 ancestors = new_ancestors)


  # (MAYBE NEED TO CHANGE THE TRANSITION) Transition prob
  log_trasition_prob  = log(0.3/(n_t_nodes))-log(0.3/n_nog_nodes)

  # Calculating the acceptance probability
  acceptance <- exp(-p_left_loglike-p_right_loglike+new_p_left_loglike+new_p_right_loglike)

  # Getting the training the left and the right index for the the grown node
  if(stats::runif(n = 1)<acceptance){

    # Erasing the terminal nodes
    tree[[p_node$left]]$ancestors <- new_ancestors
    tree[[p_node$right]]$ancestors <- new_ancestors

  } else {
    # Do nothing
  }

  # Return the new tree
  return(tree)

}

change_predictors <- function(tree = tree,
             curr_part_res = curr_part_res,
             data = data){

  # Getting the stump
  if(length(tree)==1){
    c_node <- tree$node0
    change_candidates <- which(!(1:NCOL(data$x_train) %in% c_node$ancestors))
  } else {
    # Getting the maximum index number
    max_index <- get_max_node(tree)

    # Sampling a terminal node
    terminal_nodes <- get_terminals(tree)
    n_t_nodes <- length(terminal_nodes)
    nog_nodes <- get_nogs(tree)
    n_nog_nodes <- length(nog_nodes)

    # Just in case to avoid errors
    if(n_nog_nodes==0){
      return(tree)
    }

    # Selecting a node to be pruned
    c_node_name <- sample(nog_nodes,size = 1)
    c_node <- tree[[c_node_name]]
    # Proposing a change to the stump
    change_candidates <- which(!(1:NCOL(data$x_train) %in% tree[[c_node$left]]$ancestors))
  }

  # In case there's other proposal trees (only for 1-d case)
  if(length(change_candidates)==0){
    # Gettina grown tree
    grown_tree <- grow(tree = tree,
         curr_part_res = curr_part_res,
         data = data)
    return(grown_tree)
  }

  # Doing for the stump case
  if(length(tree)==1){

    new_ancestor <- sample(change_candidates,size = 1)

    stump_loglikelihood <- nodeLogLike(curr_part_res = curr_part_res,
                                           ancestors = tree$node0$ancestors,
                                           index_node = c_node$train_index,
                                           data = data)

    new_stump_loglikelihood <- nodeLogLike(curr_part_res = curr_part_res,
                             ancestors = new_ancestor,
                             index_node = tree$node0$train_index,
                             data = data)

    # Calculating the acceptance probability
    acceptance <- exp(-stump_loglikelihood+new_stump_loglikelihood)

    # Getting the training the left and the right index for the the grown node
    if(stats::runif(n = 1)<acceptance){

      # Erasing the terminal nodes
      tree$node0$ancestors <- new_ancestor

    } else {
      # Do nothing
    }

  } else { # When we do not have an stump

    # Need to select one from the available
    single_new_ancestor <- sample(candidates,size = 1)
    replace_old_ancestor_index <- sample(1:length(tree[[c_node$left]]$ancestors),size = 1)
    new_ancestor <- tree[[c_node$left]]$ancestors
    new_ancestor[replace_old_ancestor_index] <- single_new_ancestor

    # Calculating loglikelihoods
    left_c_loglikelihood <- nodeLogLike(curr_part_res = curr_part_res,
                                        ancestors = tree[[c_node$left]]$ancestors,
                                        index_node = tree[[c_node$left]]$train_index,
                                        data = data)

    right_c_loglikelihood <- nodeLogLike(curr_part_res = curr_part_res,
                                        ancestors = tree[[c_node$right]]$ancestors,
                                        index_node = tree[[c_node$right]]$train_index,
                                        data = data)

    new_left_c_loglikelihood <- nodeLogLike(curr_part_res = curr_part_res,
                                            ancestors = new_ancestor,
                                            index_node = tree[[c_node$left]]$train_index,
                                            data = data)

    new_right_c_loglikelihood <- nodeLogLike(curr_part_res = curr_part_res,
                                            ancestors = new_ancestor,
                                            index_node = tree[[c_node$right]]$train_index,
                                            data = data)


    # Calculating the acceptance probability
    acceptance <- exp(-left_c_loglikelihood-right_c_loglikelihood + new_left_c_loglikelihood+new_right_c_loglikelihood)

    # Getting the training the left and the right index for the the grown node
    if(stats::runif(n = 1)<acceptance){

      # Erasing the terminal nodes
      tree[[c_node$left]]$ancestors <- new_ancestor
      tree[[c_node$right]]$ancestors <- new_ancestor

    } else {
      # Do nothing
    }


  }


  # Returning the new tree
  return(tree)


}


# Change a tree
change <- function(tree,
                   curr_part_res,
                   data){

  # Changing the stump
  if(length(tree)==1){
    change_stump_obj <- change_predictors(tree = tree,
                                 curr_part_res = curr_part_res,
                                 data = data)
    return(change_stump_obj)
  }

  # Sampling a terminal node
  nog_nodes <- get_nogs(tree)
  n_nog_nodes <- length(nog_nodes)
  c_node_name <- sample(nog_nodes,size = 1)
  c_node <- tree[[c_node_name]]


  valid_terminal_node <- TRUE
  valid_count <- 0


  while(valid_terminal_node){
    # Convinience while to avoid terminal nodes of 2
    # Sample a split var
    p_var <- sample(1:ncol(data$x_train),size = 1)

    # Selecting an available cutpoint from this terminal node
    valid_range_grow <- range(data$x_train[c_node$train_index,p_var])

    # Subsetting the indexes of
    valid_cutpoint <- which(data$xcut_m[,p_var]>valid_range_grow[1] & data$xcut_m[,p_var]<valid_range_grow[2])

    # When there's no valid cutpoint on the sampled terminal node
    if(length(valid_cutpoint)==0){
      return(tree)
    }

    # Getting which cutpoints are valid and sample onde index
    sample_cutpoint <- sample(valid_cutpoint,
                              size = 1)

    # Getting the left & right index
    left_index  <- data$all_var_splits[[p_var]][[sample_cutpoint]]$left_train[data$all_var_splits[[p_var]][[sample_cutpoint]]$left_train %in% c_node$train_index]
    right_index  <- data$all_var_splits[[p_var]][[sample_cutpoint]]$right_train[data$all_var_splits[[p_var]][[sample_cutpoint]]$right_train %in% c_node$train_index]

    left_test_index  <- data$all_var_splits[[p_var]][[sample_cutpoint]]$left_test[data$all_var_splits[[p_var]][[sample_cutpoint]]$left_test %in% c_node$test_index]
    right_test_index  <- data$all_var_splits[[p_var]][[sample_cutpoint]]$right_test[data$all_var_splits[[p_var]][[sample_cutpoint]]$right_test %in% c_node$test_index]



    # Verifying that the correct number was used
    if((length(left_index)+length(right_index))!=length(c_node$train_index)){
      stop("Something went wrong here --- train grown index doest match")
    }

    if((length(left_test_index)+length(right_test_index))!=length(c_node$test_index)){
      stop("Something went wrong here --- test grown index doest match")
    }

    # Avoiding having terminal nodes with just one observation
    if( (length(left_index) > data$node_min_size) & (length(right_index)>data$node_min_size)){
      # Getting out of the while
      break
    } else {

      # Adding one to the counter
      valid_count = valid_count + 1

      # Stop trying to search for a valid cutpoint
      if(valid_count > 2) {
        valid_terminal_node = FALSE
        return(tree)
      }
    }
  }

  # For convinience we are going to avoid terminal nodes less than 2
  if( (length(left_index)<2) || (length(right_index) < 2)) {
    stop("Error of invalid terminal node")
  }
  # Calculating loglikelihood for the new changed nodes and the old ones
  c_loglike_left <- nodeLogLike(curr_part_res = curr_part_res,
                                index_node = tree[[c_node$left]]$train_index,
                                ancestors = unique(tree[[c_node$left]]$ancestors),
                                data = data)


  c_loglike_right <-  nodeLogLike(curr_part_res = curr_part_res,
                                  index_node = tree[[c_node$right]]$train_index,
                                  ancestors = unique(tree[[c_node$right]]$ancestors),
                                  data = data)

  # Calculating a new ancestors left and right
  old_p_var <- tree[[c_node$left]]$node_var

  # Storing new left and right ancestors
  if(!data$all_var){
    new_left_ancestors <- tree[[c_node$left]]$ancestors
    new_left_ancestors[length(new_left_ancestors)] <- p_var

    new_right_ancestors <- tree[[c_node$right]]$ancestors
    new_right_ancestors[length(new_right_ancestors)] <- p_var
  } else {
    new_left_ancestors <- new_right_ancestors <- 1:NCOL(data$x_train)
  }

  new_c_loglike_left <-  nodeLogLike(curr_part_res = curr_part_res,
                                     index_node = left_index,
                                     ancestors = unique(new_left_ancestors),
                                     data = data)

  new_c_loglike_right <-  nodeLogLike(curr_part_res = curr_part_res,
                                      index_node = right_index,
                                      ancestors = unique(new_right_ancestors),
                                      data = data)


  # Calculating the acceptance probability
  acceptance <- exp(new_c_loglike_left+new_c_loglike_right-c_loglike_left-c_loglike_right)

  # Getting the training the left and the right index for the the grown node
  if(stats::runif(n = 1,min = 0,max = 1)<acceptance){

    # Updating the left and the right node
    # === Left =====
    tree[[c_node$left]]$node_var <- p_var
    tree[[c_node$left]]$node_cutpoint_index <- sample_cutpoint
    tree[[c_node$left]]$train_index <- left_index
    tree[[c_node$left]]$test_index <- left_test_index
    tree[[c_node$left]]$ancestors <- new_left_ancestors

    #==== Right ====
    tree[[c_node$right]]$node_var <- p_var
    tree[[c_node$right]]$node_cutpoint_index <- sample_cutpoint
    tree[[c_node$right]]$train_index <- right_index
    tree[[c_node$right]]$test_index <- right_test_index
    tree[[c_node$right]]$ancestors <- new_right_ancestors

  } else {

    # Do nothing

  }

  # Return the new tree
  return(tree)

}




# ============
# Update Betas
# ============
updateBetas <- function(tree,
                        curr_part_res,
                        data){


  # Getting the terminals
  t_nodes_names <- get_terminals(tree)


  for(i in 1:length(t_nodes_names)){


    # Select the current terminal node
    cu_t <- tree[[t_nodes_names[i]]]

    res_leaf <- matrix(curr_part_res[cu_t$train_index], ncol=1)

    # Creatinga  vector of zeros for betas_vec
    tree[[t_nodes_names[[i]]]]$betas_vec <- rep(0,ncol(data$D_train))

    # Selecting the actually parameters subsetting
    leaf_basis_subindex <- unlist(data$basis_subindex[unique(cu_t$ancestors)]) # Recall to the unique() here too
    basis_dim <- length(leaf_basis_subindex)
    D_leaf <- data$D_train[cu_t$train_index,leaf_basis_subindex, drop = FALSE]
    n_leaf <- length(cu_t$train_index)
    diag_leaf <- diag(nrow = n_leaf)
    diag_basis <- diag(nrow = basis_dim)


    #  Calculating the quantities need to the posterior of \beta
    b_ <- crossprod(D_leaf,res_leaf)
    data_tau_beta_diag <- rep(data$tau_beta, NCOL(D_leaf))
    Q_ <- (crossprod(D_leaf) + diag(data_tau_beta_diag/data$tau, nrow = NCOL(D_leaf)))
    Q_inv_ <- chol2inv(chol(Q_))
    # Q_inv_ <- solve(Q_)

    # tree[[t_nodes_names[i]]]$betas_vec[leaf_basis_subindex] <- c(keefe_mvn_sampler(b = b_,Q = Q_))
    tree[[t_nodes_names[i]]]$betas_vec[leaf_basis_subindex] <- mvnfast::rmvn(n = 1,mu = Q_inv_%*%b_,sigma = (data$tau^(-1))*Q_inv_)
  }

  # Returning the tree
  return(tree)

}

# =================
# Update \tau_betas
# =================
update_tau_betas_j <- function(forest,
                             data){


  # if(data$dif_order!=0){
  #   stop("Do not update tau_beta for peanalised version yet")
  # }

  # Setting some default hyperparameters
  # a_tau_beta <- d_tau_beta <- 0.1
  # Setting some default hyperparameters
  a_tau_beta <- data$a_tau_beta_j
  d_tau_beta <- data$d_tau_beta_j

  tau_b_shape <- 0.0
  tau_b_rate <- 0.0


  tau_b_shape <- numeric(NCOL(data$x_train))
  tau_b_rate <- numeric(NCOL(data$x_train))
  tau_beta_vec_aux <- numeric(NCOL(data$x_train))

  # Iterating over all trees
  for(i in 1:length(forest)){

    # Getting terminal nodes
    t_nodes_names <- get_terminals(forest[[i]])
    n_t_nodes <- length(t_nodes_names)

    # Iterating over the terminal nodes
    for(j in 1:length(t_nodes_names)){

      cu_t <- forest[[i]][[t_nodes_names[j]]]

      for(var_  in 1:NCOL(data$x_train)){

            # Getting ht leaf basis
            leaf_basis_subindex <- unlist(data$basis_subindex[var_]) # Recall to the unique() function here
            p_ <- length(leaf_basis_subindex)
            betas_mat_ <- matrix(cu_t$betas_vec[leaf_basis_subindex],nrow = p_)
            if(!is.null(cu_t$betas_vec)){
              tau_b_shape[var_] <- tau_b_shape[var_] + p_
              tau_b_rate[var_] <- tau_b_rate[var_] + c(crossprod(betas_mat_,crossprod(data$P[leaf_basis_subindex,leaf_basis_subindex, drop = FALSE],betas_mat_)))
            }

      }

    }


  }

  for(j in 1:NCOL(data$x_train)){
    tau_beta_vec_aux[j] <- rgamma(n = 1,
                               shape = 0.5*tau_b_shape[j] + a_tau_beta,
                               rate = 0.5*tau_b_rate[j] + d_tau_beta)

  }

  return(tau_beta_vec_aux)

}


update_tau_betas <- function(forest,
                             data){

  if(data$dif_order!=0){
    stop("Do not update tau_beta for peanalised version yet")
  }

  # Setting some default hyperparameters
  a_tau_beta <- d_tau_beta <- 0.1
  tau_b_shape <- 0.0
  tau_b_rate <- 0.0


  # Iterating over all trees
  for(i in 1:length(forest)){

    # Getting terminal nodes
    t_nodes_names <- get_terminals(forest[[i]])
    n_t_nodes <- length(t_nodes_names)

    # Iterating over the terminal nodes
    for(j in 1:length(t_nodes_names)){

      cu_t <- forest[[i]][[t_nodes_names[j]]]
      leaf_basis_subindex <- unlist(data$basis_subindex[unique(cu_t$ancestors)]) # Recall to the unique() function here

      if(!is.null(cu_t$betas_vec)){
        tau_b_shape <- tau_b_shape + length(leaf_basis_subindex)
        tau_b_rate <- tau_b_rate + c(crossprod(cu_t$betas_vec[leaf_basis_subindex]))
      }

    }


    tau_beta_vec_aux <- rgamma(n = 1,
                               shape = 0.5*tau_b_shape + a_tau_beta,
                               rate = 0.5*tau_b_rate + d_tau_beta)
  }

  return(tau_beta_vec_aux)

}


# ===================
# Updating the \delta
# ===================

# A function to get predictions
getPredictions <- function(tree,
                           data){

  # Creating the vector to hold the values of the prediction
  y_hat <- matrix(0, nrow = nrow(data$x_train), ncol = ncol(data$x_train))
  y_hat_test <- matrix(0,nrow(data$x_test), ncol = ncol(data$x_test))

  # Getting terminal nodes
  t_nodes <- get_terminals(tree = tree)
  n_t_nodes <- length(t_nodes)

  for(i in 1:n_t_nodes){

    leaf_train_index <- tree[[t_nodes[i]]]$train_index
    leaf_test_index <- tree[[t_nodes[i]]]$test_index
    leaf_ancestors <- unique(tree[[t_nodes[[i]]]]$ancestors) # recall the unique() argument here
    leaf_basis_subindex <- data$basis_subindex[leaf_ancestors]

    # Test unit
    if(length(leaf_ancestors)!=length(leaf_basis_subindex)){
      stop("Error on the getPredictions function")
    }
    # Only add the marginal effects if the variables are within that terminal node
    if(length(leaf_basis_subindex)!=0){
      for(k in 1:length(leaf_basis_subindex)){

        y_hat[leaf_train_index,leaf_ancestors[k]] <- y_hat[leaf_train_index,leaf_ancestors[k]] + data$D_train[leaf_train_index,leaf_basis_subindex[[k]], drop = FALSE]%*%tree[[t_nodes[i]]]$betas_vec[leaf_basis_subindex[[k]]]
        y_hat_test[leaf_test_index,leaf_ancestors[k]] <- y_hat_test[leaf_test_index,leaf_ancestors[k]] + data$D_test[leaf_test_index,leaf_basis_subindex[[k]], drop = FALSE]%*%tree[[t_nodes[i]]]$betas_vec[leaf_basis_subindex[[k]]]

      }
    }

  }

  # Returning both training and test set predictions
  return(list(y_train_hat = y_hat,
              y_hat_test = y_hat_test))

}

# Updating tau
update_tau <- function(y_train_hat,
                       data){

  # Sampling a tau value
  n_ <- nrow(data$x_train)
  tau_sample <- stats::rgamma(n = 1,shape = 0.5*n_+data$a_tau,rate = 0.5*crossprod((data$y_train-y_train_hat))+data$d_tau)

  return(tau_sample)

}


