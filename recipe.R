requiredPackages <- c("devtools", "shapr", "ggplot2", "tidyverse", "ggbeeswarm", "xgboost", "combinat")
for (package in requiredPackages) { 
  if (!requireNamespace(package, quietly = TRUE))
  #Installs packages if not yet installed    
  install.packages(package)
}

source("shapley_utils.R")

# Data Generation
set.seed(1234)

g_indirect_effect <- function(N){
    income <- 1 + 3*rnorm(N)
    race <- 1 + 5*rnorm(N)
    zip <- race + rnorm(N)/5
    y <- 0.5*income - 0.5*zip + rnorm(N)
    return(cbind(income, race, zip, y))
}

g_non_interventional <- function(N){
    income <- 1 + 3*rnorm(N)
    race <- 1 + 5*rnorm(N)
    favourite_show <- race + rnorm(N)/2
    y <- 0.5*income + 0.5*race + rnorm(N)
    return(cbind(income, race, favourite_show, y))
}

run_shapley <- function(generative_model, x_var, y_var, name){
    train <- generative_model(300)
    x_train <- train[, x_var]
    y_train <- train[, y_var]

    test <- generative_model(300)
    x_test <- test[, x_var]
    y_test <- test[, y_var]

    model <- xgboost(
      data = x_train,
      label = y_train,
      max.depth = 4,
      nround = 50,
      verbose = FALSE
    )
    
    return(list(
        x_var=x_var,
        x_train=x_train,
        x_test=x_test,
        y_train=y_train,
        y_test=y_test,
        model=model,
        name=name
        )
    )
}

# Other convenience functions
plot.prepare_data <- function(explanation, x_var){
  data.table::melt(
    explanation,
    measure.vars = x_var,
    variable.name = "feature",
    value.name = "shapley")
}

plot.global_shap <- function(explanation_table, title="title"){
  ggplot(data=explanation_table, aes(feature, shapley, color=factor(feature))) +
    ggtitle(title) +
    geom_beeswarm() +
    coord_flip()
}

approx_equal <- function(observed, target){
    difference <- observed-target
    return( abs(difference) < 0.00001 )
}

plot_distributions <- function(data){
  dataplot <- data[, c("income", "race", "zip", "f")]    
  means <- melt(dataplot[, lapply(.SD, mean)])
  
  ggplot(melt(dataplot),  aes(x=variable, y=value, fill=variable)) +
    geom_boxplot() +
    ggtitle("Sample distributions - features and target variable") +
    xlab("") +
    geom_text(data = means, aes(label = round(means$value, 2), y = value + 0.8))
  
  ggsave("images/21-variable-distributions.png")
}

plot_scenario <- function(p, hypothesis, alternative, hyp_coalitions, alt_coalitions, annotation, fname, title){  
  plot_data <- data.table(
      p = c(hypothesis, alternative), 
      coalition = c(rep(hyp_coalitions, length(hypothesis)), rep(alt_coalitions, length(alternative)))
  )

  ggplot(plot_data, aes(x=p, fill=coalition, color=coalition)) +
    ggtitle(title) +  
    geom_histogram(alpha=0.5, position="identity") +
    annotate("text", label=annotation, x = 4, y = 400)  
  
  ggsave(fname)
}

# Define outputs
demo.11 <- function(){
  plot.global_shap(data_non_interventional, "Shapley values - non-interventional example")

  ggsave("images/11-shapley-distribution-non-interventional.png")
}


demo.12 <- function(){
  data <- plot.prepare_data(explanation_symmetric$dt, run$x_var)
  plot.global_shap(data, "Shapley values - mediation example")
  
  ggsave("images/12-shapley-distribution-mediation.png")
}

demo.21 <- function(){
  data <- as.data.table(run$x_test)
  data[, pk := 1:.N]
  data[, f := predict(run$model, newdata = run$x_test, type = "response")]      
  
  plot_distributions(data)
}

demo.22 <- function(){
  plot(explanation_symmetric, index_x_test=test_id:test_id)
  
  ggsave("images/22-symmetric-test-example.png")
}

demo.23 <- function(){
  test_result <- t.test(p.010, p.011, paired=FALSE, alternative = "two.sided")
  
  annotation <- sprintf(
      "Two Sample t-test
      (true diff in means is 0)
      p-value=%s",
      round(test_result$p.value, 2)
  )
  
  plot_scenario(p, p.010, p.011, "{race}", "{race, zip}", annotation, "images/23-010-vs-011.png", "{race} scenario when target is zip")
}

demo.24 <- function(){
  test_result <- t.test(p.010, p.110, paired=FALSE, alternative = "two.sided")
  
  annotation <- sprintf(
      "Two Sample t-test
      (true diff in means is 0)
      p-value=%s",
      round(test_result$p.value, 2)
  )
  
  plot_scenario(p, p.010, p.110, "{race}", "{income, race}", annotation, "images/24-010-vs-110.png", "{race} scenario when target is income")
}

demo.31 <- function(){
  data <- plot.prepare_data(explanation, x_var)
  plot.global_shap(data, "Asymmetric Shapley values - mediation example")
  
  ggsave("images/31-asymmetric-global-distribution.png")
}

demo.32 <- function(){
  # Build a shap-like object - a bit hacky but I want to finish this article
  dt_asymmetric <- explanation
  dt_asymmetric[, none := p]

  explanation_asymmetric <- list(
      dt=dt_asymmetric,
      x_test=explanation_symmetric$x_test,
      p=explanation_symmetric$p)

  shapr:::plot.shapr(explanation_asymmetric, index_x_test = 100:100)
  
  ggsave("images/32-asymmetric-test-example.png")
}

# ***********************
# Execution starts here *
# ***********************

# A. Global Shapley
# Non interventional example
x_var <- c("income", "race", "favourite_show")
y_var <- "y"

run <- run_shapley(g_non_interventional, x_var, y_var, "non_interventional")

p <- mean(run$y_train)

print(p)

explainer_symmetric <- shapr(run$x_train, run$model)

explanation_symmetric <- explain(
  run$x_test,
  approach = "gaussian",
  explainer = explainer_symmetric,
  prediction_zero = p,
  seed = 2020
)

data_non_interventional <- plot.prepare_data(explanation_symmetric$dt, run$x_var)

# Mediation example
x_var <- c("income", "race", "zip")
y_var <- "y"

# Main run used throughout
run <- run_shapley(g_indirect_effect, x_var, y_var, "mediation")

# Make sure the run corresponds to the data used for the blog
stopifnot(
    approx_equal(run$y_test[[1]], -6.704256)
)

stopifnot(
    approx_equal(tail(run$y_test, 1), -3.973026)
)

p <- mean(run$y_train)

explainer_symmetric <- shapr(run$x_train, run$model)

explanation_symmetric <- explain(
  run$x_test,
  approach = "gaussian",
  explainer = explainer_symmetric,
  prediction_zero = p,
  seed = 2020
)

# B. Scenarios
# Found after manual search
test_id <- 100

# Scenarios
explainer_symmetric$feature_labels <- x_var

# Add gaussian items (covariance matrix, expectation)
explainer_symmetric <- gaussian_inputs(explanation_symmetric$x_test, explainer_symmetric, p)

combination_table <- prepare_data_gaussian(explainer_symmetric)

cnms <- colnames(explainer_symmetric$x_test)
data.table::setkeyv(combination_table, c("id", "id_combination"))

combination_table[, p_hat := predict_model(explainer_symmetric$model, newdata = .SD), .SDcols = cnms]

combination_table[id_combination == 1, p_hat := p]

p_all <- predict_model(explainer_symmetric$model, newdata = explainer_symmetric$x_test)
combination_table[id_combination == max(id_combination), p_hat := p_all[id]]

dt_res <- combination_table[, .(k = sum((p_hat * w) / sum(w))), .(id, id_combination)]
data.table::setkeyv(dt_res, c("id", "id_combination"))

dt_mat <- data.table::dcast(dt_res, id_combination ~ id, value.var = "k")
dt_mat[, id_combination := NULL]

kshap <-  t(explainer_symmetric$W %*% as.matrix(dt_mat))
dt_kshap <- data.table::as.data.table(kshap)

colnames(dt_kshap) <- c("none", cnms)

# Get distributions
p.all <- combination_table[(id==test_id)&(id_combination==8), p_hat]
p.110 <- combination_table[(id==test_id)&(id_combination==5), p_hat]
p.010 <- combination_table[(id==test_id)&(id_combination==3), p_hat]
p.011 <- combination_table[(id==test_id)&(id_combination==7), p_hat]

# C. Asymmetric Shapley
ordering_race <- get_ordering("race", explainer_symmetric, x_var, dt_res)
ordering_income <- get_ordering("income", explainer_symmetric, x_var, dt_res)
ordering_zip <- get_ordering("zip", explainer_symmetric, x_var, dt_res)

phi_zip <- ordering_zip[, list(zip=sum(weighted_v_difference)), by=test_instance_id]
phi_income <- ordering_income[, list(income=sum(weighted_v_difference)), by=test_instance_id]
phi_race <- ordering_race[, list(race=sum(weighted_v_difference)), by=test_instance_id]

explanation_colnames <- c("test_instance_id", x_var)
explanation <- cbind(phi_income, phi_race, phi_zip)[, ..explanation_colnames]

# Control outputs to generate
demo.11()
demo.12()
demo.21()
demo.22()
demo.23()
demo.24()
demo.31()
demo.32()