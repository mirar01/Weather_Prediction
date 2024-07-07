# Your Name (s1915415, mirar01)
# Place analysis code that may take too long to run every time the report.Rmd
# document is run.
# Run the code in this file with
#   source("analysis.R")
# in a fresh R session to ensure consistency of the results.

# Load function definitions
source("functions.R")

permutationTest <- function(data, nPermutation = 10000)
{
  #The first part of the function separates months into winter and summer seasons
  #and creates three data sets: one with summer, one with winter, and one with both
  #the seasons' values

  Winter_values <- data %>%
    mutate(Season = ifelse(Month==1|Month==2|Month==3|Month==10|Month==11|Month==12, "Winter", "Summer")) %>%
    filter(Element %in% c("PRCP")) %>%
    filter(Season %in% c("Winter")) %>%
    select(Name, Season, Value)

  Summer_values <- data %>%
    mutate(Season = ifelse(Month==1|Month==2|Month==3|Month==10|Month==11|Month==12, "Winter", "Summer")) %>%
    filter(Season %in% c("Summer")) %>%
    filter(Element %in% c("PRCP")) %>%
    select(Name, Season, Value)

  All_values <-data %>%
    mutate(Season = ifelse(Month==1|Month==2|Month==3|Month==10|Month==11|Month==12, "Winter", "Summer")) %>%
    filter(Element %in% c("PRCP")) %>%
    select(Name, Season, Value)

  #Next, empty vectors and lists are defined to later store results in. The reason for
  #using a list rather than a data frame is to allow for the function to produce both
  #the desired values and histograms, even though it makes the code bulkier

  T_Obs_a <- c(0)
  pVal_a <- c(0)
  T_Perm_a <- c(0)
  sd_a <- c(0)

  T_Obs_b <- c(0)
  pVal_b <- c(0)
  T_Perm_b <- c(0)
  sd_b <- c(0)

  rain_dist_w <- c()
  rain_prob_w <- c()
  rain_dist_s <- c()
  rain_prob_s <- c()

  Names <- unique(data$Name)

  plots <- vector('list', length(Names))

  #A for loop is used to separate the stations and get different values for each

  for (i in seq_along(Names)){

    #filtering out desired values per station

    T_W_S <- Winter_values  %>%
      filter(Name %in% c(Names[i]))
    T_W_vals <- T_W_S$Value

    #creating binary values depending on whether or not there was any rainfall

    W_rainy_days <- T_W_S %>%
      mutate(Value = ifelse(Value>0, 1,0))
    W_rd_vec <- W_rainy_days$Value


    T_S_S <- Summer_values %>%
      filter(Name %in% c(Names[i]))
    T_S_vals <- T_S_S$Value

    S_rainy_days <- T_S_S %>%
      mutate(Value = ifelse(Value>0, 1,0))
    S_rd_vec <- S_rainy_days$Value


    T_A_S <- All_values %>%
      filter(Name %in% c(Names[i]))
    all_Ts <- T_A_S$Value

    A_rainy_days <- T_A_S %>%
      mutate(Value = ifelse(Value>0, 1,0))
    A_rd_vec <- A_rainy_days$Value

    #Calculating the observed Test statistics

    meanT_W_vals <- mean(T_W_vals)
    meanT_S_vals <- mean(T_S_vals)
    rain_dist_w[i] <- c( meanT_W_vals)
    rain_dist_s[i] <- c( meanT_S_vals)
    T_Obs_a[i] <- abs(meanT_W_vals - meanT_S_vals)

    rain_prob_w[i] <- c(mean(W_rd_vec ))
    rain_prob_s[i] <- c(mean(S_rd_vec ))
    T_Obs_b[i] <- abs(mean(W_rd_vec )- mean(S_rd_vec ))

    # Obtaining permuted test statistics by shuffling the combined prcp values

    T_Perm_a <- replicate(nPermutation, expr = {
      shuffle <- sample(all_Ts, length(all_Ts), replace = FALSE)
      W_Sample <- shuffle[1:length(T_W_vals)]
      S_Sample <- shuffle[length(T_S_vals):length(all_Ts)]

      abs(mean(W_Sample)-mean(S_Sample))

    })

    # Plotting the histograms for the first MC permutation test

    plot_data <- data.frame(Name = replicate(nPermutation, Names[i]),
                            T_Stat_Observed = replicate(nPermutation, T_Obs_a[i]),
                            T_Permutation = T_Perm_a)


    plots[[i]] <- ggplot(data = plot_data,
             aes(x= T_Permutation)) +
             geom_histogram()+
             ggtitle(Names[i])+
             geom_vline(aes(xintercept=unique(T_Stat_Observed)),
                                                    color="blue", linetype="dashed", size=1)




    T_Perm_b <- replicate(nPermutation, expr = {
      shuffle <- sample(A_rd_vec, length(A_rd_vec), replace = FALSE)
      W_Sample <- shuffle[1:length(W_rd_vec)]
      S_Sample <- shuffle[length(S_rd_vec):length(A_rd_vec)]

      abs(mean(W_Sample)-mean(S_Sample))

    })

    #Calculating the p values and standard deviations

    pVal_a[i] <- sum(abs(T_Perm_a)>=abs(T_Obs_a[i])) / nPermutation
    sd_a[i]<- sqrt(((pVal_a[i])*(1-pVal_a[i]))/ nPermutation)

    pVal_b[i] <- sum(abs(T_Perm_b)>=abs(T_Obs_b[i])) / nPermutation
    sd_b[i]<- sqrt(((pVal_b[i])*(1-pVal_b[i]))/ nPermutation)

  }


  list(winter_rainfall_distribution = rain_dist_w,summer_rainfall_distribution = rain_dist_s,
       Test_Statistic_a=T_Obs_a, p_value_a=pVal_a,  Standard_Deviation_a= sd_a,
       winter_rainfall_probability = rain_prob_w,summer_rainfall_probability = rain_prob_s,
       Test_Statistic_b=T_Obs_b, p_value_b=pVal_b, Standard_Deviation_b = sd_b,
       Station = Names, plots)
}



mcp_tests<-permutationTest(ghcnd)
saveRDS(mcp_tests, file = "data/mcp_tests.rds")

