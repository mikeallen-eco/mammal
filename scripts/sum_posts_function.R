# Create a function to summarize posteriors
sum_posts = function(x, method = NULL, round = NULL){
  sum_posts_list = aggregate(x, by = list(rep("A", nrow(x))), FUN = function(x) quantile(x, c(0.025, 0.10, 0.5, 0.90, 0.975))) %>%
    select(-Group.1) %>%
    c()
  
  sum_posts_df = do.call("rbind", sum_posts_list) %>% 
    as.data.frame() %>%
    rename(q2.5 = 1, q10 = 2, q50 = 3, q90 = 4, q97.5 = 5) %>%
    mutate(n = 1:length(sum_posts_list),
           Method = method,
           Round = round)
  
  return(sum_posts_df)
}