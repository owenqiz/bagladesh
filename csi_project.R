##################### PART 1 COMPUTATIONS ################################
##
##
##########################################################################
# Computation for CSI
# Estimation based on truncated Empirical Distribution

# Execute this line if you have the data
library(tidyverse)
library(ggplot2)
load('dataset.Rda')
# remove ONE 0m depth tubewell from tid = 245, Tangail S. rowname = 2931
NationalSurveyData <- NationalSurveyData %>% filter(depth > 0)
thana_obs <- split(NationalSurveyData, NationalSurveyData$tid)
dist_obs <- split(NationalSurveyData, NationalSurveyData$did)

propensity_score <-function(x, y_vector, L, trunc = FALSE){
  # x is argument, i.e. oberservation to be convert
  # y_vector is a vector of ALL observations in the sub_region
  
  y_truncate <- y_vector[y_vector > L]

  # case that ALL observation less than the threshold
  if(length(y_truncate) == 0)
    return(rep(0, length(y_vector)))
  
  # when num of y_trancate == 1, then 2/(n-M-1) = infinity
  if(length(y_truncate) == 1 & trunc)
    return(max(0,1-L/x))
    
  # the base number of Eq. (2.1)
  C_nx <- max(0,1-L/x)

  # csi_new
  coeff <- ifelse(trunc, 2/length(y_truncate), 2/length(y_vector))

  y_sum <- ifelse(trunc, sum(pmin(x,y_truncate)/(x + y_truncate)), 
                          sum(pmin(x,y_vector)/(x + y_vector)))
  
  # return value of Eq. (2.1)
  return(C_nx^(coeff*y_sum))
}

csi <- function(y_vector, L = L, trunc = FALSE){
  score <- sapply(y_vector, propensity_score, L = L, y_vector = y_vector, trunc = trunc)
  
  if(all(score == 0))
    return(0)
  
  if(length(score) == 1)
    returnValue(score)

  return(mean(score))
}

logit <- function(p){return(log(p/(1-p)))}

logit_fixed <- function(p, threshold){
  if (p <= threshold) {
    return(log(threshold/(1-threshold)))
  } else if (p >= (1 - threshold)) {
    return(log((1-threshold)/threshold))
  } else {
    return(log(p/(1-p)))
  }
}

csi_adj <- function(y_vector, L = L, threshold, lgt_fixed = TRUE, trunc = FALSE){
  S <- sapply(y_vector, propensity_score, y_vector, L = L, trunc = trunc)

  # for thana has all As < L, return 0, otherwise will cause error for parameter estimation
  if(all(S == 0))
    return(0)

  n <- length(S)
  S <- S[y_vector > L]
  m <- length(S)
  
  if(lgt_fixed){
    L <- mean(sapply(S, logit_fixed, threshold = threshold))
  } else {
    L <- mean(sapply(S, logit))
  }

  return(m/(1+exp(-L))/n)
}

# district adjusted csi
district_csi <- function(district, adj = FALSE, Ls = L){
  n_total <- sum(district$size)
  n_i <- district$size
  
  # request from Dr.Sen
  threshold <- ifelse(Ls == 10, .01, .002)
  
  if(adj){
    csi_val <- district$csi_adj
  } else {
    csi_val <- district$csi
  }
  
  if(is.na(threshold)){
    L <- logit(csi_val)
  } else{
    L <- sapply(csi_val, logit_fixed, threshold = threshold)
  }
  
  T_w <- sum(n_i * L/n_total)
  return(1/(1+exp(-T_w)))
}

thana_csi <- function(data_list, L){
  n <- length(data_list)
  threshold <- ifelse(L == 10, .01, .002)
  csi_vector <- csi_adj_vector <- csi_max <- csi_min <- size <- AAs <- ALd <- rep(NA, n)
  vid <- did <- tid <- div <- dis <- tha <- rep(NA, n)

  for(i in 1:n){
  y_obs <- data_list[[i]]$as
  depth <- data_list[[i]]$depth
  
  vid[i] <- unique(data_list[[i]]$vid)
  did[i] <- unique(data_list[[i]]$did)
  tid[i] <- unique(data_list[[i]]$tid)
  div[i] <- unique(data_list[[i]]$division)
  dis[i] <- unique(data_list[[i]]$district)
  tha[i] <- unique(data_list[[i]]$thana)
  
  csi_vector[i] <- csi(y_obs, L)
  csi_adj_vector[i] <- csi_adj(y_obs, L, threshold)
  csi_max[i] <- length(y_obs[y_obs > L])/length(y_obs)
  csi_min[i] <- sum(1 - L/y_obs[y_obs > L])/length(y_obs)
  size[i] <- length(y_obs)
  AAs[i] <- mean(y_obs)
  ALd[i] <- mean(log(depth))
  }
  df <- tibble(vid, did, tid, 
               division = div, district = dis, thana = tha, size, 
               csimax = csi_max, csi = csi_vector, csi_adj = csi_adj_vector, csimin = csi_min, 
               ALd, AAs)
  return(df)
}

# thana_level csi table
thana_10 <- thana_csi(thana_obs, L = 10)
thana_50 <- thana_csi(thana_obs, L = 50)

# export computation results
write.csv(thana_10, file = 'thana_10.csv')
write.csv(thana_50, file = 'thana_50.csv')

# district_level csi table
dist_csi <- function(data_list, thana_table, L){
  n <- length(data_list)
  csi_vector <- csi_adj_vector <- csi_max <- csi_min <- size <- AAs <- ALd <- rep(NA, n)
  vid <- did <- div <- dis <- rep(NA, n)
  
  tmp_list <- split(thana_table, thana_table$did)

  for(i in 1:n){
    df <- data_list[[i]]
    y_obs <- df$as
    vid[i] <- unique(df$vid)
    did[i] <- unique(df$did)
    div[i] <- unique(df$division)
    dis[i] <- unique(df$district)
    size[i] <- nrow(df)
    csi_max[i] <- length(y_obs[y_obs > L])/length(y_obs)
    csi_min[i] <- sum(1 - L/y_obs[y_obs > L])/length(y_obs)
    csi_vector[i] <- district_csi(tmp_list[[i]], adj = FALSE, L)
    csi_adj_vector[i] <- district_csi(tmp_list[[i]], adj = TRUE, L)
    ALd[i] <- mean(log(df$depth))
    AAs[i] <- mean(y_obs)
  }
  
  df <- tibble(vid, did, 
               division = div, district = dis, size, 
               csimax = csi_max, csi = csi_vector, csi_adj = csi_adj_vector, csimin = csi_min, 
               ALd, AAs)
  return(df)
}

dist_10 <- dist_csi(dist_obs, thana_10, L = 10)
dist_50 <- dist_csi(dist_obs, thana_50, L = 50)

write.csv(dist_10, file = 'dist_10.csv')
write.csv(dist_50, file = 'dist_50.csv')

##################### PART 1 END ################################


##################### PART 2 LaTeX Tables ################################
##
##
##########################################################################
library(tidyverse)
library(xtable)

# you may load the results from PART 1, ignore warning message is OK
# thana_level CSI table
thana_10 <- read_csv("thana_10.csv", col_types = cols(X1 = col_skip(), 
                                                      did = col_integer(), size = col_integer(), 
                                                      tid = col_integer(), vid = col_integer()))

thana_50 <- read_csv("thana_50.csv", col_types = cols(X1 = col_skip(), 
                                                      did = col_integer(), size = col_integer(), 
                                                      tid = col_integer(), vid = col_integer()))


table_all <- bind_cols(thana_10, thana_50)

t1 <- table_all %>% select(c(1:11, 21:26))

colnames(t1) <- c('vid', 'did', 'tid', 'division', 'district', 'thana',
                  'n', 'csimax10', 'csi10', 'csiadj10', 'csimin10',
                  'csimax50', 'csi50', 'csiadj50', 'csimin50', 'Ald', 'AAs') 

digit_control <- c(rep(0, 4), rep(5, 8), rep(3, 2))

# print entire table
# xt <- t1[,-c(1:4)]
# xt <- xtable(xt, method = 'compact', digits = digit_control, caption = 'Thana Level CSI')
# print(xt, tabular.environment = "longtable", include.rownames = FALSE, booktabs = TRUE, caption.placement = 'top')

# print table by division (vid) for ordering
# division order is {6,5,9,7,4,3,8,1,2}
print_thana <- function(vid){
  id <- vid
  mx <-  subset(t1, t1$vid == id)
  xt <- mx[,-c(1:4)]
  xt <- xtable(xt, method = 'compact', digits = digit_control)
  
  # adding midrules for each district
  idx <- xt %>% group_by(district) %>% summarise(ids = n())
  idx <- idx$ids[-nrow(idx)]
  idx <- cumsum(idx)
  midrs <- rep('\\midrule \n', length(idx))
  
  print(xt, tabular.environment = "longtable", 
        include.rownames = FALSE, booktabs = TRUE,
        add.to.row = list(pos = as.list(idx), command = midrs))
}

# district_level CSI table

dist_10 <- read_csv("dist_10.csv", col_types = cols(X1 = col_skip(), 
                                                    did = col_integer(), size = col_integer(), 
                                                    vid = col_integer()))

dist_50 <- read_csv("dist_50.csv", col_types = cols(X1 = col_skip(), 
                                                    did = col_integer(), size = col_integer(), 
                                                    vid = col_integer()))

table_all <- bind_cols(dist_10, dist_50)
t1 <- table_all %>% select(c(1:9, 17:22))

colnames(t1) <- c('vid', 'did', 'division', 'district', 
                  'n', 'csimax10', 'csi10', 'csiadj10', 'csimin10',
                  'csimax50', 'csi50', 'csiadj50', 'csimin50', 'Ald', 'AAs') 

digit_control <- c(rep(0, 4), rep(5, 8), rep(3, 2))

# ord <- c(6,5,9,7,4,3,8,1,2)
# almost same except for line 230
print_district <- function(vid){
  id <- vid
  mx <-  subset(t1, t1$vid == id)
  xt <- mx[,-c(1:2)]
  xt <- xtable(xt, method = 'compact', digits = digit_control)
  print(xt, tabular.environment = "longtable", include.rownames = FALSE, booktabs = TRUE)
}

##################### PART 2 END ################################


##################### PART 3 PLOTTING MAPS ################################
##
##
##########################################################################
# alternative way better legend for now
library(tidyverse)
library(raster)
library(rgeos)


# loading maps from GADM.org
div_map <- readRDS(gzcon(url('https://biogeo.ucdavis.edu/data/gadm2.8/rds/BGD_adm1.rds')))
tha_map <- readRDS(gzcon(url('https://biogeo.ucdavis.edu/data/gadm2.8/rds/BGD_adm3.rds')))

# border information
load('border_chittagong_comilla.Rda')
load('border_dhaka_mymensingh.Rda')

# you may load the results from PART 1, ignore warning message is OK
thana_10 <- read_csv("thana_10.csv", col_types = cols(X1 = col_skip(), 
                                                      did = col_integer(), size = col_integer(), 
                                                      tid = col_integer(), vid = col_integer()))

thana_50 <- read_csv("thana_50.csv", col_types = cols(X1 = col_skip(), 
                                                      did = col_integer(), size = col_integer(), 
                                                      tid = col_integer(), vid = col_integer()))


tha_csi <- tibble(id = thana_10$tid, csi10 = thana_10$csi_adj, csi50 = thana_50$csi_adj)

thana_df <- fortify(tha_map)
thana_df$id <- as.numeric(thana_df$id)

t_csi <- left_join(thana_df, tha_csi,by = 'id')

### BLACK BORDER
col_pal <- c('white', 'ivory2', 'violet', 'lightblue', 'navyblue', 'yellow2', 'darkorange2', 'lightgreen', 'darkgreen','rosybrown2' ,'brown','red3')
brks <- seq(0, 1, by = .1)
size <- .65

p_10 <- ggplot() +
  geom_map(data = t_csi, map = t_csi, aes(x = long, y = lat, fill = csi10, map_id = id, color = 'NA')) +
  coord_map() + scale_fill_gradientn(colours = col_pal, breaks = brks, limits=c(0,1)) + 
  scale_color_manual('', values = 'grey30', labels = 'No data') +
  geom_polygon(data = div_map, aes(x = long, y = lat, group = group), color="red",  size = size, alpha = 0) +
  geom_path(data = chi_border, aes(x = long, y = lat), color = 'red', size = size) + 
  geom_path(data = dha_border, aes(x = long, y = lat), color = 'red', size = size) + 
  guides(fill = guide_colourbar(barwidth = 1, barheight = 18), 
         color = guide_legend(override.aes = list(labs(fill = ''), fill = 'grey50'))) +
  ggtitle('Thana Level CSI under L = 10') + xlab("Longitude") + ylab("Latitude") +  labs(fill = "CSI") +
  theme(plot.title = element_text(hjust = 0.5))

p_50 <- ggplot() +
  geom_map(data = t_csi, map = t_csi, aes(x = long, y = lat, fill = csi50, map_id = id, color = 'NA')) +
  coord_map() + scale_fill_gradientn(colours = col_pal, breaks = brks, limits=c(0,1)) + 
  scale_color_manual('', values = 'grey30', labels = 'No data') +
  geom_polygon(data = div_map, aes(x = long, y = lat, group = group), color="red",  size = size, alpha = 0) +
  geom_path(data = chi_border, aes(x = long, y = lat), color = 'red', size = size) + 
  geom_path(data = dha_border, aes(x = long, y = lat), color = 'red', size = size) + 
  guides(fill = guide_colourbar(barwidth = 1, barheight = 18), 
         color = guide_legend(override.aes = list(labs(fill = ''), fill = 'grey50'))) +
  ggtitle('Thana Level CSI under L = 50') + xlab("Longitude") + ylab("Latitude") +  labs(fill = "CSI") +
  theme(plot.title = element_text(hjust = 0.5))

pdf("csi_map.pdf", onefile = TRUE)
print(p_10)
print(p_50)
dev.off()