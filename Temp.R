# A clean slate

rm(list=ls())


# Set working directory to required folder

setwd("/home/ubuntu/EDW DATA - INPUT/17th_OCT_AWS")

# Check list of files in the directory

list.files()


# Import all the required libraries

library(dplyr)
library(data.table)
library(reshape2)
library(coop)
library(doParallel)
library(mailR)
library(zoo)
library(inline)
library(doSNOW)



# Decide states and premises to loop over

#state = c('IL','WA' , 'GA','FL','CA_S', 'CA_N', 'TX', 'NJ' , 'NY')
#state = c('CA_S', 'CA_N', 'NY')
state = c('TX','NJ') # state
prem = c('O','F') # premise



#####################Initiating All Functions Before Loop Begins######################


## Function to retain 0 at the start for store codes
## Needed till we deal with csv files

zero_format <- function(x){
  if(nchar(x) < 7){
    x = paste0(paste0(rep(0,7-nchar(x)),collapse = ""),x)
  }else{
    x
  }
}

## Haversine Function for ON & OFF Premise

gcd.hf <- function(haystack, needle, threshold){
  
  R <- 6371 # Earth mean radius [km]
  
  dlat = haystack[,3] - needle[1]
  dlon = haystack[,4] - needle[2]
  
  a <- sin(dlat/2)^2 + cos(needle[1]) * cos(haystack[,3]) * sin(dlon/2)^2
  c <- 2 * atan2(sqrt(a),sqrt(1-a))
  
  d = R * c
  
  return(structure(list(dist = min(d[d>0]) , 
                        count = sum(d>0 & d<= threshold))))
  #min_st = haystack[which.min(d[d>0]),1])))
}

new_gcd <- compiler::cmpfun(gcd.hf)



#Adding missing demographics information

demo_fix <- function(haystack, needle) {
  
  R <- 6371 # Earth mean radius [km]
  dlat = haystack[,2] - needle[1]
  dlon = haystack[,3] - needle[2]
  a <- sin(dlat/2)^2 + cos(needle[1]) * cos(haystack[,2]) * sin(dlon/2)^2
  c <- 2 * atan2(sqrt(a),sqrt(1-a))
  d = R * c
  return(haystack[which.min(d),4:ncol(demo_fix_miss_ids)])
  
}


## Fill Demographics


fill_demographics <- function(naDf, fillDf, mergeCols, fillCols) {
  fillB <- do.call(paste, c(fillDf[, mergeCols, drop = FALSE], sep="\r"))
  naB <- do.call(paste, c(naDf[, mergeCols, drop = FALSE], sep="\r"))
  m <- match(naB, fillB)
  for(col in fillCols) {
    fix <- which(is.na(naDf[,col]))
    naDf[fix, col] <- fillDf[m[fix],col]
  }
  naDf
}


## data normalization

normalise <- function(x) {
  return ((x - min(x)) / 
            (max(x) - min(x))) }

## Haversine Function for Store_Sim


store_dist <- function(long1, lat1, long2, lat2) {
  R <- 6371 # Earth mean radius [km]
  delta.long <- (long2 - long1)
  delta.lat <- (lat2 - lat1)
  a <- sin(delta.lat/2)^2 + cos(lat1) * cos(lat2) * sin(delta.long/2)^2
  c <- 2 * asin(min(1,sqrt(a)))
  d = R * c
  return(d) # Distance in km
}

## Function for Combining predicted matrix and similar stores


combine_custom_j <- function(LL1, LL2) { 
  
  bx <- rbind(LL1$bx, LL2$bx)  #sapply(LL, function(x) x$bx)
  dfs <- rbindlist(list(LL1$df, LL2$df))  #lapply(LL, function(x) x$df)
  return(list(bx = bx, df = dfs))
} 

combine_custom_i <- function(LL1, LL2) {
  
  bx <- c(LL1$bx,LL2$bx) #do.call("rbind", lapply(LL, function(x) x$bx))
  dfs <- rbindlist(list(LL1$df, LL2$df))
  return(list(bx = bx, df = dfs))
}


## Rec-Sys Function for Parallel Processing


rec_sys <- function(m,v){  
  
  stores1 = which(!is.na(prod_mat2[-m,v]),arr.ind = T)
  stores1 = as.integer(names(stores1))
  k=20
  
  if(!(is.na(length(stores1)) | length(stores1) == 0)){
    
    
    thres_stores <- which(sim >= lower & sim <= higher)
    com_stores <- thres_stores[thres_stores %in% stores1]
    rm(thres_stores)
    
    
    if(length(stores1) == 1 | length(com_stores) == 1){
      if(!identical(com_stores , integer(0)) & length(com_stores) == 1){
        stores1 = com_stores
      }
      bx = prod_mat2[stores1,v]
      df = data.frame(m,v,stores1,round(sim[stores1],3))
      
      
      return(list(bx = bx , df = df))
      rm(df)
      gc()
      
    } else {
      
      if(identical(com_stores,NULL) | identical(com_stores, integer(0))){
        
        # knn
        if(length(stores1) < k){
          k = length(stores1)
        }
        
        
        near.ind = as.numeric(names(sort(sim[1,stores1],decreasing = T)[1:k])[!names(sort(sim[1,stores1],decreasing = T)[1:k]) %in% NA])
        
        bx = crossprod(as.vector(t(prod_mat2[near.ind,v])),sim[near.ind]) / sum(sim[near.ind])
        df = data.frame(m,v,near.ind,round(sort(sim[stores1],decreasing = T)[1:k],3))
        
      }
      else{
        rm(stores1)
        com_stores = as.numeric(names(sort(sim[1,com_stores] , decreasing = T)[1:100])[!names(sort(sim[1,com_stores] , decreasing = T)[1:100]) %in% NA])
        bx = crossprod(as.vector(t(prod_mat2[com_stores,v])), sim[com_stores]) / sum(sim[com_stores])
        df = data.frame(m,v,com_stores,round(sim[com_stores],3))
        
        
      }
      
      return(list(bx = bx , df = df))
      rm(df)
      gc()
      
      
    }
  }else{
    
    bx = NaN
    prod_col = which(!is.na(dat1[,paste0(colnames(prod_mat2)[v],'_L365')]))
    req_store = as.numeric(names(sort(sim[1,prod_col] , decreasing = T)[1:100])[!names(sort(sim[1,prod_col] , decreasing = T)[1:100]) %in% NA])
    return(list(bx= bx,df= data.frame(m,v,req_store,sim[req_store])))
    rm(prod_col , req_store)
    gc()
    
    
  }  
}

##################################Running it for all States########################

# Looping over states

for(q in 1:length(state)){
  
  print(state[q])
  
  ## read both on and off data
  off_premise_rec <-read.csv(paste0(state[q],"_OFF.csv") , stringsAsFactors = F)
  on_premise_rec <-read.csv(paste0(state[q], "_ON.csv"), stringsAsFactors = F)
  
  ## Remove stores with VIP/W
  off_premise_rec = off_premise_rec[-grep('VIP|W', off_premise_rec$TDLINX_STORE_CD),]
  
  
  #Convert degrees to radians
  # Calculates the geodesic distance between two points specified by radian latitude/longitude using the
  # Haversine formula (hf)
  
  on_data <- on_premise_rec[,c("TDLINX_STORE_CD","PREMISE_TYPE_DSC","LATITUDE","LONGITUDE")]
  off_data <- off_premise_rec[,c("TDLINX_STORE_CD","PREMISE_TYPE_DSC","LATITUDE","LONGITUDE")]
  dist_col = c('LATITUDE','LONGITUDE') 
  
  
  on_data[,dist_col][] <- lapply(on_data[,dist_col] , function(x) return (x*pi/180))
  off_data[,dist_col][] <- lapply(off_data[,dist_col] , function(x) return (x*pi/180))
  
  
  distance_data = rbind(on_data,off_data)
 
  # Initiate parallel backend to calculate Haversine Distance
  
  cl = makeCluster(detectCores()-1)
  registerDoParallel(cl)
  
  
  x1 <- 
    foreach(row = iter(on_data[,3:4] , by = 'row') ,
            .combine = 'rbind')  %dopar%
            {
              new_gcd(on_data, as.numeric(row) , 0.5 )
            }
  
  x2 <- 
    foreach(row = iter(on_data[,3:4] , by = 'row') ,
            .combine = 'rbind')  %dopar%
            {
              new_gcd(off_data, as.numeric(row) , 0.5 )
            }
  
  x <- cbind(x1,x2)
  
  
  y1 <- foreach(row = iter(off_data[,3:4] , by = 'row'),
                .combine = 'rbind') %dopar%
                {
                  new_gcd(on_data , as.numeric(row), 0.5)
                }
  
  
  y2 <- foreach(row = iter(off_data[,3:4] , by = 'row'),
                .combine = 'rbind') %dopar%
                {
                  new_gcd(off_data , as.numeric(row), 0.5)
                }
  
  
  stopCluster(cl)
  
  y <- cbind(y1,y2)
  base = rbind(x,y)
  
  
  # Combining Haversine Distance for ON and OFF Premise Data
  
  distance_data = cbind(distance_data,base)
  colnames(distance_data)[5:8] <- c('Min Distance On-Premise' , 'Store Count On-Premise',
                                    'Min Distance Off-Premise', 'Store Count Off-Premise')
  
  
  distance_data <- distance_data[,c(1,5,6,7,8)]
  rm(off_data,on_data,x,y,x1,x2,y1,y2)
  gc()
  
  print("Distance Complete")
  
  
  ###############################Data Clean Up and demographics fill###########################
  
  dem_lst<-c("TDLINX_STORE_CD",
             'LATITUDE',
             'LONGITUDE',
             "WHITE_POP_PCT",
             "BLACK_POP_PCT",
             "ASIAN_POP_PCT",
             "HISP_POP_PCT",
             "LEAST_ACC_HISP_POP_PCT",
             "BI_CULTURAL_HISP_POP_PCT",
             "MOST_ACC_HISP_POP_PCT",
             "OTHER_POP_PCT",
             "INCOME_POP_L10_PCT",
             "INCOME_POP_1020_PCT",
             "INCOME_POP_2030_PCT",
             "INCOME_POP_3040_PCT",
             "INCOME_POP_4050_PCT",
             "INCOME_POP_5075_PCT",
             "INCOME_POP_75100_PCT",
             "INCOME_POP_G100_PCT",
             "MALE_POP_2124_PCT",
             "MALE_POP_2534_PCT",
             "MALE_POP_3544_PCT",
             "MALE_POP_4554_PCT",
             "MALE_POP_5564_PCT",
             "MALE_POP_6574_PCT",
             "MALE_POP_75_PCT",
             "FEMALE_POP_2124_PCT",
             "FEMALE_POP_2534_PCT",
             "FEMALE_POP_3544_PCT",
             "FEMALE_POP_4554_PCT",
             "FEMALE_POP_5564_PCT",
             "FEMALE_POP_6574_PCT",
             "FEMALE_POP_75_PCT")
  
  demo_fix_miss<-rbind(off_premise_rec[,c(dem_lst)],on_premise_rec[,c(dem_lst)])
  demo_non_miss<-demo_fix_miss[rowSums(is.na(demo_fix_miss))!=30,]
  demo_fix_miss_ids<- demo_fix_miss[rowSums(is.na(demo_fix_miss))==30,]
  
  
  
  ## Initializing Parallel Backend for Demographic Fill
  
  cl = makeCluster(detectCores()-1)
  registerDoParallel(cl)
  
  
  base = foreach(i = 1:nrow(demo_fix_miss_ids) , .combine = 'rbind') %dopar%{
    
    needle = as.numeric(demo_fix_miss_ids[i,c('LATITUDE','LONGITUDE')])
    cbind(demo_fix_miss_ids[i,1:3] , demo_fix(demo_non_miss , needle))
    
  }
  
  stopCluster(cl)
  
  ## Reading Prod Info Data
  
  prod_info <- read.csv(paste0(state[q],'_EXTRACTS.csv'), stringsAsFactors = F)
  prod_info$MASTER_PKG_SKU_CD <- as.character(prod_info$MASTER_PKG_SKU_CD)
  prod_info$TDLINX_STORE_CD[] = sapply(prod_info$TDLINX_STORE_CD , zero_format)
  
  # Removing not-required objects from memory
  
  rm(demo_fix_miss_ids,demo_fix_miss,demo_non_miss)
  gc()
  
  
  # Looping over premises
  
  for(z in 1:length(prem)){
    
    tima1 = Sys.time()
    
    if(prem[z] == 'O'){

      filtered_recom = on_premise_rec
      rm(on_premise_rec)
      
    } else {
    
      filtered_recom = off_premise_rec
      rm(off_premise_rec)
      
    }
    
    gc()
    
    print(prem[z])
    
    # Send process notifications via email to Slack and/or cbrand id
    
    sender <- "abhijeet.singh@cbrands.com"
    #recipients <- c("u3s3y1q7u1j9x0f9@informationdelivery.slack.com") #slack id
    recipients <- c("abhijeet.singh@cbrands.com", "abhishek.shetty@fractalanalytics.com") #personal mails
    send.mail(from = sender,
              to = recipients,
              subject = "Code Run Status",
              body = paste(state[q],prem[z],sep="-"),
              smtp = list(host.name = "smtp.office365.com", port = 587,
                          user.name = "abhijeet.singh@cbrands.com",           
                          passwd = "Leckie0616", tls = TRUE),
              authenticate = TRUE)    
    
    
    
    filtered_recom[,dist_col][] = lapply(filtered_recom[,dist_col], function(x) return (x*pi/180))
    
    prod_list <- colnames(filtered_recom)[which(grepl("X", colnames(filtered_recom)) & 
                                                  grepl("[[:digit:]]", colnames(filtered_recom)))]
    
    prods = gsub('X',"" , prod_list)
    
    
    ## Add distance columns to CF data
    
    final = plyr::join(filtered_recom, distance_data , type = 'left' , match = 'all' , by = 'TDLINX_STORE_CD')
    colnames(final) = gsub(".y","",colnames(final))
    rm(distance_data)
    gc()
    
    fillcols = c("WHITE_POP_PCT",
                 "BLACK_POP_PCT",
                 "ASIAN_POP_PCT",
                 "HISP_POP_PCT",
                 "LEAST_ACC_HISP_POP_PCT",
                 "BI_CULTURAL_HISP_POP_PCT",
                 "MOST_ACC_HISP_POP_PCT",
                 "OTHER_POP_PCT",
                 "INCOME_POP_L10_PCT",
                 "INCOME_POP_1020_PCT",
                 "INCOME_POP_2030_PCT",
                 "INCOME_POP_3040_PCT",
                 "INCOME_POP_4050_PCT",
                 "INCOME_POP_5075_PCT",
                 "INCOME_POP_75100_PCT",
                 "INCOME_POP_G100_PCT",
                 "MALE_POP_2124_PCT",
                 "MALE_POP_2534_PCT",
                 "MALE_POP_3544_PCT",
                 "MALE_POP_4554_PCT",
                 "MALE_POP_5564_PCT",
                 "MALE_POP_6574_PCT",
                 "MALE_POP_75_PCT",
                 "FEMALE_POP_2124_PCT",
                 "FEMALE_POP_2534_PCT",
                 "FEMALE_POP_3544_PCT",
                 "FEMALE_POP_4554_PCT",
                 "FEMALE_POP_5564_PCT",
                 "FEMALE_POP_6574_PCT",
                 "FEMALE_POP_75_PCT")
    
    
    mergeCols = c("TDLINX_STORE_CD")
    
    
    
    final_dat_two = fill_demographics(naDf = final,fillDf = base,mergeCols = mergeCols,fillCols = fillcols)
    
    
    prod_index <- colnames(final_dat_two)[which(grepl("X", colnames(final_dat_two)) & 
                                                  grepl("[[:digit:]]", colnames(final_dat_two)))]
    
    
    rm(final)
    gc()
    
    print("Demographics Complete")
    
    if(prem[z] == 'F'){
    
      keep<-c("TDLINX_STORE_CD","N_STATE_WEEKLY_VOLUME","N_STATE_SQ_FT","B_SUBCHANNEL_BARLOCAL", 
              "B_SUBCHANNEL_AUTO","B_SUBCHANNEL_BOOKS","B_SUBCHANNEL_ELEC","B_SUBCHANNEL_CRAFT", 
              "B_SUBCHANNEL_DEPT","B_CSUBHANNEL_DISCDEPT","B_SUBCHANNEL_FAMILY","B_SUBCHANNEL_HOMEIMP", 
              "B_SUBCHANNEL_HOMEGOOD","B_SUBCHANNEL_OFFICE","B_SUBCHANNEL_PET","B_SUBCHANNEL_SPORT", 
              "B_SUBCHANNEL_BANQUETHALL","B_SUBCHANNEL_CATERING","B_SUBCHANNEL_FOODSRV", 
              "B_SUBCHANNEL_CIGCONV","B_SUBCHANNEL_CNVCONV","B_SUBCHANNEL_CNVGAS", 
              "B_SUBCHANNEL_DRIGCONV","B_SUBCHANNEL_DRUGRX","B_SUBCHANNEL_EM_ATHLETIC", 
              "B_SUBCHANNEL_EM_AUTOREPAIR","B_SUBCHANNEL_EM_BEAUTY","B_SUBCHANNEL_EM_BOOK", 
              "B_SUBCHANNEL_EM_CAMPGROUND","B_SUBCHANNEL_EM_CARWASH","B_SUBCHANNEL_EM_CORRECTIONAL", 
              "B_SUBCHANNEL_EM_DELI","B_SUBCHANNEL_EM_ETHNICMKT","B_SUBCHANNEL_EM_GIFT", 
              "B_SUBCHANNEL_EM_INTERNET","B_SUBCHANNEL_EM_MARINA","B_SUBCHANNEL_EM_SPCLTYMKT", 
              "B_SUBCHANNEL_EM_NEWS","B_SUBCHANNEL_EM_OFFICE","B_SUBCHANNEL_EM_OTHER_OFF", 
              "B_SUBCHANNEL_EM_TOBACCO","B_SUBCHANNEL_EM_RESIDENTIAL","B_SUBCHANNEL_EM_GROCERY", 
              "B_SUBCHANNEL_EM_TRADING","B_SUBCHANNEL_EM_VENDING","B_SUBCHANNEL_GRCSUPERCENTER", 
              "B_SUBCHANNEL_GRCSUPERETTE","B_SUBCHANNEL_GRCCONV","B_SUBCHANNEL_GRCLIMTED", 
              "B_SUBCHANNEL_GRCNATURAL","B_SUBCHANNEL_GRCWAREHOUSE","B_SUBCHANNEL_LIQBEER", 
              "B_SUBCHANNEL_LIQCONV","B_SUBCHANNEL_LIQSUPER","B_SUBCHANNEL_LIQWINE", 
              "B_SUBCHANNEL_MASSCONV","B_SUBCHANNEL_MASSDOLLAR","B_SUBCHANNEL_MASSGENERAL", 
              "B_SUBCHANNEL_MILCLOTHING","B_SUBCHANNEL_MILCOMM","B_SUBCHANNEL_MILCVN", 
              "B_SUBCHANNEL_MILEXCH","B_SUBCHANNEL_MILLIQUOR","B_SUBCHANNEL_MILBAR", 
              "B_SUBCHANNEL_MILDINING","B_SUBCHANNEL_MILFOODSRV","B_SUBCHANNEL_MILLODGE", 
              "B_SUBCHANNEL_MILOTHER_ON","B_SUBCHANNEL_MILREC","B_SUBCHANNEL_TRANSAIRLINE", 
              "B_SUBCHANNEL_TRANSOTHER","B_SUBCHANNEL_TRANSRAIL","B_SUBCHANNEL_UNK", 
              "B_SUBCHANNEL_CLUBCONV","B_SUBCHANNEL_CLIENT","B_SUBCHANNEL_DISTRIBUTOR", 
              "B_SUBCHANNEL_WAREHOUSE","B_SUBCHANNEL_HOUSE","B_NIELSEN_SIZE_URBAN", 
              "B_NIELSEN_SIZE_85K","B_NIELSEN_SIZE_20K_85K","B_NIELSEN_SIZE_OTHER", 
              "B_NIELSEN_SIZE_UNK","B_INDUSTRY_TOP","B_INDUSTRY_MID","B_INDUSTRY_LOW", 
              "B_INDUSTRY_UNK","B_CHAIN_Y","B_CHAIN_N","B_FRANCHISE_Y","B_FRANCHISE_N", 
              #"B_GAS_Y",#"B_GAS_N",#"B_PHARMACY_Y",#"B_PHARMACY_N",
              #"B_HI_VOL_CIG_Y",#"B_HI_VOL_CIG_N",
              "B_WINE_Y","B_WINE_N","B_LIQUOR_Y","WHITE_POP_PCT","BLACK_POP_PCT","ASIAN_POP_PCT", 
              "HISP_POP_PCT","LEAST_ACC_HISP_POP_PCT","BI_CULTURAL_HISP_POP_PCT", 
              "MOST_ACC_HISP_POP_PCT","OTHER_POP_PCT","INCOME_POP_L10_PCT","INCOME_POP_1020_PCT", 
              "INCOME_POP_2030_PCT","INCOME_POP_3040_PCT","INCOME_POP_4050_PCT", 
              "INCOME_POP_5075_PCT","INCOME_POP_75100_PCT","INCOME_POP_G100_PCT", 
              "MALE_POP_2124_PCT","MALE_POP_2534_PCT","MALE_POP_3544_PCT","MALE_POP_4554_PCT", 
              "MALE_POP_5564_PCT","MALE_POP_6574_PCT","MALE_POP_75_PCT","FEMALE_POP_2124_PCT", 
              "FEMALE_POP_2534_PCT","FEMALE_POP_3544_PCT","FEMALE_POP_4554_PCT", 
              "FEMALE_POP_5564_PCT","FEMALE_POP_6574_PCT","FEMALE_POP_75_PCT", 
              "Min Distance On-Premise" ,"Store Count On-Premise","Min Distance Off-Premise",        
              "Store Count Off-Premise" 
      ) 
    } else {
      keep = c("TDLINX_STORE_CD","B_SUBCHANNEL_BARADULT","B_SUBCHANNEL_CLUBCASUAL", "B_SUBCHANNEL_BARCOUNTRY",
               "B_SUBCHANNEL_BARIRISH", "B_SUBCHANNEL_BARLOCAL", "B_SUBCHANNEL_BARPREMIUM","B_SUBCHANNEL_CLUBPREMIUM",
               "B_SUBCHANNEL_BARSPORTS","B_SUBCHANNEL_AUTO", "B_SUBCHANNEL_BOOKS","B_SUBCHANNEL_ELEC",
               "B_SUBCHANNEL_CRAFT", "B_SUBCHANNEL_DEPT", "B_CSUBHANNEL_DISCDEPT","B_SUBCHANNEL_FAMILY",
               "B_SUBCHANNEL_HOMEIMP","B_SUBCHANNEL_OFFICE","B_SUBCHANNEL_PET","B_SUBCHANNEL_SPORT",
               "B_SUBCHANNEL_BANQUETHALL","B_SUBCHANNEL_CATERING","B_SUBCHANNEL_FOODSRV","B_SUBCHANNEL_RESTCASUAL",
               "B_SUBCHANNEL_RESTFINE","B_SUBCHANNEL_RESTTHEME","B_SUBCHANNEL_EM_ATHLETIC","B_SUBCHANNEL_EM_AUTOREPAIR",
               "B_SUBCHANNEL_EM_BEAUTY","B_SUBCHANNEL_EM_BOOK","B_SUBCHANNEL_EM_CARWASH","B_SUBCHANNEL_EM_CORRECTIONAL",
               "B_SUBCHANNEL_EM_INTERNET","B_SUBCHANNEL_EM_NEWS","B_SUBCHANNEL_EM_TOBACCO","B_SUBCHANNEL_EM_RESIDENTIAL",
               "B_SUBCHANNEL_EM_TRADING","B_SUBCHANNEL_EM_VENDING","B_SUBCHANNEL_EM_AMUSEMENT","B_SUBCHANNEL_EM_COFFEE",
               "B_SUBCHANNEL_EM_CONCESSIONAIRE","B_SUBCHANNEL_EM_GAMBLING","B_SUBCHANNEL_EM_HOTEL","B_SUBCHANNEL_EM_INSTITUTION",
               "B_SUBCHANNEL_EM_OTHERON","B_SUBCHANNEL_EM_RESTAURANT","B_SUBCHANNEL_EM_SPECIALEVNT","B_SUBCHANNEL_EM_SPORTCLUB",
               "B_SUBCHANNEL_EM_WINERY","B_SUBCHANNEL_LODGFULL","B_SUBCHANNEL_LODGLUXURY","B_SUBCHANNEL_LODGRESORT",
               "B_SUBCHANNEL_MILCLOTHING","B_SUBCHANNEL_MILCOMM","B_SUBCHANNEL_MILCVN","B_SUBCHANNEL_MILEXCH",
               "B_SUBCHANNEL_MILLIQUOR","B_SUBCHANNEL_MILBAR","B_SUBCHANNEL_MILDINING","B_SUBCHANNEL_MILFOODSRV",
               "B_SUBCHANNEL_MILLODGE","B_SUBCHANNEL_MILOTHER_ON","B_SUBCHANNEL_MILREC","B_SUBCHANNEL_RECBILLBOWL",
               "B_SUBCHANNEL_RECGAMBLING","B_SUBCHANNEL_RECCRUISE","B_SUBCHANNEL_RECGOLF","B_SUBCHANNEL_RECOTHER",
               "B_SUBCHANNEL_RECPRIVATE", "B_SUBCHANNEL_RECSTADIUM","B_SUBCHANNEL_RECTHEATER","B_SUBCHANNEL_RECTHEMEPARK",
               "B_SUBCHANNEL_TRANSAIRLINE","B_SUBCHANNEL_TRANSOTHER","B_SUBCHANNEL_TRANSRAIL","B_SUBCHANNEL_UNK",
               "B_SUBCHANNEL_CLIENT","B_SUBCHANNEL_DISTRIBUTOR","B_SUBCHANNEL_WAREHOUSE","B_SUBCHANNEL_HOUSE",
               "B_FOOD_TYPE_AMERICAN","B_FOOD_TYPE_ASIAN","B_FOOD_TYPE_CARIBBEAN","B_FOOD_TYPE_DELI",
               "B_FOOD_TYPE_DINER","B_FOOD_TYPE_EUROPEAN","B_FOOD_TYPE_FRENCH","B_FOOD_TYPE_HEALTH",
               "B_FOOD_TYPE_ITALIAN","B_FOOD_TYPE_LATIN","B_FOOD_TYPE_MIDDLEEAST","B_FOOD_TYPE_MEXICAN",
               "B_FOOD_TYPE_PIZZA","B_FOOD_TYPE_BBQ","B_FOOD_TYPE_SEAFOOD","B_FOOD_TYPE_STEAK","B_FOOD_TYPE_VARIED",
               "B_FOOD_TYPE_SOUTHWEST","B_FOOD_TYPE_CHINESE","B_FOOD_TYPE_JAPANESE","B_FOOD_TYPE_KOREAN",
               "B_FOOD_TYPE_INDIAN","B_FOOD_TYPE_VIETNAMESE","B_FOOD_TYPE_THAI","B_FOOD_TYPE_CAJUN",
               "B_FOOD_TYPE_BAGEL","B_FOOD_TYPE_BAKERY","B_FOOD_TYPE_CHICKEN","B_FOOD_TYPE_COOKIE","B_FOOD_TYPE_DONUT",
               "B_FOOD_TYPE_HAMBURGER","B_FOOD_TYPE_HOTDOG","B_FOOD_TYPE_ROASTBEEF","B_FOOD_TYPE_SANDWICH",
               "B_FOOD_TYPE_SOUPSALAD","B_FOOD_TYPE_UNK","B_NIELSEN_SIZE_URBAN","B_NIELSEN_SIZE_85K",
               "B_NIELSEN_SIZE_20K_85K","B_NIELSEN_SIZE_OTHER","B_NIELSEN_SIZE_UNK","B_INDUSTRY_TOP",
               "B_INDUSTRY_MID","B_INDUSTRY_LOW","B_INDUSTRY_UNK","B_CHAIN_Y","B_CHAIN_N","B_FRANCHISE_Y","B_FRANCHISE_N",
               "B_WINE_Y","B_WINE_N","B_LIQUOR_Y","B_LIQUOR_N","WHITE_POP_PCT","BLACK_POP_PCT",
               "ASIAN_POP_PCT","HISP_POP_PCT","LEAST_ACC_HISP_POP_PCT","BI_CULTURAL_HISP_POP_PCT",
               "MOST_ACC_HISP_POP_PCT","OTHER_POP_PCT","INCOME_POP_L10_PCT","INCOME_POP_1020_PCT",
               "INCOME_POP_2030_PCT","INCOME_POP_3040_PCT","INCOME_POP_4050_PCT","INCOME_POP_5075_PCT",
               "INCOME_POP_75100_PCT","INCOME_POP_G100_PCT","MALE_POP_2124_PCT","MALE_POP_2534_PCT",
               "MALE_POP_3544_PCT","MALE_POP_4554_PCT","MALE_POP_5564_PCT","MALE_POP_6574_PCT",
               "MALE_POP_75_PCT","FEMALE_POP_2124_PCT","FEMALE_POP_2534_PCT","FEMALE_POP_3544_PCT",
               "FEMALE_POP_4554_PCT","FEMALE_POP_5564_PCT","FEMALE_POP_6574_PCT","FEMALE_POP_75_PCT",
               "Min Distance On-Premise" ,"Store Count On-Premise" , "Min Distance Off-Premise",
               "Store Count Off-Premise"
      )
      
      
    }
    
    keep = keep[keep %in% names(final_dat_two)]
    
    # Keeping only required columns for further processing
    
    final_dat_two<-final_dat_two[,c(keep,prod_index)]
    
    
    ## Normalizing Store & Distance Columns
    
    
    final_dat_two[,c("Min Distance On-Premise" ,       
                     "Store Count On-Premise",        
                     "Min Distance Off-Premise",       
                     "Store Count Off-Premise" )][] = lapply(final_dat_two[,c("Min Distance On-Premise" ,       
                                                                              "Store Count On-Premise",        
                                                                              "Min Distance Off-Premise",       
                                                                              "Store Count Off-Premise" )],as.numeric)
    
    
    final_dat_two[,c("Min Distance On-Premise" ,       
                     "Store Count On-Premise",        
                     "Min Distance Off-Premise",       
                     "Store Count Off-Premise" )][]<-lapply(final_dat_two[,c("Min Distance On-Premise" ,       
                                                                             "Store Count On-Premise",        
                                                                             "Min Distance Off-Premise",       
                                                                             "Store Count Off-Premise")],
                                                            normalise)
    
    
    
    
    final_dat_two<-final_dat_two[c(keep,prod_index)]
    log.ind = sapply(final_dat_two, is.logical)
    final_dat_two[,log.ind] = sapply(final_dat_two[,log.ind],as.numeric)
    rm(log.ind , base)
    gc()
    
    
    
    ## Combining subchannels as per new mapping
    
    Sub_mapping <-read.csv('Subchannel_Mapping.csv',header = TRUE,stringsAsFactors = FALSE)
    Subchannel_grp<- match(names(final_dat_two), Sub_mapping$SUBCHANNEL)
    names(final_dat_two)[!is.na(Subchannel_grp)] <- as.character(Sub_mapping$GROUPS[na.omit(Subchannel_grp)])
    
    
    ## Combining foodtype as per new mapping
    
    Food_mapping <-read.csv('Foodtype_Mapping.csv',header = TRUE,stringsAsFactors = FALSE)
    Foodtype_grp<- match(names(final_dat_two), Food_mapping$FOODTYPE)
    names(final_dat_two)[!is.na(Foodtype_grp)] <- as.character(Food_mapping$GROUPS[na.omit(Foodtype_grp)])
    
    rm(Food_mapping , Sub_mapping)
    
    
    ## Combining columns with same name
    
    final_data_numeric <- final_dat_two[,sapply(final_dat_two, is.numeric)]
    names(final_data_numeric) <- gsub("\\..*", "", names(final_data_numeric))
    final_data_numeric <- as.data.frame(do.call(cbind, by(t(final_data_numeric),INDICES=names(final_data_numeric),FUN=colSums)))
    
    
    final_dat_two <- cbind(final_dat_two$TDLINX_STORE_CD,final_data_numeric)
    colnames(final_dat_two)[1] <- 'TDLINX_STORE_CD'
    rownames(final_dat_two) <- 1:nrow(final_dat_two)
    
    print("Mapping Complete")
    
    ###############################Colaborative Filtering################################
    
    newd<-final_dat_two
    colnames(newd)[-1] = gsub('X',"",colnames(newd)[-1])
    
    # Creating similarity matrix
    
    final_dat_two[is.na(final_dat_two)] <- 0
    sim1 <- cosine(as.matrix(t(final_dat_two[,-1])))
    rm(final_dat_two,final_data_numeric)
    gc()
    newd[,1] = as.character(newd[,1])
    
    
    ## Extract Table for required Premise
    
    prod_info1 <- subset(prod_info,PREMISE_TYPE_CD == prem[z])
    prod_info1 = prod_info1[prod_info1$MASTER_PKG_SKU_CD %in% prods,]
    
    vec <- c('PREMISE_TYPE_CD','TDLINX_STORE_CD','MASTER_PKG_SKU_CD','L365_TY_QTY')
    prod_info2 = prod_info1[,vec]
    prod_info2 = unique(prod_info2)
    
    prod_info2 = prod_info2[!prod_info2$TDLINX_STORE_CD %in% c(NA,"") ,]
    prod_info2 = prod_info2[prod_info2$MASTER_PKG_SKU_CD %in% gsub("X","",prod_list),]
    
    new_prod <- dcast(prod_info2,PREMISE_TYPE_CD + TDLINX_STORE_CD ~ MASTER_PKG_SKU_CD , value.var = 'L365_TY_QTY')
    
    colnames(new_prod)[grep('^[0-9]',colnames(new_prod))] = paste0(colnames(new_prod)[grep('^[0-9]',colnames(new_prod))],'_L365')
    new_prod$TDLINK_STORE_CD[] = sapply(new_prod$TDLINK_STORE_CD,zero_format)
    dat1 = plyr::join(newd , new_prod,type = 'left', match = 'all', by = 'TDLINX_STORE_CD')
    L365_names = colnames(dat1)[which(grepl("[[:digit:]]" , colnames(dat1)) & grepl('_L365', colnames(dat1)))]
    dat1 = dat1[,c('TDLINX_STORE_CD' , L365_names)]
    prod_remove = names(which(colSums(dat1[,colnames(dat1)[grep('^[0-9]',colnames(dat1))]], na.rm = T) == 0))
    dat1 = dat1[,!colnames(dat1) %in% prod_remove]
    prod_list = names(newd)[grep('^[0-9]' , names(newd)[-1])+1]
    rmprods = which(colSums(newd[,prod_list],na.rm  = T) == 0)
    not_req = unique(names(rmprods)[!names(rmprods) %in% gsub('_L365',"",names(dat1)[-1])])
    
    
    rm(prod_info2,new_prod);gc()  
    
    # Creating over product matrix for Collaborative Filtering Loop
    
    train_data = as.matrix(data.matrix(newd[,-1]))
    rownames(train_data) = 1:nrow(train_data)
    prod_list = prod_list[!prod_list %in% not_req]
    prod_mat1 = train_data[,prod_list]
    prod_mat2 = data.matrix(as(prod_mat1,'matrix'))
    rm(train_data , prod_mat1)
    rownames(sim1) = 1:dim(sim1)[1]
    colnames(sim1) = 1:dim(sim1)[2]
    gc()
    
    # Assigning values for threshold and nearest neighbors
    
    higher = 0.999999  
    lower = 0.8  
    nn=20  
    k=20 
    
    
    ## Parallelized Colaborative Filtering Loop
    
    print('Initiating Colab Loop')  
    
 
    ## Creating Iteration variables & initiating parallel backend
   
    col = ncol(prod_mat2)
    row1 = nrow(prod_mat2)    
    mcoptions = list(preschedule = TRUE)
    gc()
    cl = makeCluster(detectCores()-8)
    registerDoSNOW(cl)
    
  

    res <-  foreach(v=icount(col), 
                    .combine = combine_custom_j,
                    .packages = c('foreach','coop','iterators','data.table'),
                    .options.multicore = mcoptions) %dopar% 
                    {
                      foreach(m = icount(row1), 
                              sim = iter(sim1,by = 'row'),
                              .combine = combine_custom_i,
                              .packages = c('foreach','coop','iterators','data.table')
                      ) %do% 
                        
                      { 
                        rec_sys(m,v)
                      }                                                
                    }
    
    stopCluster(cl)
    
    
    final_mat <- data.frame(t(res[[1]]))
    store_sim = data.table(res[[2]])  
    rm(res)
    gc()
    
    print("Colaborative Filtering Complete")
    
    
    colnames(final_mat) = prod_list
    
    names(store_sim) <- c('TDLINX_STORE_CD' , 
                          'MASTER_PKG_SKU_CD' , 
                          'TDLINX_STORE_SIM_CD' , 
                          'STORE_SIM_PROD_COSINE')
    
    gc()
    
    
    #### Replication Code for New Products
    
    print('New Products')
    
    ## Initiating parallel backend for New Product Replication
    
    cl = makeCluster(detectCores() - 1)
    registerDoParallel(cl)
    
    
    prod_store_sim = foreach(sim = iter(sim1 , by = 'row') ,
                             .combine = function(...) rbindlist(list(...)),
                             .multicombine=TRUE,
                             .packages = c('coop','foreach','iterators'),
                             .options.multicore = mcoptions) %dopar%
      
                             { 
                               data.frame(
                                 names(sort(sim[,-1],decreasing = T)[1:100]),
                                 sort(sim[,-1],decreasing = T)[1:100])
                               
                             }
    
    
    
    stopCluster(cl)
    rm(sim1)
    gc()
    
    
    store_len = rep(1:nrow(prod_mat2),each = 100)[1:nrow(prod_store_sim)]     
    prod_store_sim = cbind.data.frame(store_len , prod_store_sim)
    newprod = rbindlist(replicate(length(not_req) , coredata(prod_store_sim) , simplify = F))
    req_cols = sort(replicate(nrow(newd),not_req))
    newprod$MASTER_PKG_SKU_CD = req_cols
    names(newprod)[c(1:3)] <- c("TDLINX_STORE_CD",'TDLINX_STORE_SIM_CD','STORE_SIM_PROD_COSINE')
    
    print('Creating new matrix')
    
    newprod_final_mat <- data.frame(matrix(0 , 
                                           nrow = nrow(prod_mat2) , 
                                           ncol = length(not_req) , 
                                           dimnames = list(NULL,not_req)))
    
    colnames(newprod_final_mat) <- gsub('X',"",colnames(newprod_final_mat))
    gc()
    
    
    print("New Product Complete")
    
    
    ########################## Store and Product Mapping##########################################
    
    store_map <- cbind.data.frame(1:nrow(newd),newd[,1])
    prod_map <- cbind.data.frame(1:length(colnames(prod_mat2)),colnames(prod_mat2))
    store_map[,2] = as.character(store_map[,2])
    prod_map[,2] = as.character(prod_map[,2])
    
    newprod$TDLINX_STORE_CD = store_map[,2][match(newprod$TDLINX_STORE_CD , store_map[,1])]
    newprod$TDLINX_STORE_SIM_CD = store_map[,2][match(newprod$TDLINX_STORE_SIM_CD , store_map[,1])]
    
    store_sim$TDLINX_STORE_CD = store_map[,2][match(store_sim$TDLINX_STORE_CD,store_map[,1])]
    store_sim$MASTER_PKG_SKU_CD = prod_map[,2][match(store_sim$MASTER_PKG_SKU_CD , prod_map[,1])]
    store_sim$TDLINX_STORE_SIM_CD <- store_map[,2][match(store_sim$TDLINX_STORE_SIM_CD , store_map[,1])] 
    
    
    store_sim = rbindlist(list(store_sim , newprod) , use.names = T)
    rm(newprod , prod_store_sim)
    gc()
    
    store_count = store_sim[,.(Num_Store = .N) , by = .(TDLINX_STORE_CD , MASTER_PKG_SKU_CD)]
    
    gc()
    
    ###################### Numbr of Products column ######################################
    
    num_prod = data.table(cbind.data.frame(newd[,1],rowSums(!is.na(newd[,prod_list]))))
    colnames(num_prod) = c("TDLINX_STORE_CD",'NUM_PROD')
    store_sim = num_prod[store_sim , on = 'TDLINX_STORE_CD']
    gc()
    rm(num_prod)
    
    ######################################################################################
    
    store_unique_sim = unique(store_sim[,c('TDLINX_STORE_CD','TDLINX_STORE_SIM_CD'),with = F])
    
    filtered_recom1 <- as.data.table(filtered_recom[,c('TDLINX_STORE_CD','MKT_CD',"PREMISE_TYPE_CD",'LATITUDE','LONGITUDE')])
    store_unique_sim = filtered_recom1[store_unique_sim , on = 'TDLINX_STORE_CD']
    rm(filtered_recom1)
    gc()
    
    ######################################################################################
    
    filtered_recom1 <- as.data.table(filtered_recom[,c('TDLINX_STORE_CD',"INDUSTRY_VOL_CD",'LATITUDE','LONGITUDE')])
    
    store_unique_sim = filtered_recom1[store_unique_sim , on = c('TDLINX_STORE_CD' = 'TDLINX_STORE_SIM_CD')]
    colnames(store_unique_sim)[1] <- 'TDLINX_STORE_SIM_CD'
    colnames(store_unique_sim)[colnames(store_unique_sim) == 'i.TDLINX_STORE_CD'] <- 'TDLINX_STORE_CD'
    gc()
    
    ######################################################################################
    
    
    names(store_unique_sim)[names(store_unique_sim) %in% 
                              c('LATITUDE','LONGITUDE','i.LATITUDE','i.LONGITUDE')] <- c('Sim_Store_Lat',
                                                                                         'Sim_Store_Lon',
                                                                                         'Store_Lat',
                                                                                         'Store_Lon')
    
    
    ## Finding Distance with Similar Stores
    
    
    store_unique_sim$Dist= mcmapply(store_dist , 
                                    store_unique_sim$Store_Lon,
                                    store_unique_sim$Store_Lat,
                                    store_unique_sim$Sim_Store_Lon,
                                    store_unique_sim$Sim_Store_Lat,mc.cores = detectCores()-1)
    
    
    # To kill zombie processes if any created by using mcmapply
    
    includes <- '#include <sys/wait.h>'
    code <- 'int wstat; while (waitpid(-1, &wstat, WNOHANG) > 0) {};'
    wait <- cfunction(body=code, includes=includes, convention='.C')
    wait()
    
    
    store_unique_sim[,':='(Store_Lat = NULL , 
                           Store_Lon = NULL , 
                           Sim_Store_Lat = NULL , 
                           Sim_Store_Lon = NULL)]
    gc()
    
    #Find Minimum distance for each store
    
    min_store = store_unique_sim[,.(Min_dist = min(Dist , na.rm = T)) , by = .(TDLINX_STORE_CD)]
    store_unique_sim = min_store[store_unique_sim , on = 'TDLINX_STORE_CD']
    gc()
    rm(min_store)
    
    
    store_unique_sim[,CLOSEST_STORE_SIM_IND := ifelse(Dist == Min_dist , 1,0)]
    store_unique_sim[,':='(Min_dist = NULL , Dist = NULL)]
    store_unique_sim[,LAST_UPD_DT := Sys.time()]
    
    fin = data.table(unique(filtered_recom[,c('MODEL_NM','MODEL_NO','SCEN_NM','MKT_CD')]))
    store_unique_sim = fin[store_unique_sim , on = 'MKT_CD']
    
    store_sim = store_unique_sim[store_sim , on = c('TDLINX_STORE_CD' , 'TDLINX_STORE_SIM_CD')]
    rm(store_unique_sim)
    gc()
    
    saveRDS(store_sim , paste0(state[q],'_',prem[z], '_SIM STORE.rds'))
    rm(store_sim);gc()
    
    print("Store Sim Complete")
    
    
    ########################### Predicted Ratings Unpivot #################################################
    
    
    final_mat = cbind(final_mat , newprod_final_mat)
    final_mat$TDLINX_STORE_CD = 1:nrow(final_mat)
    final_mat$TDLINX_STORE_CD <- store_map[,2][match(final_mat[,ncol(final_mat)] ,store_map[,1])]
    final_mat = melt(final_mat , id.vars = 'TDLINX_STORE_CD')
    colnames(final_mat)[-1] <- c('MASTER_PKG_SKU_CD','PRED_TPS')
    final_mat = data.table(final_mat)
    final_mat[, NEW_PROD_FLG := ifelse(MASTER_PKG_SKU_CD %in% not_req , 1, 0)]
    final_mat[, MASTER_PKG_SKU_CD := as.character(MASTER_PKG_SKU_CD)]
    
    
    ##############################################################################################
  
    
    prod1 = data.table(prod_info1[,c('TDLINX_STORE_CD','MASTER_PKG_SKU_CD','TPS')])
    prod1$MASTER_PKG_SKU_CD <- as.character(prod1$MASTER_PKG_SKU_CD)
    
    final_mat = prod1[final_mat , on = c('TDLINX_STORE_CD',"MASTER_PKG_SKU_CD")]
    final_mat = unique(final_mat) 
    
    
    filtered_recom <- data.table(filtered_recom[,c('TDLINX_STORE_CD','MKT_CD',"PREMISE_TYPE_CD","FOOD_TYPE_GROUP_DSC","CHANNEL_GROUP_DSC")])
    filtered_recom$TDLINX_STORE_CD[] = sapply(filtered_recom$TDLINX_STORE_CD,zero_format) 
    final_mat = filtered_recom[final_mat , on = 'TDLINX_STORE_CD']
    rm(filtered_recom)
    gc()
    
   
   ## Merging POD data
    
    req_pod = data.table(prod_info1[,c('TDLINX_STORE_CD','MASTER_PKG_SKU_CD','POD')])
    final_mat = req_pod[final_mat , on = c('TDLINX_STORE_CD','MASTER_PKG_SKU_CD')]
    rm(req_pod);gc()
    
    
    
    ## Adding L90_TY_QTY column
    
    prod_L90 <- data.table(prod_info1[,c('TDLINX_STORE_CD','MASTER_PKG_SKU_CD','L90_TY_QTY')])
    final_mat = prod_L90[final_mat , on = c('TDLINX_STORE_CD','MASTER_PKG_SKU_CD')]
    rm(prod_L90)
    
    
    
    ############################################################################
    
    final_mat$POD[is.na(final_mat$POD) | is.nan(final_mat$POD)] = 0
    final_mat$PRED_TPS[is.na(final_mat$PRED_TPS) | is.nan(final_mat$PRED_TPS)] = 0
    final_mat$TPS[is.na(final_mat$TPS) | is.nan(final_mat$TPS)] = 0
    
    
    ## New Additions - POD
    
    if(!state[q] %in% c('CA_S', 'CA_N')){
      
      final_mat[,PRED_TPS :=  (ifelse(NEW_PROD_FLG <  1 , PRED_TPS, PRED_TPS + POD))]
      final_mat[,TPS := ifelse(NEW_PROD_FLG < 1 , TPS  , TPS+ POD)]
      
    }
    
    
    
    ## Bayesian Ranking
    
    ## New Additions - average rating
    
    avg_rating = final_mat[NEW_PROD_FLG > 0,.(AVG_PRED_TPS = mean(PRED_TPS,na.rm = T)),by = .(TDLINX_STORE_CD)]    
    final_mat = avg_rating[final_mat , on = 'TDLINX_STORE_CD']
    final_mat = store_count[final_mat , on = c('TDLINX_STORE_CD','MASTER_PKG_SKU_CD')]
    final_mat$Num_Store[is.na(final_mat$Num_Store)] <- 0
    
    
    ## Bayesian New Formula
    
    final_mat$RECMND_RANK = (final_mat$Num_Store / (final_mat$Num_Store + 10)) * final_mat$PRED_TPS + (5/(final_mat$Num_Store + 10)) *final_mat$AVG_PRED_TPS
    final_mat[,RECMND_RANK := normalise(RECMND_RANK) , by = TDLINX_STORE_CD]
    final_mat$RECMND_RANK[is.na(final_mat$RECMND_RANK) | is.nan(final_mat$RECMND_RANK)] <- 0
    
    #rm(prod_info)
    
    final_mat$IMPACT_CD = ifelse(final_mat$RECMND_RANK > 0.7 & final_mat$RECMND_RANK <= 1,"H",
                                 ifelse(final_mat$RECMND_RANK > 0.4 & final_mat$RECMND_RANK <= 0.7, "M","L"))
    
    final_mat = fin[final_mat , on = 'MKT_CD']
    final_mat[,LAST_UPD_DT := Sys.time()]
    #final_mat[,NEW_PROD_FLG := Sys.time()]
    
    final_mat[,':='(Num_Store = NULL , POD = NULL , L90_TY_QTY = NULL , AVG_PRED_TPS = NULL)]
    
    
    saveRDS(final_mat , paste0(state[q],'_',prem[z],'_TPS VIEW.rds'))
    rm(final_mat)
    rm(new_prod , store_count , train_data , avg_rating , not_req , prod_info1)
    gc()
    
    print("TPS View Complete")
    tima2 = Sys.time()
    print(paste0('Run time for ',state[q]," ",prem[z]," - ",tima2-tima1))
    gc()
  }
}
