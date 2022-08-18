
#function to install packages needed 
check_pkg_deps <- function() {
  if(!require(formattable)) {
    message("installing the 'formattable' package")
    install.packages("formattable")
  }
  if(!require(tidyverse)) {
    message("installing the 'tidyverse' package")
    install.packages("tidyverse")
  }
  if(!require(Hmisc)) {
    message("installing the 'Hmisc' package")
    install.packages("Hmisc")
  }
  if(!require(here)) {
    message("installing the 'here' package")
    install.packages("here")
  }
  if(!require(htmltools)) {
    message("installing the 'htmltools' package")
    install.packages("htmltools")
  }
  if(!require(webshot)) {
    message("installing the 'webshot' package")
    install.packages("webshot")
    #webshot
    webshot::install_phantomjs(force = FALSE)
  }
}

#sum weights function 
sum_Weights <- function(weight) {
  if(!sum(weight)==1) {
    message("Scaling weights")
    scaled_weight <- weight / sum(weight)
  } else {
    scaled_weight <- weight 
  }
}

## effective sample size
get_sample_size <- function(treatment=NULL, weight) {
  scaled_weight <- sum_Weights(weight)
  if(is.null(treatment)) {
    samplesize <- (1/sum(scaled_weight^2))*sum(scaled_weight)^2
  } else {
    samplesize <- (1/sum(scaled_weight[treatment]^2))*sum(scaled_weight[treatment])^2
  }
}

#export baseline table
export_formattable <- function(f, file, width = NULL, height = NULL, 
                               background = "white", delay = 0.2)
{
  w <- as.htmlwidget(f, width = width, height = height)
  path <- html_print(w, background = background, viewer = NULL)
  url <- paste0("file:///", gsub("\\\\", "/", normalizePath(path)))
  webshot(url,
          file = file,
          selector = ".formattable_widget",
          delay = delay)
}

comma <- function(x) {
  format(x, digits = 2, big.mark = ",")
}

comma(4000)

# Documentation
# dat name of dataframe
# treatment: vector that contains the treatment information, where 0=not treated, 1=treated
# cont: vector of the continuous variables you wish to use: e.g. c("age", "bmi", etc...)
# cat: vector of the dichotomous variables you wish to use: e.g. c("sex","CHF", etc...)
# weight: vector containing the weights
# out_type: 1=just overall, 2=just treated and not treated (DEFAULT), 3=overall and treated/not treated, will default to treated / not_treated
# width: scale of the output png, I prefer "70%"

baseline_table <- function(dat,treatment,cont,cat,weight,out_type=2,width=NULL){
  
  #check to see if necessary packages are installed
  check_pkg_deps()
  
  # first is set dat as tibble
  dat <- as_tibble(dat)
  
  #making the dataframes of the continuous and dichotomous variables
  xcont <- dat %>% select(all_of(cont))
  xcat <- dat %>% select(all_of(cat))
  
  #label names for final table
  labelnames <- c(cont, cat)
  
  #setting directory to location of program
  here::here()
  
  
  
  #effective sample size function which calls sum weights function
  
  ss.treated <- round(get_sample_size(treatment=treatment==1,weight=weight), 0)
  ss.not_treated <- round(get_sample_size(treatment=treatment==0,weight=weight), 0)
  ss.ovrall <- ss.treated+ss.not_treated
  
  #scale the weights for the descriptive table
  scaled_w <-  treatment*weight*mean(treatment)/sum(treatment*weight) + (1-treatment)*weight*mean(1-treatment)/sum((1-treatment)*weight)
  
  #Treated
  med.cont.treated <- map_dbl(xcont[treatment==1,], wtd.quantile, probs=0.50, weights=scaled_w[treatment==1], type='i/n', normwt=FALSE)
  lowest.cont.treated <- map_dbl(xcont[treatment==1,], wtd.quantile, probs=0.05, weights=scaled_w[treatment==1], type='i/n', normwt=FALSE)
  low.cont.treated <- map_dbl(xcont[treatment==1,], wtd.quantile, probs=0.25, weights=scaled_w[treatment==1], type='i/n', normwt=FALSE)
  up.cont.treated <- map_dbl(xcont[treatment==1,], wtd.quantile, probs=0.75, weights=scaled_w[treatment==1], type='i/n', normwt=FALSE)
  highest.cont.treated <- map_dbl(xcont[treatment==1,], wtd.quantile, probs=0.95, weights=scaled_w[treatment==1], type='i/n', normwt=FALSE)
  cont.treated <-  as_tibble(round(cbind(first=med.cont.treated,lowest=lowest.cont.treated, second=low.cont.treated,third=up.cont.treated, highest=highest.cont.treated),1))
  cont.treated <- cont.treated %>% mutate(Treated = paste(first, " ", "(", lowest, ", ", second, ", ", third, ", ", highest, ")", sep="")) %>%
    select(Treated)
  
  prop.cat.treated <- map_dbl(xcat[treatment==1,], wtd.mean, weights=scaled_w[treatment==1],normwt="ignored")
  n.cat.treated <- round(prop.cat.treated*ss.treated, 0)
  
  cat.treated <- as_tibble(cbind(first=round(n.cat.treated,0),second=prop.cat.treated))
  cat.treated <- cat.treated %>% mutate(Treated = paste(comma(first)," ", "(",percent(second,1),")", sep="")) %>%
    select(Treated)
  treated <- as_tibble(rbind(cont.treated,cat.treated))
  
  #Not Treated
  med.cont.not_treated <- map_dbl(xcont[treatment==0,], wtd.quantile, probs=0.50, weights=scaled_w[treatment==0], type='i/n', normwt=FALSE)
  lowest.cont.not_treated <- map_dbl(xcont[treatment==0,], wtd.quantile, probs=0.05, weights=scaled_w[treatment==0], type='i/n', normwt=FALSE)
  low.cont.not_treated <- map_dbl(xcont[treatment==0,], wtd.quantile, probs=0.25, weights=scaled_w[treatment==0], type='i/n', normwt=FALSE)
  up.cont.not_treated <- map_dbl(xcont[treatment==0,], wtd.quantile, probs=0.75, weights=scaled_w[treatment==0], type='i/n', normwt=FALSE)
  highest.cont.not_treated <- map_dbl(xcont[treatment==0,], wtd.quantile, probs=0.95, weights=scaled_w[treatment==0], type='i/n', normwt=FALSE)
  cont.not_treated <-  as_tibble(round(cbind(first=med.cont.not_treated,lowest=lowest.cont.not_treated, second=low.cont.not_treated,third=up.cont.not_treated, highest=highest.cont.not_treated),1))
  cont.not_treated <- cont.not_treated %>% mutate(Not_treated = paste(first, " ", "(", lowest, ", ", second, ", ", third, ", ", highest, ")", sep="")) %>%
    select(Not_treated)
  
  prop.cat.not_treated <- map_dbl(xcat[treatment==0,], wtd.mean, weights=scaled_w[treatment==0],normwt="ignored")
  n.cat.not_treated <- round(prop.cat.not_treated*ss.not_treated, 0)
  
  cat.not_treated <- as_tibble(cbind(first=round(n.cat.not_treated,0),second=prop.cat.not_treated))
  cat.not_treated <- cat.not_treated %>% mutate(Not_treated = paste(comma(first)," ", "(",percent(second,1),")", sep="")) %>%
    select(Not_treated)
  not_treated <- as_tibble(rbind(cont.not_treated,cat.not_treated))
  not_treated <- not_treated %>% rename("Not Treated" = Not_treated)
  
  #overall
  med.cont.ovrall <- map_dbl(xcont, wtd.quantile, probs=0.50, weights=scaled_w, type='i/n', normwt=FALSE)
  lowest.cont.ovrall <- map_dbl(xcont, wtd.quantile, probs=0.05, weights=scaled_w, type='i/n', normwt=FALSE)
  low.cont.ovrall <- map_dbl(xcont, wtd.quantile, probs=0.25, weights=scaled_w, type='i/n', normwt=FALSE)
  up.cont.ovrall <- map_dbl(xcont, wtd.quantile, probs=0.75, weights=scaled_w, type='i/n', normwt=FALSE)
  highest.cont.ovrall <- map_dbl(xcont, wtd.quantile, probs=0.95, weights=scaled_w, type='i/n', normwt=FALSE)
  cont.ovrall <-  as_tibble(round(cbind(first=med.cont.ovrall,lowest=lowest.cont.ovrall, second=low.cont.ovrall,third=up.cont.ovrall, highest=highest.cont.ovrall),1))
  cont.ovrall <- cont.ovrall %>% mutate(Overall = paste(first, " ", "(", lowest, ", ", second, ", ", third, ", ", highest, ")", sep="")) %>%
    select(Overall)
  
  #just need to sum the treated/non treated for this and do numerator/sample size
  
  #so i need to sum n.cat.treated and n.cat.not_treated for number overall treated
  #then I need to divide this by the ss.ovrall to get the proportion
  cat.ovrall.num <- n.cat.treated + n.cat.not_treated
  cat.ovrall.prop <- cat.ovrall.num / ss.ovrall
  
  cat.ovrall <- as_tibble(cbind(first=cat.ovrall.num,second=cat.ovrall.prop))
  cat.ovrall <- cat.ovrall %>% mutate(Overall = paste(comma(first)," ", "(",percent(second,1),")", sep="")) %>%
    select(Overall)
  overall <- as_tibble(rbind(cont.ovrall,cat.ovrall))
    
  #if/then to print based on out_type
  
  if (near(out_type,1)) { #overall
    
    ess.1 <- tibble(Characteristic = "Effective Sample Size", Overall = comma(ss.ovrall))
    
    total <- overall %>% mutate(Characteristic = labelnames) %>%
      select(Characteristic, Overall)
    
    final <- rbind(ess.1, total)
    
    final_table <- formattable(final,
                               align=c("r","c"),
                               list(Characteristic = formatter(
                                 "span", style = ~ style(color = "grey",font.weight = "bold")) 
                               ))
    export_formattable(final_table,"Baseline Table.png",width=width)
    write.csv(final_table, file="Baseline_table.csv", row.names = F)
    return(final_table)
    
  } else if (near(out_type,2)) {  #treated and not treated
    
    ess.2 <- tibble(Characteristic = "Effective Sample Size", Treated = comma(ss.treated), `Not Treated` = comma(ss.not_treated))
    
    total <- as_tibble(cbind(treated, not_treated))
    total <- total %>% mutate(Characteristic = labelnames) %>%
      select(Characteristic, Treated, `Not Treated`)
    
    final <- rbind(ess.2, total)
    
    final_table <- formattable(final,
                               align=c("r","c","c"),
                               list(Characteristic = formatter(
                                 "span", style = ~ style(color = "grey",font.weight = "bold")) 
                               ))
    export_formattable(final_table,"Baseline Table.png",width=width)
    write.csv(final_table, file="Baseline_table.csv", row.names = F)
    return(final_table)
    
  } else { #overall and treated/not treated
    
    ess.3 <- tibble(Characteristic = "Effective Sample Size", Overall = comma(ss.ovrall), Treated = comma(ss.treated), `Not Treated` = comma(ss.not_treated))
    
    total <- as_tibble(cbind(overall, treated, not_treated))
    total <- total %>% mutate(Characteristic = labelnames) %>%
      select(Characteristic, Overall, Treated, `Not Treated`)
    
    final <- rbind(ess.3, total)
    
    final_table <- formattable(final,
                               align=c("r","c","c","c"),
                               list(Characteristic = formatter(
                                 "span", style = ~ style(color = "grey",font.weight = "bold")) 
                               ))
    export_formattable(final_table,"Baseline Table.png",width=width)
    write.csv(final_table, file="Baseline_table.csv", row.names = F)
    return(final_table)
  }
  
  
}
