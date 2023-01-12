#===============================================================================
# Description: Develop the best climate-fire model
#===============================================================================

#===============================================================================
# 1). Preliminary -----
#===============================================================================

# Clean up
rm(list = ls())
graphics.off()
gc()

# Packages
wants <- c("R.matlab")
needs <- wants[!(wants %in% installed.packages()[, "Package"])]
if (length(needs))
  install.packages(needs)
lapply(wants, function(i)
  require(i, character.only = TRUE))
rm(needs, wants)

# Directories
dir <- list()
dir$root <- '~/Dropbox/estcena/scripts/fires_california/'
dir$data = '~/Dropbox/estcena/scripts/fires_california/data_def/'
# dir$out = paste0(dir$data, '/climate-fire-models/')

# Misc
misc <- list()
misc$years_model <- 1971:2021

#===============================================================================
# 2). load data -----
#===============================================================================

## load data
dum = readMat(paste0(dir$data, "fires/frap_forest_sierra_ncoast_year.mat"))
BAS = as.numeric(dum$FIRE) #1971:2020
dum = readMat(paste0(dir$data, "nclimgrid/nclimgrid.mat"))
tsmax = dum$nclimgrid[1, (11 * 12 + 1):dim(dum$nclimgrid)[2]] #period 1971:2021
prec = dum$nclimgrid[2, (10 * 12 + 1):dim(dum$nclimgrid)[2]] #period 1970:2021

## find the best predictor
num_iter = 5005
corr = array(NA, dim = c(num_iter))
sig = array(NA, dim = c(num_iter))
best_start_pr = array(NA, dim = c(num_iter))
best_stop_pr = array(NA, dim = c(num_iter))
best_start_tx = array(NA, dim = c(num_iter))
best_stop_tx = array(NA, dim = c(num_iter))
best_x1 = array(NA, dim = c(num_iter, length(misc$years_model)))
best_x2 = array(NA, dim = c(num_iter, length(misc$years_model)))
best_x3 = array(NA, dim = c(num_iter, length(misc$years_model)))

best_x1_1 = array(NA, dim = c(num_iter, length(misc$years_model)))
best_x1_12 = array(NA, dim = c(num_iter, length(misc$years_model)))
best_x1_13 = array(NA, dim = c(num_iter, length(misc$years_model)))
best_x1_123 = array(NA, dim = c(num_iter, length(misc$years_model)))

best_x2_12 = array(NA, dim = c(num_iter, length(misc$years_model)))
best_x2_123 = array(NA, dim = c(num_iter, length(misc$years_model)))

best_x3_13 = array(NA, dim = c(num_iter, length(misc$years_model)))
best_x3_123 = array(NA, dim = c(num_iter, length(misc$years_model)))

# years = misc$years_model
# years = scale(years_model)

# climdiv = array(data = NA,dim = c(12,dim(MyData)[1])) #spi1,2,3,4,5,6,9,12,24,tasmax.tasmin,tas
#iscC = c(3, 6, 8)


k = 0
for (im in 1:10) {
  tmp_stop = im:10
  for (istop in 1:length(tmp_stop)) {
    tx = rep(NA, length(misc$years_model))
    for (iyear in 1:length(misc$years_model)) {
      i1 = (iyear - 1) * 12 + im
      i2 = (iyear - 1) * 12 + tmp_stop[istop]
      tx[iyear] = mean(tsmax[i1:i2], na.rm = TRUE)
    }
    
    for (im_pr in 10:22) {
      tmp_stop_pr = im_pr:22
      for (istop_pr in 1:length(tmp_stop_pr)) {
        k = k + 1
        print(k)
      
        pr = rep(NA, length(misc$years_model))
        for (iyear in 1:length(misc$years_model)) {
          i1 = (iyear - 1) * 12 + im_pr
          i2 = (iyear - 1) * 12 + tmp_stop_pr[istop_pr]
          pr[iyear] = mean(prec[i1:i2], na.rm = TRUE)
        }
        
        
        pre1 = vector()
        pre12 = vector()
        pre13 = vector()
        pre123 = vector()
        
        for (iy in 1:length(misc$years_model)) {
          BAS_train = BAS[-iy]
          SPI_train = pr[-iy]
          tx_train = tx[-iy]
          years_train = misc$years_model[-iy]
          
          SPI_test = pr[iy]
          tx_test = tx[iy]
          years_test = misc$years_model[iy]
          
          mydata_train = data.frame(
            "y" = log(BAS_train),
            "x1" = (years_train),
            "x2" = (tx_train),
            "x3" = (SPI_train)
          )
          
          mydata_test = data.frame(
            "x1" = (years_test),
            "x2" = (tx_test),
            "x3" = (SPI_test)
          )
          
          fit1 <- lm(y ~ x1 , data = mydata_train)
          fit12 <- lm(y ~  x1 + x2 , data = mydata_train)
          fit13 <- lm(y ~  x1 + x3, data = mydata_train)
          fit123 <- lm(y ~  x1 + x2 + x3, data = mydata_train)
          
          pre1[iy] = fit1$coefficients[1] + fit1$coefficients[2] * mydata_test$x1 
          pre12[iy] = fit12$coefficients[1] + fit12$coefficients[2] * mydata_test$x1 + fit12$coefficients[3] * mydata_test$x2 
          pre13[iy] = fit13$coefficients[1] + fit13$coefficients[2] * mydata_test$x1 + fit13$coefficients[3] * mydata_test$x3
          pre123[iy] = fit123$coefficients[1] + fit123$coefficients[2] * mydata_test$x1 + fit123$coefficients[3] * mydata_test$x2 + fit123$coefficients[4] * mydata_test$x2 
          
          # best_x1_1[k, iy] = fit1$coefficients[2]
          # best_x1_12[k, iy] = fit12$coefficients[2]
          # best_x1_13[k, iy] = fit13$coefficients[2]
          # best_x1_123[k, iy] = fit123$coefficients[2]
          #
          # best_x2_12[k, iy] = fit12$coefficients[3]
          # best_x2_123[k, iy] = fit123$coefficients[3]
          #
          # best_x3_13[k, iy] = fit13$coefficients[3]
          # best_x3_123[k, iy] = fit123$coefficients[4]
          
          
        }
        
        rho1 = cor.test(log(BAS),
                        pre1,
                        use = "pairwise.complete.obs",
                        alternative = "greater")
        rho12 = cor.test(log(BAS),
                         pre12,
                         use = "pairwise.complete.obs",
                         alternative = "greater")
        rho13 = cor.test(log(BAS),
                         pre13,
                         use = "pairwise.complete.obs",
                         alternative = "greater")
        rho123 = cor.test(log(BAS),
                          pre123,
                          use = "pairwise.complete.obs",
                          alternative = "greater")
        if (rho1$estimate > rho12$estimate &
            rho1$estimate > rho13$estimate &
            rho1$estimate > rho123$estimate) {
          corr[k] = rho1$estimate
          sig[k] = rho1$p.value
          # best_x1[k, ] = best_x1_1[k, ]
        }
        if (rho12$estimate > rho1$estimate &
            rho12$estimate > rho13$estimate &
            rho12$estimate > rho123$estimate) {
          corr[k] = rho12$estimate
          sig[k] = rho12$p.value
          # best_x1[k, ] = best_x1_12[k, ]
          # best_x2[k, ] = best_x2_12[k, ]
          best_start_tx[k] = im
          best_stop_tx[k] = tmp_stop[istop]
        }
        if (rho13$estimate > rho1$estimate &
            rho13$estimate > rho12$estimate &
            rho13$estimate > rho123$estimate) {
          corr[k] = rho13$estimate
          sig[k] = rho13$p.value
          # best_x1[k, ] = best_x1_13[k, ]
          # best_x3[k, ] = best_x3_13[k, ]
          best_start_pr[k] = im_pr
          best_stop_pr[k] = tmp_stop_pr[istop_pr]
        }
        if (rho123$estimate > rho1$estimate &
            rho123$estimate > rho12$estimate &
            rho123$estimate > rho13$estimate) {
          corr[k] = rho123$estimate
          sig[k] = rho123$p.value
          # best_x1[k, ] = best_x1_123[k, ]
          # best_x2[k, ] = best_x2_123[k, ]
          # best_x3[k, ] = best_x3_123[k, ]
          best_start_tx[k] = im
          best_stop_tx[k] = tmp_stop[istop]
          best_start_pr[k] = im_pr
          best_stop_pr[k] = tmp_stop_pr[istop_pr]
          
        }
        
      }
    }
  }
}


sig = p.adjust(sig, method = "fdr")
corr[sig > 0.01] = NA



dum = max((corr), na.rm = TRUE)
idx = which(abs(corr) == dum, arr.ind = TRUE)




best_corr = corr[idx]
best_start_pr_fin = best_start_pr[idx]
best_stop_pr_fin = best_stop_pr[idx]
best_start_tx_fin = best_start_tx[idx]
best_stop_tx_fin = best_stop_tx[idx]
best_sig = sig[idx]

# best_x1_fin = best_x1[idx,]
# best_x2_fin = best_x2[idx,]
# best_x3_fin = best_x3[idx,]


# save(best_x1_fin, file = file.path(dir_out, paste0('best_x1_fin_', model_name, '.RData')))
# save(best_x2_fin, file = file.path(dir_out, paste0('best_x2_fin_', model_name, '.RData')))
# save(best_x3_fin, file = file.path(dir_out, paste0('best_x3_fin_', model_name, '.RData')))



rm(tx)
tx = rep(NA, length(misc$years_model))
pr = rep(NA, length(misc$years_model))
pr2=pr
for (iyear in 1:length(misc$years_model)) {
  i1 = (iyear - 1) * 12 + best_start_tx_fin
  i2 = (iyear - 1) * 12 + best_stop_tx_fin
  tx[iyear] = mean(tsmax[i1:i2], na.rm = TRUE)
  i1 = (iyear ) * 12 + best_start_tx_fin
  i2 = (iyear ) * 12 + best_stop_tx_fin
  pr2[iyear] = mean(prec[i1:i2], na.rm = TRUE)
  i1 = (iyear -1) * 12 + best_start_pr_fin
  i2 = (iyear -1) * 12 + best_stop_pr_fin
  pr[iyear] = mean(prec[i1:i2], na.rm = TRUE)
}

cor.test((tx-mean(tx)),(pr2-mean(pr2)))
plot.ts(pr2)
plot((tx-mean(tx)),(pr-mean(pr)))
cor.test((tx-mean(tx)),(pr-mean(pr)))
## all the series
mydata_reg = data.frame(
  "y" = log(BAS / 1000),
  "x1" = misc$years_model,
  "x2" = tx,
  "x3" = pr
)
fit <- lm(y ~ x1 + x2 + x3, data = mydata_reg)
# confint.default(fit)
summary(fit)
pred = fit$coefficients[1] + fit$coefficients[2] *
  mydata_reg$x1 + fit$coefficients[3] *
  mydata_reg$x2 + fit$coefficients[4] *
  mydata_reg$x3 

# dev.off()
# par(mfrow=c(1,3))
# plot(log(BAS))
# plot.ts(log(BAS))
plot.ts(log(BAS / 1000))
# lines(log(BAS))
lines((pred), col = "red")
print(cor.test(log(BAS), pred))


AIC(fit)

fit <- lm(y ~ x2, data = mydata_reg)
# confint.default(fit)
summary(fit)
AIC(fit)

fit <- lm(y ~ x2 + x3, data = mydata_reg)
# confint.default(fit)
summary(fit)

fit <- lm(y ~ x2 + x2*x3, data = mydata_reg)
# confint.default(fit)
summary(fit)
