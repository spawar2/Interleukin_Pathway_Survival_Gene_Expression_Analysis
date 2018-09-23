library(Biobase)
library(GEOquery)
library(survminer)
library(survival)

#Pancreatic Cancer
gset <- getGEO("GSE21501", GSEMatrix =TRUE, getGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL4133", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]
p <- pData(gset)

status <- na.omit(as.numeric(as.matrix(p[,50])))
time <- na.omit(as.numeric(as.matrix(p[,51])))

# Low-expression pancreatic cancer patients
except IL1A 23, 56
IL1B 33
IL1R1 28, 33, 36, 39, 42, 47, 64
IL1R2 Nothing
IL1RAP 5, 15, 57, 86, 87
IL1RN 1, 2, 3, 4, 5, 38, 76, 81
MYD88 32, 33, 35, 38, 39, 41, 43, 46, 51, 55, 56, 61, 62, 63

Less_exprs_Status <- status[c(23, 56, 33, 28, 36, 39, 42, 47, 64, 5, 15, 57, 86, 87, 1, 2, 3, 4, 5, 38, 76, 81, 32, 33, 35, 38, 39, 41, 43, 46, 51, 55, 56, 61, 62, 63)]
Less_exprs_Time <- time[c(23, 56, 33, 28, 36, 39, 42, 47, 64, 5, 15, 57, 86, 87, 1, 2, 3, 4, 5, 38, 76, 81, 32, 33, 35, 38, 39, 41, 43, 46, 51, 55, 56, 61, 62, 63)]

More_exprs_Status <- status[-c(23, 56, 33, 28, 36, 39, 42, 47, 64, 5, 15, 57, 86, 87, 1, 2, 3, 4, 5, 38, 76, 81, 32, 33, 35, 38, 39, 41, 43, 46, 51, 55, 56, 61, 62, 63)]
More_exprs_Time <- time[-c(23, 56, 33, 28, 36, 39, 42, 47, 64, 5, 15, 57, 86, 87, 1, 2, 3, 4, 5, 38, 76, 81, 32, 33, 35, 38, 39, 41, 43, 46, 51, 55, 56, 61, 62, 63)]

m <- exprs(gset)
data_Less <- m[,c(23, 56, 33, 28, 36, 39, 42, 47, 64, 5, 15, 57, 86, 87, 1, 2, 3, 4, 5, 38, 76, 81, 32, 33, 35, 38, 39, 41, 43, 46, 51, 55, 56, 61, 62, 63)]
data_More <- m[,-c(23, 56, 33, 28, 36, 39, 42, 47, 64, 5, 15, 57, 86, 87, 1, 2, 3, 4, 5, 38, 76, 81, 32, 33, 35, 38, 39, 41, 43, 46, 51, 55, 56, 61, 62, 63)]

data_Less <- as.data.frame(data_Less)
data_More <- as.data.frame(data_More)

surv_object <- Surv(More_exprs_Time, More_exprs_Status)
fit <- survfit(surv_object ~ More_exprs_Status, data=data_More)
plot(fit, col = "red", xlab="Survival time (days)", ylab="Survival Probability")

surv_object <- Surv(Less_exprs_Time, Less_exprs_Status)
fit2 <- survfit(surv_object ~ Less_exprs_Status, data=data_Less)
lines(fit2, col="green")

# P value is insignificant########
fit.list <- list(Normal = fit2, Tumor = fit)
surv_pvalue(fit.list, combine = TRUE)

#Lung Cancer
gset <- getGEO("GSE37745", GSEMatrix =TRUE, getGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL570", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]
p <- pData(gset)

status <- as.matrix(p[,45])
status[status == "no"] <- 0
status[status == "yes"] <- 1
status <- as.numeric(status)
time <- as.numeric(as.matrix(p[,43]))
# convert days to months
f <- function(x) x*0.032855
time <- f(time)

surv_object <- Surv(time, status)
fit2 <- survfit(surv_object ~ status)
lines(fit2, col="blue")

#Ovarian Cancer
gset <- getGEO("GSE14764", GSEMatrix =TRUE, getGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL96", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]
p <- pData(gset)

status <- as.numeric(as.matrix(p[,40]))
time <- as.numeric(as.matrix(p[,41]))

surv_object <- Surv(time, status)
fit2 <- survfit(surv_object ~ status)
lines(fit2, col="orange")

#Head and Neck Cancer
gset <- getGEO("GSE65858", GSEMatrix =TRUE, getGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL10558", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]
p <- pData(gset)

status <- as.matrix(p[,108])
status[status == "FALSE"] <- 0
status[status == "TRUE"] <- 1
status <- as.numeric(status)
time <- as.numeric(as.matrix(p[,107]))
# convert days to months
f <- function(x) x*0.032855
time <- f(time)

surv_object <- Surv(time, status)
fit2 <- survfit(surv_object ~ status)
lines(fit2, col="black")
legend("topright", legend=c("Pancreatic: High", "Pancreatic: Low", "Lung", "Ovarian", "Head and Neck"),
       col=c("red", "Green", "blue", "orange", "black"), lty=1:2, cex=0.8)

























