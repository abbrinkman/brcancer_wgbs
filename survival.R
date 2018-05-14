#!/usr/bin/env Rscript

# parse the options
suppressPackageStartupMessages(library("optparse"))
option_list <- list(
    make_option(c("-i", "--input"), type="character", action="store", dest="infile",
        help="xlsx file with clustering membership info"))
opt <- parse_args(OptionParser(option_list=option_list))

# check if all arguments are specified, and if infile exists
if (class(opt$infile)=="NULL" ) {
  print_help(OptionParser(option_list=option_list))
  stop("incorrect arguments")
}

if( file.access(opt$infile) == -1) {
  stop(sprintf("Specified file ( %s ) does not exist", opt$infile))
}


suppressPackageStartupMessages(library(survival))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(openxlsx))

clusters <- read.xlsx(opt$infile)

pdfname <- basename(opt$infile)
pdfname <- paste0(gsub("\\.xlsx$","", pdfname), "_survival.pdf")

source("~/tools/BASIS_common_functions.txt.R")

# these columns are relevant:
# ===========================
# donor_vital_status
# donor_survival_time_in_DAYS
# donor_interval_of_last_follow.up_in_DAYS

# fix factor class into integer
clinj$donor_survival_time_in_DAYS <- as.integer(as.character(clinj$donor_survival_time_in_DAYS))
clinj$donor_interval_of_last_follow.up_in_DAYS <- as.integer(as.character(clinj$donor_interval_of_last_follow.up_in_DAYS))

# for the patients that have no value in "donor_survival_time_in_DAYS" -> 
#   take "donor_interval_of_last_follow.up_in_DAYS" i.c.w. "donor_vital_status"
s <- is.na(clinj$donor_survival_time_in_DAYS)
clinj$donor_survival_time_in_DAYS[s] <- clinj$donor_interval_of_last_follow.up_in_DAYS[s]

# take out patients for which vital status is missing
clinj <- clinj[!is.na(clinj$donor_vital_status),]
clinj <- clinj[!is.na(clinj$donor_survival_time_in_DAYS),]

# use only samples that were isolated BEFORE treatment
clinj <- clinj[clinj$sample_removed_pre_or_post.treatment=="pre",]

##create survival object. We have 'right censored' data
#time: Time to death or on-study time, weeks
#delta: Death indicator (0=alive, 1=dead)
yr = 365.25 #days
svl <- data.frame(
  'patient'=clinj$sample_name,
  'OS.yrs'=clinj$donor_survival_time_in_DAYS/yr,
  'died'=ifelse(clinj$donor_vital_status=="deceased", TRUE, FALSE),
  'ER'=clinj$ER,
  'PR'=clinj$PR,
  'HER2'=clinj$HER2,
  'grade'=clinj$tumour_grade,
  'AIMS'=clinj$AIMS,
  'cluster'=clusters$cluster[match(gsub("$","a",clinj$sample_name), clusters$patient)]
)

survivalFunc = function(DF, GROUPS, MAX_YRS) {
  colors <- c("magenta", "blue", "darkgreen","orange","darkred", "red", "green","purple","pink","brown")
  if (GROUPS=="AIMS") {
    colors <- c("Basal"="red", "Her2"="purple", "LumA"="blue", "LumB"="lightblue", "Normal"="green")
  }
  DF <- DF[DF$OS.yrs <= MAX_YRS,]
  Surv(DF$OS.yrs, DF$died) -> surv.object
  survfit(surv.object ~ DF[,GROUPS]) -> fit
  survdiff(surv.object ~ DF[,GROUPS]) -> diff
  1 - pchisq(diff$chisq, length(diff$n) - 1) -> pvalue
  plot(fit, col=colors, xlab="overall survival (years)", ylab="survival probability", 
      main=GROUPS, lwd=1, mark.time=T)
  groupnames <- sapply(strsplit(names(fit$strata), "="), function(X) {X[2]})
  legend('bottomleft', bty='n', legend=groupnames, text.col=colors, cex=1)
  legend('topright', bty='n', legend=paste('p =', signif(pvalue,3)))
}
pdf(pdfname)
par(mfrow=c(2,2))
survivalFunc(svl, "cluster", 5) #### Supplemental Figure 6C ####
survivalFunc(svl, "cluster", 10)
survivalFunc(svl, "cluster", 25)
dev.off()



