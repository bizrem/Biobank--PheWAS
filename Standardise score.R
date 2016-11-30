#Read in GRS generated in plink that has NonEuros, recommended and relateds already excluded
GRS <- read.table("Y:/Biobank/GRS_sex_mismatch_not_removed/bb_GRS_from_final_corrected.profile", header=TRUE)

#remove unwanted columns: IID, pheno, CNT, CNT2; leaving FID and score
GRS <- GRS[,-c(2:5)]

#create mean and SD of score as vectors
mean <- mean(GRS$SCORE)
sd <- sd(GRS$SCORE)

#create standardized score using z transformation - putting in as a column in GRS df
GRS$SDscore <- ((GRS$SCORE - mean)/sd)

#histogram of SDscore - for look-see-ies
hist(GRS$SDscore) 

#rename FID to id - to match phewas file
library(reshape)
GRS <- rename(GRS, c(FID="id"))

#save - to be able to load it later on
save(GRS, 
     file = "O:/EBI fellowship/Biobank/GRS/SDscore.RData")


save(phewas, file="O:/EBI fellowship/Biobank/analysis")

###Remove sex-mismatch
sex_mismatch_IDs <- read.table("Y:/Biobank/PRS_biobank/biobank_exclusions/sex_mismatch.txt")
sexmisid <- sex_mismatch_IDs$V1

for(i in sexmisid){
  phewas[!phewas$id == i, ]
}


# read in linker file and match the phenotype and genotype ids
linker = read.csv("Y:/Biobank/ukb6531.csv", header=T, dec=",")
phewas$fid <- linker$app.8786[match(phewas$id, linker$app.10074)] #generates a new variable matching the pheno and geno ids
phewas$id<-NULL # delete the original id in phews
phewas<- rename(phewas, c(fid= "id")) #rename the new variable to id so can merge with genetic data
phewas <- merge(GRS, phewas, by= "id") #merge the GRS df in to phewas df

###Remove sex-mismatch
sex_mismatch_IDs <- read.table("Y:/Biobank/PRS_biobank/biobank_exclusions/sex_mismatch.txt")
sexmisid <- sex_mismatch_IDs$V1

for(i in sexmisid){
  phewas[!phewas$id == i, ]
}



library(readstata13)
genetic_vars <- read.dta13("Y:/Biobank/biobank_genotype_supp_NMD_150417.dta")
genetic_vars$chip <- ifelse(genetic_vars$n_22000_0_0<0, 1,
                            ifelse(genetic_vars$n_22000_0_0<2000, 2, NA))

phewas$chip <- genetic_vars$chip[match(phewas$id, genetic_vars$n_eid)]



###Need to remove sex mismatch and also withdrawn consent


save(GRS, bd, ncillness, phewas, NCIllness, nonccodes, noncnames, linker, pheno_glm_output, ncn,
     file = "O:/EBI fellowship/Biobank/GRS/SDscore.RData")




vars <- names(phewas) # generates a list of the names of the columns in the data frame to analyse
vars <- vars[-1:-3] # removes aln, qlet, SDscore from the list so glm won't be performed on those.
ncn <- lapply(vars, function(x){
  glm(phewas[,x] ~ phewas$SDscore, na.action=na.exclude, family=binomial(link= "logit"), maxit=100)
})


#creates a data frame which will be an intermediatery in the loop - creates it using the first glm output as correct columns in there - 5 columns
library(broom)
glm_out <- tidy(ncn[[1]])
#create a final output dataframe which will include all statistical findings I want from the glm 
#(only the second row from tidy -ed transformed dataframe - so will have 1 row per variable in the final table 'pheno_glm_output3)
pheno_glm_output <- data.frame(matrix(NA, nrow = length(vars), ncol = 5))
#loop through the glm_output list and transform the statistical findings into a dataframe and extract the second row and put it in 'glm_output3'
for (i in 1:length(ncn)){
  pheno_glm_output2 <- tidy(ncn[[i]])
  pheno_glm_output[i,]<- pheno_glm_output2[2,]
}

#Tidy up the data frame
colnames(pheno_glm_output)<-colnames(pheno_glm_output2) #set the column names as the statistical findings from tidy()
row.names(pheno_glm_output) <- vars #set the row names as the variables that were analysed
pheno_glm_output$term <- pheno_glm_output$statistic <- NULL #get rid of the statistics 'term'-just the SDscore names and 'statistic' - the t value


library(xlsx)
write.xlsx(pheno_glm_output, file= "O:/EBI fellowship/Biobank/ncn.xlsx")

fam$sex.c <- fam$sex

# Get list of IDs with mismatch between reported (sex.r) and genetic (sex.g) sex
fam$sex.r <- fam$sex 
gen_sex <- read.table("/panfs/panasas01/shared-biobank/data/derived/genetic_sex/data.txt")
fam$sex.g <- gen_sex$V3[match(fam$FID, gen_sex$V1)]
sex_mismatch_IDs <- fam$FID[which (fam$sex.r==1 & fam$sex.g==0 | fam$sex.r==2 & fam$sex.g==1)]
length(sex_mismatch_id)
# 182 IDs
sex_mismatch <- as.data.frame(sex_mismatch_IDs)
sex_mismatch$V2 <- sex_mismatch$sex_mismatch_IDs
write.table(sex_mismatch, file="../data/sex_mismatch.txt", row.names=F, col.names=F)


library(readstata13)
genetic_vars <- read.dta13("../data/biobank_genotype_supp_NMD_150417.dta")
genetic_vars$chip <- ifelse(genetic_vars$n_22000_0_0<0, 1,
                            ifelse(genetic_vars$n_22000_0_0<2000, 2, NA))

fam$chip <- genetic_vars$chip[match(fam$FID, genetic_vars$n_eid)]

fam <- fam[c("FID", "IID", "PID", "MID", "sex", "ECZ", "sex.c", "chip")]
write.table(fam, file="../data/ukbiobank_ecz.pheno", row.names=F, quote=F)
