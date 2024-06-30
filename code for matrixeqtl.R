#Practice Matrix eQTL

install.packages("MatrixEQTL")

library("MatrixEQTL")
base.dir = find.package("MatrixEQTL");

#Set the parameters such as selected linear model
useModel = modelLINEAR; #modelANOVA or modelLINEAR or modelLINEA_CROSS

#Set names of genotype and expression data files
SNP_file_name = paste(base.dir, "/data/SNP.txt", sep="");
expression_file_name = paste(base.dir, "/data/GE.txt", sep="");

#A seperate file can be provided for covariates
covariates_file_name = paste(base.dir, "/data/Covariates.txt", sep="");
output_file_name = tempfile();

#Set the p-val threshold. For larger datasets the threshold should be lower otherwise the output will be excessive.
pvOutputThreshold = 1e-2;

#define the covariance matric for the error term. If the covariance matric is a multiple of identity, set it to numeric().
errorCovariance = numeric();

#Load genotype data
snps = SlicedData$new();
snps$fileDelimiter = "\t"; #the TAB character
snps$fileOmitCharacters = "NA"; #denote missing values
snps$fileSkipRows = 1 #one row of column labels
snps$fileSkipColumns = 1 #one column of row labels
snps$fileSliceSize = 2000 #read file in pieces of 2000 rows
snps$LoadFile(SNP_file_name)

#Load gene expression data
gene = SlicedData$new();
gene$fileDelimiter = "\t";
gene$fileOmitCharacters = "NA";
gene$fileSkipRows = 1;
gene$fileSkipColumns = 1;
gene$fileSliceSize = 2000;
gene$LoadFile(expression_file_name);

#Load covariates
cvrt = SlicedData$new();
cvrt$fileDelimiter = "\t";
cvrt$fileOmitCharacter = "NA";
cvrt$fileSkipRows = 1;
cvrt$fileSkipColumns = 1
if(length(covariates_file_name)>0) {
  cvrt$LoadFile(covariates_file_name);
  }


me = Matrix_eQTL_engine(
  snps = snps,
  gene = gene,
  cvrt = cvrt,
  output_file_name = output_file_name,
  pvOutputThreshold = pvOutputThreshold,
  useModel = useModel,
  errorCovariance = errorCovariance,
  verbose = TRUE,
  pvalue.hist = TRUE,
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = FALSE);

unlink(output_file_name);

#Results
cat('Analysis done in: ', me$time.in.sec, 'seconds', '\n');
cat('Detected eQTLs:', '\n');
show(me$all$eqtls)
plot(me)

#Each significant gene-SNP association is recorded in a separate line in the output file and in "me".
#In case of Cis/Trans eQTL analysis, two output files are produced. Every record contains a SNP name, transcript name, estimate of effect size, t- or F statistic and FDR.
#cis and trans eQTL analysis requires more information
#pvOutputThreshold.cis - p-value threshold for cis-eqtls
#output_file_name.cis - detected cis-eqtls are saved in this file
#cisDist - max distance at which gene-SNP pair is considered local
#snpspos - data frame with info about SNP locations. must hahve 3 columns: SNP name, chromosome and position
#genepos - data frame with info about gene locations. must have 4 columnd: the name, chromosome, positions of left and right ends








