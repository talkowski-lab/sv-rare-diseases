suppressPackageStartupMessages({
    require(edgeR,quietly=T)
    require(bsseq,quietly=T)
    require(optparse,quietly=T)
    require(BiocParallel,quietly=T)
})

# Define the command line options
option_list <- list(
    make_option(c("-d", "--dir"), type="character", default=NULL, 
              help="Directory path that houses the bismark-formatted input files", metavar="character"),
    make_option(c("-o", "--outputdir"), type="character", default=NULL, help="Output directory", metavar="character"),
    make_option(c("-p", "--pedfile"), type="character", default=NULL, help="PED file for dataset", metavar="character"),
    make_option(c("-t", "--threads"), type="integer", default=1, help="Number of threads to use", metavar="integer"))

# Parse the command line arguments
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

message("Reading PED file")
# read in the PED file and process it to get affected sample(s)
ped <- read.table(opt$pedfile, header=F, stringsAsFactors=F)
affected <- ped[ped$V6==2,2] # get individual sample names of affected samples

message("Changing directory and reading methylation input files")
# change to the input file directory and find the bismark-formatted files
setwd(opt$dir)
file_list <- list.files(pattern = ".bismark.cov")

message("file list to read:")
print(file_list)

# read in the bismark-formatted files and set their attributes
bsobj.raw <- read.bismark(files=file_list)
message("completed reading bismark files")
bsobj <- orderBSseq(bsobj.raw) # need to make sure it's ordered for the smoothing function

samples=sub("\\.bismark.cov","",file_list)
sampleNames(bsobj) <- samples
message("sample names derived from file names:")
print(samples)
bsseq::pData(bsobj) <- data.frame(status = ifelse(samples %in% affected, "affected", "unaffected"))

message("running smoothing function")
bsobj.smooth <- BSmooth(bsobj, BPPARAM= MulticoreParam(workers = opt$threads, progressbar = TRUE), verbose=T)
message("completed smoothing")

sampleNames(bsobj.smooth)<- samples
BS.cov <- getCoverage(bsobj)
dmrs_df <- data.frame()

# switch working directory to the output directory
setwd(opt$outputdir)

# save the smoothed data for followup analysis
message("saving smoothed data to file")
saveRDS(bsobj.smooth,file="smoothed_bsobj.rds")
