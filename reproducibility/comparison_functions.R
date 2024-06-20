# if (!require("remotes")) {
#   install.packages("remotes")
# }
# if (!require("xbioc")) {
#   remotes::install_github("renozao/xbioc")
# }
# if (!require("SCDC")) {
#   remotes::install_github("meichendong/SCDC")
# }

library(reshape)

if (!file.exists('config.csv')) {
  stop('config.csv must be created with columns "username" and "token"')
}

# https://doi.org/10.1093/bib/bbz166
# https://meichendong.github.io/SCDC/
SCDC <- function(sce, groupby, bulk_df, sample_col="") {
  
  if (sample_col=="") {
    sce$sample = 'sample'
    sample_col='sample'
  }

  rlang::check_installed("SCDC", reason = "to run SCDC deconvolution.")
  
  # In an ExpressionSet, assayData, columns are samples and rows are genes
  bulk_eset = Biobase::ExpressionSet(assayData=as.matrix(bulk_df))
  
  sc_metadata = as.data.frame(cbind(sce@colData[groupby], sce@colData[sample_col]))
  sc_eset = Biobase::ExpressionSet(assayData=as.matrix(normcounts(sce)), phenoData = as(sc_metadata, "AnnotatedDataFrame"))
                                   
 # Removed ct.sub = c("alpha","beta","delta","gamma","ductal","acinar")
  celltypes = unique(sce@colData[[groupby]])
  samples = unique(sce@colData[[sample_col]])
  print(celltypes)
  print(samples)
  
  if (length(samples) == 1) {
    sc_qc = SCDC::SCDC_qc_ONE(sc_eset, ct.varname = groupby, sample = sample_col, qcthreshold = 0.7, ct.sub=celltypes, generate.figure=FALSE)
    results <- SCDC::SCDC_prop_ONE(
      bulk.eset = bulk_eset, 
      sc.eset = sc_qc$sc.eset.qc, 
      ct.varname = groupby, 
      ct.sub=celltypes,
      sample = sample_col)
    return(results)
  } else {
    sc_qc = SCDC::SCDC_qc(sc_eset, ct.varname = groupby, sample = sample_col, qcthreshold = 0.7, ct.sub=celltypes, generate.figure=FALSE)
    results <- SCDC::SCDC_prop(
      bulk.eset = bulk_eset, 
      sc.eset = sc_qc$sc.eset.qc, 
      ct.varname = groupby, 
      ct.sub=celltypes,
      sample = sample_col)
    return(results)
  }
                                   
}

Tonkin_Hill_Mixture_Model <- function(sce, groupby, bulk_df) {
  
}

MUSIC <- function(sce, groupby, bulk_df) {
  
}

CibersortX <- function(sce, groupby, bulk_df, n_threads = 1) {
  
  config = read.csv('config.csv')
  cibersortx_username = config$username
  cibersortx_token = config$token
  
  rownames(bulk_df) = gsub(x=rownames(bulk_df), pattern="-", replacement = "")
  rownames(bulk_df) = gsub(x=rownames(bulk_df), pattern="_", replacement = "")
  
  colnames(bulk_df) = gsub(x=colnames(bulk_df), pattern="-", replacement = "")
  colnames(bulk_df) = gsub(x=colnames(bulk_df), pattern="_", replacement = "")
  
  bulk_df = cbind(genes = rownames(bulk_df), bulk_df, row.names = NULL)
  
  cibersort_sc = normcounts(sce)
  cibersort_sc = as.data.frame(as.matrix(cibersort_sc))
  rownames(cibersort_sc) = gsub(x=rownames(cibersort_sc), pattern="-", replacement = "")
  rownames(cibersort_sc) = gsub(x=rownames(cibersort_sc), pattern="_", replacement = "")
  colnames(cibersort_sc) = gsub(x=colnames(cibersort_sc), pattern="-", replacement = "")
  colnames(cibersort_sc) = gsub(x=colnames(cibersort_sc), pattern="_", replacement = "")
  colnames(cibersort_sc) = paste0("c", sce@colData[[groupby]])
  cibersort_sc = cbind(data.frame(genes = rownames(cibersort_sc)), cibersort_sc, row.names = NULL)
  
  max_genes_to_use = min(500, length(sce@rowRanges)-100)
  
  if (!all(colnames(normcounts(sce)) == rownames(sce@colData))) {
    stop('gene names do not match')
  }
  
  # create new directory
  temp_dir = fs::path('temp_blase_cibersortx')
  print(temp_dir)
  
  if (dir.exists(temp_dir) == FALSE) {
    dir.create(temp_dir, recursive = TRUE)
  }
  
  original_wd = getwd()
  setwd(temp_dir)
  
  # save files to wd (as 'reference.txt' and 'bulk.txt')
  write.table(bulk_df, file='bulk.txt', quote=FALSE, sep='\t', col.names = TRUE, row.names=FALSE)
  write.table(cibersort_sc, file='reference.txt', quote=FALSE, sep='\t', col.names = TRUE, row.names=FALSE)
  
  # submit to docker machine here
  print(paste('Running: docker', "pull cibersortx/fractions"))
  res = system2("docker", args="pull cibersortx/fractions")
  print(res)
  
  cwd = getwd()
  get_fractions_args = paste0('run --platform linux/amd64 -v ', cwd, ':/src/data -v ', cwd,':/src/outdir cibersortx/fractions --single_cell TRUE --username ', cibersortx_username, ' --token ', cibersortx_token, ' --refsample reference.txt --mixture bulk.txt --verbose TRUE --perm 200 --G.max ', max_genes_to_use)
  print(paste('Running: docker', get_fractions_args))
  res = system2('docker', args=get_fractions_args)
  print(res)
  
  # read in results 
  result_path = fs::path('CIBERSORTx_Results', ext = 'txt')
  cibersort_results = read.table(result_path, header=TRUE)
  setwd(original_wd)
  
  # tidy up cibersortx directory
  fs::dir_delete(temp_dir)

  # return results
  return(cibersort_results)
  
}
