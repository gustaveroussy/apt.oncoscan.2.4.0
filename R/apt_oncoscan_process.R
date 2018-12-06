## Performs OS CEL pairs processing
apt.oncoscan.process <- function(ATChannelCel = NULL, GCChannelCel = NULL, samplename = NULL, dual.norm = TRUE, out.dir = getwd(), temp.files.keep = FALSE, force.OS = NULL, apt.build = "na33.r2") {

  # ## TEMP
  # setwd("/home/job/svn/genomics/CGH/R/00_PIPELINE/TEST_ZONE/OS/")
  # ATChannelCel = "17P00005_A.CEL"
  # GCChannelCel = "17P00005_C.CEL"
  # samplename = "17P00005"
  # dual.norm = TRUE
  # out.dir = getwd()
  # temp.files.keep = FALSE
  # force.OS = NULL
  # apt.build = "na33.r2"

  if (any(is.null(c(ATChannelCel, GCChannelCel)))) stop("A CEL file is required !")
  if (is.null(samplename)) stop("A samplename is required !")
  if(!file.exists(ATChannelCel)) stop(paste0("Could not find ATChannelCel file ", ATChannelCel, "!"))
  if(!file.exists(GCChannelCel)) stop(paste0("Could not find GCChannelCel file ", GCChannelCel, "!"))
  if (ATChannelCel == GCChannelCel) stop("ATChannelCel and GCChannelCel files identical !")
  if (!dir.exists(out.dir)) stop(paste0("Output directory [", out.dir, "] does not exist !"))
  
  out.dir <- tools::file_path_as_absolute(out.dir)
  ATChannelCel <- tools::file_path_as_absolute(ATChannelCel)
  GCChannelCel <- tools::file_path_as_absolute(GCChannelCel)
  
  
  ## Checking build compatibility
  knownbuilds <- c("na33.r2", "na33.r4", "na36.r1")
  if (!(tolower(apt.build) %in% knownbuilds)) warning(paste0(" WARNING : The requested build ", apt.build, " is not in the validated list. Program may crash / fail !"))
  
  ## Checking apt-copynumber-cyto-ssa package loc
  tool.version <- "2.4.0"
  self.pkg.name <- paste0("apt.oncoscan.", tool.version)
  bin.dir <- system.file("apt/bin/", package = self.pkg.name)
  
  ## Checking apt-copynumber-onco annotation package loc
  sup.arrays <- c("OncoScan", "OncoScan_CNV")
  chptA <- affxparser::readCelHeader(filename = ATChannelCel)$chiptype
  chptC <- affxparser::readCelHeader(filename = GCChannelCel)$chiptype
  if(!(chptA %in% sup.arrays)) stop(paste0("Array type not supported for ", ATChannelCel, " !"))
  if(!(chptC %in% sup.arrays)) stop(paste0("Array type not supported for ", GCChannelCel, " !"))
  if (chptA != chptC) stop("ATChannelCel and CGChannelCel files correspond to different array types !")
  if (chptA == "OncoScan") pkg.root <- "OncoScan" else if (chptA == "OncoScan_CNV") pkg.root <- "OncoScanCNV"
  res.pkg.name <- paste0(pkg.root, ".", tolower(apt.build))
  if (!(res.pkg.name %in% utils::installed.packages())) stop(paste0("Package ", res.pkg.name, " not found !"))
  res.dir <- system.file("apt/res/", package = res.pkg.name)
  suppressPackageStartupMessages(require(res.pkg.name, character.only = TRUE))
  apt.files <- annotation.set.describe()
  
  ## Checking annotation files availability
  for (f in names(apt.files)) { if (!file.exists(paste0(res.dir, "/", apt.files[[f]]))) stop(paste0("File ", apt.files[[f]], " is not available for ", apt.build, " !")) }
  
  ## Checking the OS
  os.list <- c("linux", "windows", "osx")
  my.os <- get.os()
  tmsg(paste0("OS is reported as : ", my.os))
  if (!is.null(force.OS)) {
    if (!(force.OS %in% os.list)) stop("Specified forced OS is not supported !")
    my.os <- force.OS
    tmsg(paste0("WARNING : Forcing OS to : ", my.os))
  } else if (!(my.os %in% os.list)) stop(paste0("Current OS [", my.os, "] not supported ! If you are sure of your OS support, use force.OS option with any of 'linux', 'windows', 'osx'"))
  
  if (my.os == "windows") my.os <- paste0(my.os, ".exe")
  
  oridir <- getwd()
  out.dir.p <- paste0(out.dir, "/", samplename)
  out.dir.w <- paste0(out.dir.p, "/temp")
  dir.create(path = out.dir.w, recursive = TRUE)
  setwd(out.dir.p)
  
  pfile <- paste0(out.dir.p, "/", samplename, "_pairs.txt")
  write.table(data.frame(ATChannelCel = ATChannelCel, GCChannelCel = GCChannelCel, stringsAsFactors = FALSE), file = pfile, quote = FALSE, sep = "\t", row.names = FALSE)

  apt.cmd <- c(paste0(bin.dir, paste0("/apt-copynumber-onco-ssa_", my.os, " ")),
               paste0("--log-file ", out.dir.w, "/apt-copynumber-onco-ssa.log "),
               paste0("--out-dir ", out.dir.w, " "),
               paste0("--temp-dir ", out.dir.w, "/temp "),
               "--force false ",
               "--keep-intermediate-data false ",
               "--keep-combined-cels false ",
               paste0("--snp-qc-snp-list ", res.dir, "/", apt.files$snplist, " "),
               paste0("--reference-file ", res.dir, "/", apt.files$refmodel, " "),
               "--y-target -0.8 ",
               "--nddqc-min-count 2000 ",
               paste0("--set-analysis-name ", samplename, " "),
               "--brlmmp-HARD 3 ",
               "--brlmmp-SB 0.750000 ",
               "--brlmmp-CM 2 ",
               "--brlmmp-bins 100 ",
               "--brlmmp-mix 1 ",
               "--brlmmp-bic 2 ",
               "--brlmmp-CSepPen 0.1 ",
               "--brlmmp-CSepThr 16 ",
               "--brlmmp-lambda 1.0 ",
               "--brlmmp-wobble 0.20 ",
               "--brlmmp-copyqc 0.000000 ",
               "--brlmmp-copytype -1 ",
               "--brlmmp-ocean 0.0000000000000001 ",
               "--brlmmp-clustertype 2 ",
               "--brlmmp-transform mva ",
               "--brlmmp-MS 0.06 ",
               paste0("--annotation-file ", res.dir, "/", apt.files$annotdb, " "),
               paste0("--x-probes-file ", res.dir, "/", apt.files$xprobes, " "),
               paste0("--y-probes-file ", res.dir, "/", apt.files$yprobes, " "),
               paste0("--cel-pairs-file ", pfile, " "),
               "--alpha-cn-calibrate 0.780931 ",
               "--alpha-X-cn-calibrate 0.800507 ",
               "--alpha-Y-cn-calibrate 0.800507 ",
               "--beta-cn-calibrate 1 ",
               "--beta-X-cn-calibrate 1 ",
               "--beta-Y-cn-calibrate 1 ",
               "--agr-denominator-source SNP ",
               "--agr-denominator-percentile 0.75 ",
               "--igender-female-threshold 0.44 ",
               "--igender-male-threshold 0.63 ",
               "--nddicov-min-nd-probe-count-per-bin 30 ",
               "--nddicov-use-nd-signal-markers off ",
               paste0("--dual-channel-normalization ", tolower(as.character(dual.norm)), " "),
               "--ad-outlier-trim 3.0 ",
               "--normal-diploid-detection-node:enable=true ",
               "--ndd-weights-type FLD ",
               "--ndd-point-count 128 ",
               "--ndd-bandwidth 0.25 ",
               "--ndd-cutoff 0.05 ",
               "--ndd-step 40 ",
               "--ndd-window 275 ",
               "--ndd-min-markers 10 ",
               "--ndd-log2-threshold 0.28 ",
               "--ndd-th1 0,0.4 ",
               "--ndd-th2 0.7,1.3 ",
               "--ndd-th3 0.7,1.3 ",
               "--ndd-th4 0,0.4 ",
               "--ndd-min-percent-nd 0 ",
               "--ndd-min-number-nd 2000 ",
               "--ndd-med-l2r-low-percentile-ND 0.015 ",
               "--ndd-med-l2r-threshold-ND 0.20 ",
               "--ndd-iqr-factor-coarse 5.0 ",
               "--ndd-iqr-factor-fine 3.0 ",
               "--sig-gender-reference-chromosome 2 ",
               "--sig-gender-xx-cutoff 0.85 ",
               "--sig-gender-xx-cutoff-high 1.22 ",
               "--sig-gender-y-cutoff 0.31 ",
               "--hpf-mini-block-rows 8 ",
               "--hpf-mini-block-cols 8 ",
               "--hpf-global-smooth-weight 256 ",
               "--hpf-local-smooth-weight 64 ",
               "--hpf-converged 0.000100 ",
               "--nddlrcov-min-nd-marker-count-per-bin 30 ",
               "--nddlrcov-use-nd-lr-markers on ",
               "--log2-ratio-adjuster-wave-node:enable=true ",
               "--wave-bandwidth 101 ",
               "--wave-bin-count 25 ",
               "--wave-count 6 ",
               "--wave-smooth true ",
               "--median-autosome-median-normalization true ",
               "--shrinkage-v2-node:enable=true ",
               "--shrink-lprec 0.25,2.0 ",
               "--outlier-shrink-lprec 0.25,2.0,8.0 ",
               "--shrink-converge 0.000001 ",
               "--outlier-shrink-converge 0.1 ",
               "--shrink-downweight-maxiter 25 ",
               "--outlier-shrink-large-stdev 125.0 ",
               "--cn-gender-cutoff 0.5 ",
               "--cn-state-threshold 0.0 ",
               "--baf-transform linear ",
               "--baf-trim 10.0 ",
               "--kernel-sigma-span 50 ",
               "--sig-snp-min-call-rate 0.6 ",
               "--sig-snp-min-compare-rate 0.65 ",
               "--sig-snp-min-concordance 0.60 ",
               "--sig-snp-min-strength-sketch-rank 0.305 ",
               paste0("--sig-snp-report-filename ", out.dir.w, "/", chptA, ".CelPairCheckReport.txt "),
               "--ap-ad-step 20 ",
               "--ap-ad-window 100 ",
               "--ap-ad-point-count 128 ",
               "--ap-ad-bandwidth 0.250000 ",
               "--ap-ad-cutoff 0.050000 ",
               "--ap-ad-threshold 0.350000 ",
               "--ap-ad-height-threshold 0.000000 ",
               "--ap-ad-height-threshold-bound 1024.000000 ",
               "--ap-ad-symmetry true ",
               "--ap-baf-step 20 ",
               "--ap-baf-window 100 ",
               "--ap-baf-point-count 129 ",
               "--ap-baf-bandwidth 0.2 ",
               "--ap-baf-cutoff 0.00000 ",
               "--ap-baf-threshold 0.30000 ",
               "--ap-baf-height-threshold 3.000000 ",
               "--ap-baf-height-threshold-bound 1.476981 ",
               "--ap-baf-symmetry true ",
               "--loh-error-rate 0.05 ",
               "--loh-beta 0.001 ",
               "--loh-alpha 0.01 ",
               "--loh-separation 1000000 ",
               "--loh-min-marker-count 10 ",
               "--loh-no-call-threshold 0.05 ",
               "--loh-min-genomic-span 2500000 ",
               "--loh-ss-thresh 3.59 ",
               "--loh-ss-begin-slope 1.6416 ",
               "--loh-baf-upper-thresh 0.5 ",
               "--loh-baf-lower-thresh 0.5 ",
               "--loh-slope-thresh 0.25 ",
               "--loh-alpha-cn-calibrate 0.780931 ",
               "--loh-alpha-X-cn-calibrate 0.800507 ",
               "--loh-alpha-Y-cn-calibrate 0.800507 ",
               "--loh-beta-cn-calibrate 1 ",
               "--loh-beta-X-cn-calibrate 1 ",
               "--loh-beta-Y-cn-calibrate 1 ",
               "--cn-neutral-loh-node:enable=false ",
               "--tuscan-m3-low-cost-band 0.030000 ",
               "--tuscan-nd-snp-qc-threshold 19.000000 ",
               "--tuscan-mapd-qc-threshold 0.4000000 ",
               "--tuscan-m1-shift 0.000000 ",
               "--tuscan-m1-trim 0.010000 ",
               "--tuscan-m1-threshold 0.001050 ",
               "--tuscan-m2-lr-cost 0.500000 ",
               "--tuscan-m2-min-seg-size 151 ",
               "--tuscan-m2-min-alt-seg-size 19 ",
               "--tuscan-m2-min-fine-seg-size 19 ",
               "--tuscan-m2-exclude-chr 24,25 ",
               "--tuscan-consecutive-hom-count 100 ",
               "--tuscan-no-break-vicinity 19 ",
               "--lambda 50.000000 ",
               "--tuscan-m3-rho-min 0.200000 ",
               "--tuscan-m3-rho-max 1.000000 ",
               "--tuscan-m3-rho-N 100 ",
               "--tuscan-m3-psi-min 1.000000 ",
               "--tuscan-m3-psi-max 6.000000 ",
               "--tuscan-m3-psi-N 100 ",
               "--tuscan-m3-low-cost-weight 0.010000 ",
               "--tuscan-m3-max-percent-hom 85.000000 ",
               "--tuscan-m3-min-count-het 10 ",
               "--tuscan-m3-ploidy-factor 3.000000 ",
               "--tuscan-m3-ploidy-sigma 0.350000 ",
               "--alt-mode-bin-size 0.333333 ",
               "--nGoF-bin-size 0.500000 ",
               "--MedL2Rcutoff 0.080000 ",
               "--MedBAFcutoff 0.000000 ",
               "--min-nGoF-length 1001 ",
               "--nGoF4cutoff 1000.000000 ",
               "--gamma 0.800000 ",
               "--tuscan-round-ploidy true ",
               "--tuscan-m1-use-raw-baf false ",
               "--infMarkL2RMin 0.06 ",
               "--infMarkerThreshold 1000 ",
               "--MinChrYthreshold -3 ",
               "--MedL2Rcutoff-no-cn2 0.08 ",
               "--MedBAFcutoff-no-cn2 0 ",
               "--min-nGoF-length-no-cn2 51 ",
               "--nGoF4cutoff-no-cn2 5000 ",
               "--tuscan-fp-exponent 2.0 ",
               "--tuscan-fp-min-med-l2r 0.1 ",
               "--tuscan-fp-max-med-l2r 0.25 ",
               "--tuscan-use-fp-exponent-in-normal-mode false ",
               "--tuscan-cn-bin '[0,6)[6,10)[10,20)[20,50)[50,10000]' ",
               "--tuscan-cn-bin-values 'as-is-2,6-10,10-20,20-50,50+' ",
               "--tuscan-run-once false ",
               "--tuscan-maxCN2Weight 0.025 ",
               "--tuscan-maxNonCN2Increase 1 ",
               "--tuscan-FN_l2r_floor 0.095 ",
               "--tuscan-FN_exp 2.46 ",
               "--tuscan-FN_low_l2r_threshold 0.08 ",
               "--tuscan-FN_low_BAF_threshold 0.48 ",
               "--tuscan-REL_l2r_weight 0.5 ",
               "--tuscan-REL_CN2_gauss_P 1.1 ",
               "--tuscan-REL_CN2_gauss_sigma 0.35 ",
               "--tuscan-REL_CN2_gauss_mu 2 ",
               "--tuscan-REL_CN4_gauss_P 1.023 ",
               "--tuscan-REL_CN4_gauss_sigma 0.2 ",
               "--tuscan-REL_CN4_gauss_mu 4 ",
               "--tuscan-REL_threshold 0.1 ",
               "--tuscan-REL_threshold1 0.76 ",
               "--tuscan-REL_threshold2_hi 0.9 ",
               "--tuscan-REL_threshold2_low 0.7 ",
               "--tuscan-REL_maxMinimumIndex 2 ",
               "--tuscan-REL_noSearchIfminWeight true ",
               "--tuscan-infMarkTrimBand 0.05 ",
               "--tuscan-m3-ploidy4-factor 1.25 ",
               "--tuscan-m3-ploidy4-sigma 0.2 ",
               "--tuscan-m3-ploidy4-mu 4.0 ",
               "--tuscan-lowPloidy1 0 ",
               "--tuscan-lowPloidy2 1.76 ",
               "--tuscan-lowPloidyTB 1 ",
               "--tuscan-lowPloidy_nGoF4 200 ",
               "--tuscan-hiPloidy1 2.3 ",
               "--tuscan-hiPloidy2 2.5 ",
               "--tuscan-hiPloidyTB 1 ",
               "--tuscan-hiPloidy_nGoF4 200 ",
               "--tuscan-centromere-breakpoints 1,121535434,2,92326171,3,90504854,4,49660117,5,46405641,6,58830166,7,58054331,8,43838887,9,47367679,10,39254935,11,51644205,12,34856694,13,-1,14,-1,15,-1,16,35335801,17,22263006,18,15460898,19,24681782,20,26369569,21,-1,22,-1,24,58632012,25,10104553 ",
               "--tuscan-centromere-min-seg-size 19 ",
               "--tuscan-report-internal-ploidy true ",
               paste0("--report-file ", out.dir.w, "/QCReport.txt "),
               "--cychp-output false ",
               "--text-output false ",
               "--oschp-output true ",
               "--mangle-probeset-names false ",
               "--waviness-block-size 50 ",
               "--waviness-genomic-span 0 ")

  tmsg("Running APT ...")
  # cmdres <- system(command = paste0(apt.cmd, collapse = ""), intern = TRUE)
  cmdres <- try(system(command = paste0(apt.cmd, collapse = ""), intern = TRUE))
  
  oscf <- list.files(path = out.dir.w, pattern = "\\.oschp$", full.names = TRUE, recursive = FALSE, ignore.case = TRUE)
  
  if (!file.exists(oscf)) return(cmdres)
  
  logf <- list.files(path = out.dir.w, pattern = "\\.log$", full.names = TRUE, recursive = FALSE, ignore.case = TRUE)
  qcf <- list.files(path = out.dir.w, pattern = "QCReport\\.txt$", full.names = TRUE, recursive = FALSE, ignore.case = TRUE)
  pcf <- list.files(path = out.dir.w, pattern = "\\.CelPairCheckReport\\.txt$", full.names = TRUE, recursive = FALSE, ignore.case = TRUE)
  
  ## Renaming files
  tmsg("Renaming OSCHP ...")
  new.oscf <- paste0(out.dir.p, "/", samplename, "_", tool.version, "_", apt.build, ".oschp")
  new.logf <- paste0(out.dir.p, "/", samplename, "_", tool.version, "_", apt.build, ".log")
  new.qcf <- paste0(out.dir.p, "/", samplename, "_", tool.version, "_", apt.build, ".qc.txt")
  new.pcf <- paste0(out.dir.p, "/", samplename, "_", tool.version, "_", apt.build, ".paircheck.txt")
  file.rename(from = oscf[1], to = new.oscf)
  file.rename(from = logf[1], to = new.logf)
  file.rename(from = qcf[1], to = new.qcf)
  file.rename(from = pcf[1], to = new.pcf)
  
  ## Cleaning
  if(!temp.files.keep) {
    tmsg("Removing temporary files ...")
    # system(paste0("rm -Rf ", out.dir.w))
    unlink(x = out.dir.w, recursive = TRUE, force = TRUE)
  }
  setwd(oridir)
  
  tmsg("Done.")
  return(new.oscf)
}


apt.oncoscan.process.batch <- function(pairs.file = NULL, nthread = 1, cluster.type = "PSOCK", ...) {
  ## Checking the pairs.file
  if (!file.exists(pairs.file)) stop("Could not find pairs.file !")
  message("Reading and checking pairs.file ...")
  mypairs <- read.table(file = pairs.file, header = TRUE, sep="\t", check.names = FALSE, as.is = TRUE)
  head.ok <- c("ATChannelCel", "GCChannelCel", "SampleName")
  head.chk <- all(colnames(pairs.file) == head.ok)
  if (!head.chk) {
    message("Invalid header in pairs.file !")
    message(paste0("EXPECTED : ", head.ok))
    message(paste0("FOUND : ", colnames(mypairs)))
    stop("Invalid header.")
  }
  
  at.chk <- file.exists(mypairs$ATChannelCel)
  gc.chk <- file.exists(mypairs$GCChannelCel)
  if (!all(at.chk) || !all(at.chk)) {
    message("Some CEL files from the pairs.file could not be found (wrong path ?) !")
    message("Missing AT CEL files :")
    message(mypairs$ATChannelCel[which(!at.chk)])
    message("Missing GC CEL files :")
    message(mypairs$GCChannelCel[which(!gc.chk)])
    stop("Missing CEL.")
  }
  sn.chk <- duplicated(mypairs$SampleName)
  if (any(sn.chk)) {
    message("pairs.file contains duplicated SampleNames !")
    message(mypairs$SampleName[which(duplicated(mypairs$SampleName))])
    stop("Duplicated SampleNames.")
  }
  if(any(mypairs$ATChannelCel == mypairs$GCChannelCel)) {
    message("Some ATChannelCel and GCChannelCel are the same for at least one sample !")
    stop("Identical CEL files for AT and GC.")
  }
  
  ## Adjusting cores/threads
  message("Adjusting number of cores if needed ...")
  if (is.null(nthread)) nthread <- parallel::detectCores(logical = TRUE) -1
  if (nrow(mypairs) < nthread) nthread <- nrow(mypairs)
  
  message("Running APT ...")
  `%dopar%` <- foreach::"%dopar%"
  cl <- parallel::makeCluster(spec = nthread, type = cluster.type, outfile = "")
  doParallel::registerDoParallel(cl)
  
  p <- 0
  osres <- foreach(p = seq_len(nrow(mypairs)), .inorder = FALSE, .errorhandling = "stop") %dopar% {
    apt.oncoscan.process(ATChannelCel = mypairs$ATChannelCel[p], GCChannelCel = mypairs$GCChannelCel[p], samplename = mypairs$SampleName[p], ...)
  }
  message("Stopping cluster ...")
  parallel::stopCluster(cl)
  message("Done.")
}

## Print thread-tagged message
tmsg <- function(text = NULL) { message(paste0(" [", Sys.info()[['nodename']], ":", Sys.getpid(), "] ", text)) }

## A more robust way to get machine OS type
get.os <- function(){
  sysinf <- Sys.info()
  if (!is.null(sysinf)){
    os <- sysinf['sysname']
    if (os == 'Darwin')
      os <- "osx"
  } else { ## mystery machine
    os <- .Platform$OS.type
    if (grepl("^darwin", R.version$os))
      os <- "osx"
    if (grepl("linux-gnu", R.version$os))
      os <- "linux"
  }
  return(tolower(os))
}
