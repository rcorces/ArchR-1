# Bulk Tracks Methods -----------------------------------------------------------

## Create fragment files -----------------------------------------------------------

keepFilteredChromosomes <- function (x, remove = c("chrM"), underscore = TRUE, standard = TRUE, pruning.mode = "coarse") {
  if (standard) {
    x <- GenomeInfoDb::keepStandardChromosomes(x, pruning.mode = pruning.mode)
  }
  seqNames <- seqlevels(x)
  chrRemove <- c()
  if (underscore) {
    chrRemove <- c(chrRemove, which(grepl("_", seqNames)))
  }
  chrRemove <- c(chrRemove, which(seqNames %in% remove))
  if (length(chrRemove) > 0) {
    chrKeep <- seqNames[-chrRemove]
  }
  else {
    chrKeep <- seqNames
  }
  seqlevels(x, pruning.mode = pruning.mode) <- chrKeep
  return(x)
}

.getGroupFragsFromProj <- function(
  ArchRProj = NULL,
  groupBy = NULL,
  outDir = "Shiny/fragments/"
){
  
  dir.create(outDir, showWarnings = FALSE) 
  
  # find barcodes of cells in that groupBy.
  groups <- getCellColData(ArchRProj, select = groupBy, drop = TRUE)
  cells <- ArchRProj$cellNames
  cellGroups <- split(cells, groups)
  
  # outputs unique cell groups/clusters. 
  clusters <- names(cellGroups)
  
  for (cluster in clusters){
    cat("Making fragment file for cluster:", cluster,"\n")
    # get GRanges with all fragments for that cluster 
    cellNames = cellGroups[[cluster]]
    fragments <- getFragmentsFromProject(ArchRProj = ArchRProj, cellNames = cellNames)
    fragments <- unlist(fragments, use.names = FALSE) 
    # filter Fragments
    ATACFragments <- keepFilteredChromosomes(fragments)
    saveRDS(ATACFragments, paste0(outDir, cluster, "_fragments.rds"))
  }
}

## Create coverage objects -----------------------------------------------------------

addSeqLengths <- function (gr, genome){
  gr <- validGRanges(gr)
  genome <- validBSgenome(genome)
  stopifnot(all(as.character(seqnames(gr)) %in% as.character(seqnames(genome))))
  seqlengths(gr) <- seqlengths(genome)[as.character(names(seqlengths(gr)))]
  return(gr)
}

validGRanges <- function(gr = NULL){
  stopifnot(!is.null(gr))
  if(inherits(gr, "GRanges")){
    return(gr)
  }else{
    stop("Error cannot validate genomic range!")
  }
}

.getClusterCoverage <- function(
  ArchRProj = NULL,
  tileSize = 100, 
  groupBy = "Clusters",
  outDir = "Shiny/coverage/"
  # geneAnnotation = getGeneAnnotation(ArchRProj),
){
  fragfiles = list.files(path = "./Shiny/fragments", full.names = TRUE)
  dir.create(outDir, showWarnings = FALSE) 
  
  # find barcodes of cells in that groupBy.
  groups <- getCellColData(ArchRProj, select = groupBy, drop = TRUE)
  cells <- ArchRProj$cellNames
  cellGroups <- split(cells, groups)
  
  # outputs unique cell groups/clusters. 
  clusters <- names(cellGroups)
  
  chrRegions <- getChromSizes(ArchRProj)
  genome <- getGenome(ArchRProj)
  
  clusteridx=0
  for(file in fragfiles){
    clusteridx = clusteridx + 1
    ATACFragments <- readRDS(file) 
    
    #fragmentsToInsertions()
    left <- GRanges(seqnames = seqnames(ATACFragments), 
                    ranges = IRanges(start(ATACFragments), width = 1))
    right <- GRanges(seqnames = seqnames(ATACFragments), 
                     ranges = IRanges(end(ATACFragments), width = 1))
    # call sort() after sortSeqlevels() to sort also the ranges in addition 
    # to the chromosomes. 
    insertions <- c(left, right) %>% sortSeqlevels() %>% 
      sort()
    
    #binnedCoverage
    message("creating bins for cluster ",clusters[clusteridx], "...")
    bins <- unlist(slidingWindows(chrRegions, width = tileSize, step = tileSize))
    message("counting overlaps for cluster ",clusters[clusteridx], "...")
    bins$reads <- countOverlaps(bins, insertions, maxgap = -1L, minoverlap = 0L, type = "any")
    addSeqLengths(bins, genome)
    message("creating binned coverage for cluster ",clusters[clusteridx], "...")
    #each value is multiplied by that weight.
    binnedCoverage <- coverage(bins, weight = bins$reads)
    saveRDS(binnedCoverage, paste0(outDir,clusters[clusteridx], "_cvg.rds"))
  }
}


#' Export a Shiny App based on ArchRProj
#' 
#' Generate all files required for an autonomous Shiny app to display your browser tracks.
#'
#' @param ArchRProj An `ArchRProject` object loaded in the environment. Can do this using: loadArchRProject("path to ArchRProject/")
#' @param threads The number of threads to use for parallel execution.
#' @param verbose A boolean value that determines whether standard output should be printed.
#' @param logFile The path to a file to be used for logging ArchR output.
#' @export
exportShinyArchR <- function(
  ArchRProj = NULL,
  outDir = "Shiny",
  #where Shiny files are now. OFFLINE. 
  #myDir = 
  groupBy = "Clusters",
  tileSize = 100,
  threads = getArchRThreads(),
  verbose = TRUE,
  logFile = createLogFile("exportShinyArchR")
){
  
  .validInput(input = ArchRProj, name = "ArchRProj", valid = c("ArchRProj"))
  # .validInput(input = features, name = "features", valid = c("granges", "grangeslist", "null"))
  # .validInput(input = loops, name = "loops", valid = c("granges", "grangeslist", "null"))
  # .validInput(input = minCells, name = "minCells", valid = c("integer"))
  # .validInput(input = baseSize, name = "baseSize", valid = c("integer"))
  # .validInput(input = borderWidth, name = "borderWidth", valid = c("numeric"))
  # .validInput(input = tickWidth, name = "tickWidth", valid = c("numeric"))
  # .validInput(input = facetbaseSize, name = "facetbaseSize", valid = c("numeric"))
  # geneAnnotation <- .validGeneAnnotation(geneAnnotation)
  # .validInput(input = browserTheme, name = "browserTheme", valid = c("character"))
  .validInput(input = threads, name = "threads", valid = c("integer"))
  .validInput(input = verbose, name = "verbose", valid = c("boolean"))
  .validInput(input = logFile, name = "logFile", valid = c("character"))
  
  .startLogging(logFile=logFile)
  .logThis(mget(names(formals()),sys.frame(sys.nframe())), "exportShinyArchR Input-Parameters", logFile = logFile)
  
  .requirePackage("shiny", installInfo = 'install.packages("shiny")')
  .requirePackage("rhandsontable", installInfo = 'install.packages("rhandsontable")')

# Make directory for Shiny App 
  if(!dir.exists(outDir)) {
    
    dir.create(outDir)
  # if(length(dir(outDir,  all.files = TRUE, include.dirs = TRUE, no.. = TRUE)) > 0){
  #   stop("Please specify a new or empty directory")
  # }

  filesUrl <- c(
    "https://raw.githubusercontent.com/paupaiz/ArchR/Shiny_export/R/Shiny/app.R"
    "https://raw.githubusercontent.com/paupaiz/ArchR/Shiny_export/R/Shiny/global.R"
    "https://raw.githubusercontent.com/paupaiz/ArchR/Shiny_export/R/Shiny/server.R"
    "https://raw.githubusercontent.com/paupaiz/ArchR/Shiny_export/R/Shiny/ui.R"
  )
  
  downloadFiles <- lapply(seq_along(filesUrl), function(x){
    download.file(
      url = filesUrl[x], 
      destfile = file.path(outDir, basename(filesUrl[x]))
    )
    })
  
  }else{
    message("Using existing Shiny files...")
  }

  # Create a copy of the ArchRProj object
  ArchRProjShiny <- ArchRProj
  # Add metadata to ArchRProjShiny
  
  if (is.na(paste0("ArchRProj$", groupBy))) {
    stop("groupBy is not part of cellColData")
  } else if ((any(is.na(paste0("ArchRProj$", groupBy))))) {
    stop("incomplete data. some NA observations for groupBy")
  } else {
    ArchRProjShiny@projectMetadata[["groupBy"]] <- groupBy
  }
  
  ArchRProjShiny@projectMetadata[["tileSize"]] <- tileSize
  
  saveArchRProject(ArchRProj = ArchRProjShiny, outputDirectory = "Save-ArchRProjShiny")
  
# Create fragment files 
.getGroupFragsFromProj(ArchRProj = ArchRProjShiny, groupBy = groupBy)

# Create coverage objects
.getClusterCoverage(ArchRProj = ArchRProjShiny, tileSize = tileSize, groupBy = groupBy)

  ## ready to launch ---------------------------------------------------------------
  message("App created! To launch, 
          ArchRProjShiny <- loadArchRProject('path to ArchRProject/') and 
          run shiny::runApp('", outDir, "') from parent directory")
#  runApp("myappdir")
}
