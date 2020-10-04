validate_metadata <- function(metadata="thesauri/metadata.xlsx") {
  
  metadata <- read_excel(metadata)
  total_df <- metadata[-c(1:4),] #remove references, basisOfRecord and coordinates info
  id_df <- total_df %>% filter(category == "trait") %>% select(traitName = names, identifier)
  total_df[,c("identifier", "relatedTerm", "broaderTerm", "narrowerTerm", "comments")] <- NULL #remove identifier info
  datasets <- list.files("raw_datasets", pattern = ".csv$", full.names = TRUE) #list with file path and names
  dataset_names <- gsub("raw_datasets/|.csv", "", datasets) #names of all datasets
  metadata_df <- metadata[1:4,] %>% select(dataset_names)
  '%notin%' <- Negate('%in%')
  clean = TRUE
  
  for (i in 1:length(datasets)) {
    
    # Load thesaurus
    loop_thesaurus <- total_df %>% select(category, verbatim = dataset_names[i]) %>% drop_na()
    
    # Check if column headers from the form and dataset match
    raw_dataset <- fread(input = datasets[i])
    
    stop = FALSE
    for (k in loop_thesaurus$verbatim) {
      if (k %notin% colnames(raw_dataset)) {
        message(paste0('"',k, '" was not found as column header in "', dataset_names[i], '".'))
        stop = TRUE
        clean = FALSE
      }
    }
    if (stop) {next()}
    
    # Check for non-numeric values in trait columns
    dataset <- fread(input = datasets[i], select = loop_thesaurus$verbatim)
    traits <- loop_thesaurus %>% filter(category == "trait")
    trait_cols <- dataset %>% select(traits$verbatim)
    trait_names <- as.vector(names(trait_cols))
    
    for (j in colnames(trait_cols)) {
      num_check <- class(trait_cols[[j]])
      if (num_check %in% c("numeric","integer") == FALSE) {
        message(paste0('The column "', j, '" in the dataset "', dataset_names[i],'" contains non-numeric values.'))
        clean = FALSE
      }
    }    
  }
  if (clean) {
    message("The datasets and metadata thesaurus are validated.")
    } else {
      message("\nThe mistakes above have to be fixed first before resuming the workflow.")
    }
}