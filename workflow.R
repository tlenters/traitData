#----Install and load packages----
ls <- c("readxl","data.table","dplyr","reshape2","traitdataform","stringr","rgdal","Taxonstand","tidyverse","rangeBuilder")
new.packages <- ls[!(ls %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library(readxl)
library(data.table)
library(dplyr)
library(reshape2)
library(traitdataform)
library(stringr)
library(rgdal)
library(Taxonstand)
library(tidyverse)

setwd("~/data_repository")


#----Load metadata and units thesauri----
metadata <- read_excel("thesauri/metadata.xlsx")
units <- read_excel("thesauri/units.xlsx")


#---Data integration----
total_df <- metadata[-c(1:4),] #remove references, basisOfRecord and coordinates info
id_df <- total_df %>% filter(category == "trait") %>% select(traitName = names, identifier) #df with trait names and URIs
total_df$identifier <- NULL #remove identifier info
unit_df <- units[-1,] #remove references

datasets <- list.files("raw_datasets", pattern = ".csv$", full.names = TRUE) #list with file path and names
dataset_names <- gsub("raw_datasets/|.csv", "", datasets) #names of all datasets
metadata_df <- metadata[1:4,] %>% select(dataset_names) #df with metadata information

ext_df <- data.frame() #empty df used for row binding

for (i in 1:length(datasets)){
  
  #thesaurus subsets
  loop_thesaurus <- total_df %>% select(category, classes, names, verbatim = dataset_names[i]) %>% drop_na()
  numeric_names <- loop_thesaurus %>% filter(classes == "numeric")
  character_names <- loop_thesaurus %>% filter(classes == "character")
  
  #read classes from dataset
  numeric_df <- fread(input = datasets[i], 
                      select = numeric_names$verbatim, 
                      col.names = numeric_names$names,
                      colClasses = "numeric",
                      encoding = "UTF-8")
  character_df <- fread(input = datasets[i], 
                        select = character_names$verbatim,
                        col.names = character_names$names,
                        colClasses = "character",
                        encoding = "UTF-8")
  
  #bind for full dataset with standardized names
  dataset <- cbind(numeric_df, character_df)
  
  #melt data_full to long table
  taxon_occ_meas <- loop_thesaurus %>% filter(category != "trait")
  traits <- loop_thesaurus %>% filter(category == "trait")
  data_melt <- melt(data = dataset, id.vars = taxon_occ_meas$names, measure.vars = traits$names, value.name = "traitValue", variable.name = "traitName", variable.factor = FALSE) %>% drop_na(traitValue)
  
  #add references, basisOfRecord and CRS
  data_melt$references <- metadata_df[1,i] #add reference as column
  data_melt$basisOfRecord <- metadata_df[2,i] #add basisOfRecord column
  data_melt$verbatimCoordinateSystem <- metadata_df[3,i] #add verbatimCoordinateSystem as column
  data_melt$verbatimSRS <- metadata_df[4,i] #add verbatimSRS as column
  
  #add unit column
  loop_unit <- unit_df %>% select(traitName = names, traitUnit = dataset_names[i]) %>% drop_na() #relate unit to traitName
  data_melt <- full_join(data_melt, loop_unit, by = "traitName") #add unit column by merging
  
  #standardize measurement units
  data_melt$traitValue <- ifelse(data_melt$traitUnit == "m", data_melt$traitValue*100, data_melt$traitValue)
  data_melt$traitValue <- ifelse(data_melt$traitUnit == "mm", data_melt$traitValue/10, data_melt$traitValue)
  data_melt$traitUnit <- ifelse(data_melt$traitUnit == "m" | data_melt$traitUnit == "mm", "cm", data_melt$traitUnit)
  
  #add verbatim traitname column
  verba_name <- select(traits, traitName = names, verbatimTraitName = verbatim)  #df to relate traitname to verbatim traitname
  data_melt <- full_join(data_melt, verba_name, by = "traitName") #add verbatim traitname column by merging
  
  #add trait-identifier column
  data_melt <- full_join(data_melt, id_df, by = "traitName")
  data_melt <- subset(data_melt, is.na(traitValue) == FALSE)
  
  #output df (used for core table and extensions)
  ext_df <- bind_rows(ext_df, data_melt)
  
}


#----Taxon extension----
ext_df$scientificName <- ifelse(ext_df$scientificName %in% NA, gsub("[[:punct:]]+", "", paste(ext_df$genus, ext_df$specificEpithet)), gsub("[[:punct:]]+", "", ext_df$scientificName)) #make scientific name if not present
ext_df$genus <- ifelse(ext_df$genus %in% NA, word(ext_df$scientificName, 1), ext_df$genus) #make genus name if not present
ext_df$specificEpithet <- ifelse(ext_df$specificEpithet %in% NA, word(ext_df$scientificName, 2), ext_df$specificEpithet) #make specific epithet if not present

taxa_std <- standardise_taxa(ext_df) #add taxonomic info columns
taxon_names <- subset(metadata, category == "taxon") #subset only columns with names of taxon category
taxon_names <- c(taxon_names$names) #vector with taxon names
taxon_subset <- select(taxa_std, one_of(taxon_names)) #only select columns with names in taxon vector
taxa_df <- taxa_std %>% select(taxonID = taxonID,
                               scientificNameStd = scientificNameStd, 
                               kingdom = kingdom,
                               phylum = phylum,
                               class = class,
                               order = order,
                               family = family) %>% mutate_each(list(as.character))
taxa_bind <- cbind(taxa_df, taxon_subset) #combine columns
Taxon <- distinct(taxa_bind) #only select unique rows
Taxon$taxonID2 <- seq(1:nrow(Taxon)) #add more specific taxonID

ext_df2 <- full_join(Taxon, taxa_std, by = c(colnames(Taxon[!colnames(Taxon) %in% "taxonID2"]))) #add taxonID2 column

tpl <- TPL(splist = unique(Taxon$scientificName), corr = TRUE) #add taxonomic info form The Plant List
tpl_merge <- left_join(Taxon, tpl, by = c("scientificName" = "Taxon")) #merge with Taxon df

temp_zip <- tempfile() #make temporary files to story WCVP zip and text files
temp_txt <- tempfile()
download.file(url = "http://sftp.kew.org/pub/data-repositories/WCVP/wcvp_v2_jun_2020.zip", destfile = temp_zip) #download WCVP database
unzip(zipfile = temp_zip, exdir = temp_txt) #unzip and save database as text file
WCVP <- fread(file.path(temp_txt, "wcvp_export.txt")) #read WCVP text file
unlink(c(temp_zip, temp_txt)) #remove temporary files

wcvp_sub <- subset(WCVP, family %in% Taxon$family & rank == "SPECIES") #subset WCVP database to only included family names present in the datasets and only species names
wcvp_sub$accepted <- ifelse(wcvp_sub$taxonomic_status == "Accepted", wcvp_sub$taxon_name, wcvp_sub$taxon_name) #make column for all accepted scientific names
wcvp_merge <- left_join(tpl_merge, wcvp_sub, by = c("scientificName" = "taxon_name")) #merge with tpl_merge

Taxon2 <- data.frame(taxonID = wcvp_merge$taxonID2, #df for Taxon extension 
                     verbatimScientificName = wcvp_merge$scientificName,
                     scientificNameGBIF = wcvp_merge$scientificNameStd,
                     scientificNameTPL = ifelse(wcvp_merge$Plant.Name.Index == "FALSE" & wcvp_merge$Higher.level == "FALSE", NA ,paste(wcvp_merge$New.Genus, wcvp_merge$New.Species)),
                     scientificNameWCVP = wcvp_merge$accepted,
                     verbatimInfraspecificEpithet = wcvp_merge$infraspecificEpithet,
                     kingdom = wcvp_merge$kingdom,
                     phylum = wcvp_merge$phylum,
                     class = wcvp_merge$class,
                     order = wcvp_merge$order,
                     family = wcvp_merge$family.x,
                     genus = wcvp_merge$genus.x,
                     stringsAsFactors = F,
                     wcvp_merge[setdiff(taxon_names, c("specificEpithet","infraspecificEpithet","scientificName","genus"))],
                     GBIFID = wcvp_merge$taxonID,
                     TPLID = ifelse(wcvp_merge$New.ID == "", NA,
                                    paste0("http://www.theplantlist.org/tpl",wcvp_merge$TPL.version,"/record/", wcvp_merge$New.ID)),
                     WCVPID = ifelse(wcvp_merge$accepted_kew_id == "", 
                                     paste0("http://powo.science.kew.org/taxon/urn:lsid:ipni.org:names:", wcvp_merge$kew_id), 
                                     paste0("http://powo.science.kew.org/taxon/urn:lsid:ipni.org:names:", wcvp_merge$accepted_kew_id)))
fwrite(Taxon2, "output_tables/Taxon_ext.csv", row.names = FALSE) #write taxon extension file


#----Measurement or Fact extension----
subset_meas <- subset(taxon_occ_meas, category == "measurement") #subset only columns with names of measurement category
vector_meas <- c(subset_meas$names) #vector with measurement names
Measurement <- data.frame(ext_df2[vector_meas], #df with all measurement or fact info
                          basisOfRecord = ext_df2$basisOfRecord,
                          references = ext_df2$references)
Measurement <- distinct(Measurement) #select only distinct rows
Measurement$measurementID <- seq(nrow(Measurement)) #add measurementID
Measurement <- Measurement %>% select(measurementID, everything()) #measurementID as first column
write.csv(Measurement, "output_tables/Measurement_or_Fact_ext.csv", fileEncoding = "Latin1", row.names = FALSE) #write measurement or fact extension file

ext_df3 <- full_join(Measurement, ext_df2, by = c(colnames(Measurement[!colnames(Measurement) %in% "measurementID"])))


#----Occurrence extension----
EPSG_df <- make_EPSG() #download EPSG database
EPSG_df$verbatimSRS <- gsub("# ","", EPSG_df$note) #remove symbols
EPSG_df[,c("note","prj4")] <- NULL #remove note and prj4 columns
EPSG_df <- EPSG_df[!duplicated(EPSG_df[,"verbatimSRS"]),] #remove duplicate verbatimSRS values
ext_df3 <- left_join(ext_df3, EPSG_df, by = "verbatimSRS") #merge by verbatimSRS
ext_df3$geodeticDatum <- ext_df3$code #rename EPSG code column
ext_df3$code <- NULL
ext_df3$verbatimSRS <- NULL
ext_df3$geodeticDatum <- ifelse(ext_df3$geodeticDatum %in% NA, "unknown", as.character(paste0("EPSG:",ext_df3$geodeticDatum))) #add "EPSG:" tag
ext_df3$country <- rangeBuilder::standardizeCountry(ext_df3$verbatimCountry, fuzzyDist = 5) #standardize country names

subset_occ <- subset(taxon_occ_meas, category == "occurrence") #subset only columns with names of occurrence category
vector_occ <- c(subset_occ$names) #vector with occurrence names
Occurrence <- data.frame(ext_df3[vector_occ],
                         country = ext_df3$country,
                         geodeticDatum = ext_df3$geodeticDatum, stringsAsFactors = FALSE) #df with all occurrence info
Occurrence <- distinct(Occurrence) #only select unique rows
Occurrence$occurrenceID <- seq(1:nrow(Occurrence)) #add occurrenceID
Occurrence <- Occurrence %>% select(occurrenceID, everything()) #move ID to first position
fwrite(Occurrence, "output_tables/Occurrence_ext.csv", row.names = FALSE) #write occurrence extension file

ext_df4 <- full_join(Occurrence, ext_df3, by = c(colnames(Occurrence[!colnames(Occurrence) %in% "occurrenceID"]))) #occurrenceID column


#----Core table----
core_table <- ext_df4 %>% select(scientificName = scientificNameStd,#df with all core values and ID's
                               verbatimScientificName = scientificName,
                               verbatimTraitName = verbatimTraitName,
                               traitName = traitName,
                               traitValue = traitValue,
                               traitUnit = traitUnit,
                               traitID = identifier,
                               taxonID = taxonID2,
                               measurementID = measurementID,
                               occurrenceID = occurrenceID)
fwrite(core_table, "output_tables/core_table.csv", row.names = FALSE) #write file with core values and ID's

