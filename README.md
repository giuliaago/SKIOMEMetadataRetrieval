# SKIOME Metadata Retrieval

## Download SRAdb

Install SRAdb package if not present

Use SRAdb package to download the SRAdb database containing all the metadata associated with publicly available datasets.

Download the database (>36 Gb):
```{r Download SRAdb}
library(SRAdb)

# Download SRAdb:
sqlfile <- getSRAdbFile()
```

## PART 1 - Search studies of interest with SRAdb and create a dataframe with all the metadata retrived with SRAdb and EDirect

## SRAdb querying and metadata download:

Query of the SRAdb database with a full-text search to recover a list of studies of interest with their associated metadata.

```{r Connect SRAdb and query}
library(SRAdb)

# Connect to the db:
sra_con <- dbConnect(SQLite(),sqlfile)

SRA_dataframe <- getSRA( search_terms ='human skin microbiome OR human skin microbiota OR human skin metagenome', 
                    out_types = c("col_desc", "experiment", "fastq", "metaInfo", "run", "sample", "sra", "sra_ft",
                                  "sra_ft_content", "sra_ft_segdir", "sra_ft_segments", "study", "submission"), sra_con )
dim(SRA_dataframe)
```

Extract from the dataframe created with SRAdb the list of studies
The recovered list can be used to download further information on the run with EDirect tool
```{r Extract the list of studies}
library(readr)
auto_list <- unique(SRA_dataframe$study)
auto_list <- as.data.frame(studies_list)
write_tsv(auto_list, col_names = F, path = "SRA_studies_list.txt")
# Note: it is advisable to store this file in a sub-directory where to run EDirect and store all the metadata files
```

Use the generated file to download run_info file with EDirect.

## EDirect metadata import for the automatic dataset list

```{r Create a list with the names of the EDirect files} 
# convert to vector of characters
auto_list <- as.character(auto_list$studies_list)
# append the _runinfo to each character to create a vector with the names of the file you need to import
auto_list <- paste(auto_list, "_runinfo", sep = "")

```

Import metadata obtained from EDirect tool into a single dataframe
```{r Direct metadata import from automatically searched studies}
# set the directory to where the EDirect runifo files are located
setwd("/mnt/storageN/davide/R_SRAdb/SRA_automatic_studies_metadata")

# create a empty dataframe:
auto_edirect <- data.frame(matrix(ncol = 47, nrow = 0))
# import one file as a temporary variable
tmp <- read.csv(auto_list[1])
# use the imported file to name the dataframe columns
colnames(auto_edirect) <- colnames(tmp)
# remove the temporary file
rm(tmp)

# import all the datasets:
for (k in 1:length(auto_list)){
  tmp <- read.csv(auto_list[k])
  # remove redundant headers
  for (n in nrow(tmp):1){
    if (tmp[n, 1]=="Run"){
      tmp <- tmp[-c(n),]
    }
  }
  # merge datasets
  auto_edirect <- rbind(auto_edirect, tmp)
  rownames(auto_edirect) <- 1:nrow(auto_edirect)
}

setwd("/mnt/storageN/davide/R_SRAdb")
```

## Join SRAdb and E-Direct information in a single dataframe:
```{r Join SRAdb and EDirect dataframes}
library(dplyr)
colnames(SRA_dataframe)[2] <- "Run"
auto_df <- merge(SRA_dataframe, auto_edirect, by = c("Run"), all.x = T, all.y = T)
```

## Add release date column
```{r Create a release date column}
release <- auto_df$ReleaseDate
for (w in 1:length(release)){
  release[w] <- substr(release[w], 1, 4)
}
auto_df <- cbind(auto_df, release)
colnames(auto_df)[117] <- "Year_of_release"
```

## Reorganize the columns:
```{r Columns removal and reorganization}
curated_df <- auto_df %>% select(Run, Year_of_release, ReleaseDate, LoadDate, Submission, spots.y, bases.y, spots_with_mates, avgLength, size_MB, download_path, Experiment, LibraryStrategy, LibraryLayout, InsertSize, Platform, Model, experiment_title, design_description, library_name, experiment_attribute, Sample, BioSample, TaxID, ScientificName, Sex, Body_Site, sample_alias, description, sample_attribute, CenterName, SRAStudy, BioProject, center_project_name, sradb_updated, study_description, study_abstract)

# rename columns
colnames(curated_df) <- c("Run", "Year_of_release", "Release_Date", "Load_Date", "Submission", "Spots", "Bases", "Spots_with_mates", "AvgLength", "Size_MB", "Download_path", "Experiment", "Library_Strategy", "Library_Layout", "Insert_Size", "Platform", "Model", "Experiment_title", "Design_description", "Library_name", "Experiment_attribute", "Sample_ID", "BioSample", "TaxID", "Scientific_Name", "Sex", "Body_Site", "Sample_alias", "Description", "Sample_attribute", "Center_Name", "Study_ID", "BioProject", "Center_project_name", "Sradb_updated", "Study_description", "Study_abstract")
```


## PART 2 - Search studies of interest both manually and with SRAdb and create a single dataframe with all the retrieved studies and the associated metadata retrived with SRAdb, EDirect and manually.

For the automatic search of SRAdb we will use the same dataframe created in part 1

Here we include the manual dataframe with his metadata

## Import the manual dataframe:
```{r Import manual dataframe}
manual_df <- read.csv("Human Skin Datasets.csv")
```

Use SRAdb to specifically search for metadata associated with the manually retrieved list of studies
To this end we use as query the study IDs:
```{r Connect SRAdb and query with the list of studies ID obtained by the manual search}
library(SRAdb)

# Connect to the db:
sra_con <- dbConnect(SQLite(),sqlfile)

SRA_manual <- getSRA( search_terms ='ERP126714 OR SRP144464 OR SRP091029 OR SRP185615 OR DRP003136 OR SRP071579 OR SRP067274 OR SRP062726 OR
                DRP000933OR DRP004449 OR SRP091027 OR ERP017003 OR ERP019566 OR ERP024081 OR ERP016977 OR ERP014521 OR SRP255681 OR 
                SRP090974 OR ERP104625 OR SRP135853 OR ERP109379 OR ERP005296 OR ERP013063 OR SRP262012 OR SRP276976 OR SRP22682 OR 
                SRP285536 OR ERP014863 OR ERP107880 OR SRP109593 OR SRP247381 OR SRP063707 OR SRP187334 OR SRP213896 OR ERP013768 OR 
                SRP252817 OR SRP276001 OR ERP109520 OR SRP034607 OR SRP185750 OR SRP279904 OR SRP109057 OR SRP293771 OR SRP067489 OR 
                ERP114401 OR SRP276479 OR SRP267292 OR SRP167680 OR ERP112507 OR ERP001059 OR SRP230017 OR SRP275714 OR SRP188413 OR 
                SRP050070 OR SRP186628 OR ERP105323 OR SRP057641 OR ERP016529 OR ERP022958 OR SRP113602 OR ERP016629 OR SRP214545 OR 
                ERP018577 OR ERP005182 OR SRP051059 OR SRP133226 OR ERP016280 OR SRP002480', 
                out_types = c("col_desc", "experiment", "fastq", "metaInfo", "run", "sample", "sra", "sra_ft",
                                  "sra_ft_content", "sra_ft_segdir", "sra_ft_segments", "study", "submission"), sra_con )
dim(SRA_manual)
```
Note: Not all the manually found studies are found by SRAdb

## Merge the manual SRAdb dataset and the automatic SRAdb dataset and remove redundant rows (duplicated samples)
```{r Merging the two SRAdb dataframe}
unique(SRA_dataframe$study)
unique(SRA_manual$study)
setdiff(SRA_manual$study, SRA_dataframe$study)
setdiff(SRA_dataframe$study, SRA_manual$study)

# bind the "manual" dataset with the "sradb search" dataset
merged_df <- rbind(SRA_dataframe, SRA_manual)
unique(merged_df$study)

# Remove duplicate rows based on certain columns (variables):
# Remove duplicated rows based on Sepal.Length
merged_df <- merged_df %>% distinct(run, .keep_all = TRUE)
length(unique(merged_df$study))
```


Extract from the dataframes created with SRAdb automatic search and the manual search the list of studies
The recovered list can be used to download further information on the run with EDirect tool
```{r Extract the list of automatic studies}
library(readr)
auto_list <- unique(SRA_dataframe$study)
auto_list <- as.data.frame(auto_list)
manual_list <- unique(manual_df$Study_alias)
manual_list <- as.data.frame(manual_list)
colnames(auto_list)[1] <- "study"
colnames(manual_list)[1] <- "study"
edirect_list <- rbind(auto_list, manual_list)
# remove redundant rows
edirect_list <- edirect_list %>% distinct(study, .keep_all = TRUE)
# create a list with the study_ID of the samples you still haven't downloaded from EDirect
# you will use this list to perform a new EDirect search:
list_for_edirect <- setdiff(manual_list$study, auto_list$study)
list_for_edirect <- as.data.frame(list_for_edirect)
write_tsv(list_for_edirect, col_names = F, path = "SRA_automatic_studies_metadata/list_for_edirect.txt")
# Note: it is advisable to store this file in a sub-directory where to run EDirect and store all the metadata files
```


```{r Create another list with the names of the EDirect files} 
# convert to vector of characters
edirect_list <- as.character(edirect_list$study)
# append the _runinfo to each character to create a vector with the names of the file you need to import
edirect_list <- paste(edirect_list, "_runinfo", sep = "")

```

## EDirect metadata import for the manual dataset list

Import metadata obtained from EDirect tool into a sigle dataframe
```{r EDirect metadata import from manually searched studies}

# set the directory to where the EDirect runifo files are located
setwd("/mnt/storageN/davide/R_SRAdb/SRA_automatic_studies_metadata")

# create a empty dataframe:
edirect_df <- data.frame(matrix(ncol = 47, nrow = 0))
# import one file in a temporary file
tmp <- read.csv(edirect_list[1])
# use the imported file to name the dataframe columns
colnames(edirect_df) <- colnames(tmp)
# remove the temporary file
rm(tmp)

# import all the datasets:
for (k in 1:length(edirect_list)){
  tmp <- read.csv(edirect_list[k])
  # remove redundant headers
  for (n in nrow(tmp):1){
    if (tmp[n, 1]=="Run"){
      tmp <- tmp[-c(n),]
    }
  }
  # merge datasets
  edirect_df <- rbind(edirect_df, tmp)
  rownames(edirect_df) <- 1:nrow(edirect_df)
}

setwd("/mnt/storageN/davide/R_SRAdb")
```

## Merging the dataframes:
```{r Merging the dataframes}
library(dplyr)
# join SRAdb and E-Direct information in a single dataset:
colnames(SRA_dataframe)[2] <- "Run"
colnames(SRA_manual)[2] <- "Run"
merged_df_2 <- rbind(SRA_dataframe, SRA_manual)

merged_df_2 <- merge(merged_df_2, edirect_df, by = c("Run"), all.x = T, all.y = T)
# merge sparse study_id columns:
merged_df_2 <- merged_df_2 %>%
  mutate(study = coalesce(study,SRAStudy))
# include manual retrieved study-specific information:
merged_df_2 <- left_join(merged_df_2, manual_df, by = c("study"))
```

## Add release date column
```{r Create a release date column - part 2}
release <- merged_df_2$ReleaseDate
for (w in 1:length(release)){
  release[w] <- substr(release[w], 1, 4)
}
merged_df_2 <- cbind(merged_df_2, release)
colnames(merged_df_2)[135] <- "Year_of_release"
```

## Reorganize the columns:
```{r Columns removal and reorganization -part 2}
curated_df_2 <- merged_df_2 %>% select(Run, Year_of_release, ReleaseDate, LoadDate, Submission, spots.y, bases.y, spots_with_mates, avgLength, size_MB, download_path, Experiment, LibraryStrategy, LibraryLayout, InsertSize, Platform, Model, experiment_title, design_description, library_name, experiment_attribute, Sample, BioSample, TaxID, ScientificName, Sex, Body_Site, sample_alias, description, sample_attribute, CenterName, SRAStudy, BioProject, center_project_name, sradb_updated, collection_method, vRegion, clustering_used, rSequences.OTUs.ASVs., taxon_database.used_in_the_study., taxon_database_version, coocurrence_study, disease.condition, Location, MGnify.analysis, DOI, Year_Of_Publication, study_description, study_abstract, Notes., Manual_Validation)

# rename columns:
colnames(curated_df_2) <- c("Run", "Year_of_release", "Release_Date", "Load_Date", "Submission", "Spots", "Bases", "Spots_with_mates", "AvgLength", "Size_MB", "Download_path", "Experiment", "Library_Strategy", "Library_Layout", "Insert_Size", "Platform", "Model", "Experiment_title", "Design_description", "Library_name", "Experiment_attribute", "Sample_ID", "BioSample", "TaxID", "Scientific_Name", "Sex", "Body_Site", "Sample_alias", "Description", "Sample_attribute", "Center_Name", "Study_ID", "BioProject", "Center_project_name", "Sradb_updated", "Collection_method", "Region_16S", "Clustering_used", "nSequences_OTUs_ASVs", "Taxon_database_used_in_the_study", "Taxon_database_version", "Coocurrence_study", "Disease_condition", "Location", "MGnify_analysis", "DOI", "Year_Of_Publication", "Study_description", "Study_abstract", "Notes", "Manual_Validation")
```

## Subset only manually retrieved studies:
```{r Subset only manually retrieved studies}
manual_dataframe <- filter(curated_df_2, Manual_Validation == "Manually_Validated")
length(unique(manual_dataframe$Study_ID))
```

## Remove rows and correct errors:
```{r Export the dataframes}
# correct two errors in the library strategy metadata:
# correction applied based on the information on the published article
library(dplyr)
manual_dataframe <- manual_dataframe %>% mutate(Library_Strategy = ifelse(Study_ID == "ERP126714" & Library_Strategy == "RNA-Seq", "AMPLICON", Library_Strategy))

manual_dataframe <- manual_dataframe %>% mutate(Library_Strategy = ifelse(Study_ID == "ERP024081" & Library_Strategy == "WGS", "AMPLICON", Library_Strategy))

# remove WGS rows:
manual_dataframe <- manual_dataframe[manual_dataframe$Library_Strategy != "WGS",]
```

## Export the dataframes:
```{r Export the dataframes to csv}
write.csv(curated_df, "Dataframe_1_Automatic.csv", row.names = FALSE)
write.csv(curated_df_2, "Dataframe_2.csv", row.names = FALSE)
write.csv(manual_dataframe, "Dataframe_3_Manual.csv", row.names = FALSE)
```

# Data Frames exploration:
```{r Exporing the dataframes}
# describe the datasets
dim(curated_df)
length(unique(curated_df$Study_ID))

dim(curated_df_2)
length(unique(curated_df_2$Study_ID))

dim(manual_dataframe)
length(unique(manual_dataframe$Study_ID))

# compare manual vs automatic search
length(setdiff(unique(curated_df$Study_ID), unique(manual_dataframe$Study_ID)))
length(setdiff(unique(manual_dataframe$Study_ID), unique(curated_df$Study_ID)))

# Year
unique(curated_df$Year_of_release)
unique(manual_dataframe$Year_of_release)
unique(manual_dataframe$Year_Of_Publication)

# collection method:
unique(manual_dataframe$Collection_method)
length(which(manual_df$collection_method == "swab"))
length(which(manual_df$collection_method == "biopsy"))
length(which(manual_df$collection_method == "scrubs buffer washes"))
length(which(manual_df$collection_method == "swab and tape-strip"))
length(which(manual_df$collection_method == "swab and biopsy"))
length(which(manual_df$collection_method == "swab, pore strip, glue strip"))

# 16S region:
# studies:
unique(manual_df$vRegion)
length(which(manual_df$vRegion == "Many (V1-V3 & V1-V8)"))
length(which(manual_df$vRegion == "V1-V2"))
length(which(manual_df$vRegion == "V1-V2 and v4"))
length(which(manual_df$vRegion == "V1-V3"))
length(which(manual_df$vRegion == "V1-V4"))
length(which(manual_df$vRegion == "V2-V3"))
length(which(manual_df$vRegion == "V3-V4"))
length(which(manual_df$vRegion == "V3-V5"))
length(which(manual_df$vRegion == "V4"))
length(which(manual_df$vRegion == "V4-V5"))
# samples:
unique(manual_dataframe$Region_16S)
length(which(manual_dataframe$Region_16S == "Many (V1-V3 & V1-V8)"))
length(which(manual_dataframe$Region_16S == "V1-V2"))
length(which(manual_dataframe$Region_16S == "V1-V2 and v4"))
length(which(manual_dataframe$Region_16S == "V1-V3"))
length(which(manual_dataframe$Region_16S == "V1-V4"))
length(which(manual_dataframe$Region_16S == "V2-V3"))
length(which(manual_dataframe$Region_16S == "V3-V4"))
length(which(manual_dataframe$Region_16S == "V3-V5"))
length(which(manual_dataframe$Region_16S == "V4"))
length(which(manual_dataframe$Region_16S == "V4-V5"))

# seq platform:
# manual dataframe
unique(manual_dataframe$Platform)
length(which(manual_dataframe$Platform == "LS454"))
length(which(manual_dataframe$Platform == "ILLUMINA"))

unique(manual_dataframe$Model)
length(which(manual_dataframe$Model == "Illumina MiSeq"))
length(which(manual_dataframe$Model == "Illumina MiniSeq"))
length(which(manual_dataframe$Model == "Illumina HiSeq 2000"))
length(which(manual_dataframe$Model == "Illumina HiSeq 2500"))
length(which(manual_dataframe$Model == "Illumina NovaSeq 6000"))
length(which(manual_dataframe$Model == "454 GS FLX Titanium"))
length(which(manual_dataframe$Model == "454 GS FLX+"))
length(which(manual_dataframe$Model == "454 GS FLX"))
length(which(manual_dataframe$Model == "454 GS"))
length(which(manual_dataframe$Model == "unspecified"))

# manual + automatic dataframe
unique(curated_df_2$Platform)
length(which(curated_df_2$Platform == "LS454"))
length(which(curated_df_2$Platform == "ILLUMINA"))
length(which(curated_df_2$Platform == "ION_TORRENT"))
length(which(curated_df_2$Platform == "CAPILLARY"))

unique(curated_df_2$Model)
length(which(curated_df_2$Model == "Illumina MiSeq"))
length(which(curated_df_2$Model == "Illumina MiniSeq"))
length(which(curated_df_2$Model == "Illumina HiSeq 2000"))
length(which(curated_df_2$Model == "Illumina HiSeq 2500"))
length(which(curated_df_2$Model == "Illumina NovaSeq 6000"))
length(which(curated_df_2$Model == "454 GS FLX Titanium"))
length(which(curated_df_2$Model == "454 GS FLX+"))
length(which(curated_df_2$Model == "454 GS FLX"))
length(which(curated_df_2$Model == "454 GS"))
length(which(curated_df_2$Model == "unspecified"))

# Clustering used:
unique(manual_df$clustering_used)
length(which(manual_df$clustering_used == "OTUs"))
length(which(manual_df$clustering_used == "ASVs"))
length(which(manual_df$clustering_used == "RSVs"))
sum(is.na(manual_df$clustering_used))

# Taxonomic database used:
unique(manual_df$taxon_database.used_in_the_study.)
length(which(manual_df$taxon_database.used_in_the_study. == "Greengenes"))
length(which(manual_df$taxon_database.used_in_the_study. == "SILVA"))
length(which(manual_df$taxon_database.used_in_the_study. == "RDP"))
length(which(manual_df$taxon_database.used_in_the_study. == "NCBI"))
length(which(manual_df$taxon_database.used_in_the_study. == "EzTaxon-e"))
length(which(manual_df$taxon_database.used_in_the_study. == "HOMD and Greengenes"))
sum(is.na(manual_df$taxon_database.used_in_the_study.))

# Co-occurrence analysis:
length(which(manual_df$coocurrence_study == "Yes"))

# Disease/condition:
unique(manual_df$disease.condition)
length(which(manual_df$disease.condition == "Healthy"))

length(which(manual_df$disease.condition == "Atopic Dermatitis"))
length(which(manual_df$disease.condition == "Healthy & Atopic Dermatitis"))
length(which(manual_df$disease.condition == "Atopic Dermatitis and Psoriasis"))
length(which(manual_df$disease.condition == "Psoriasis"))
length(which(manual_df$disease.condition == "Parapsoriasis"))

length(which(manual_df$disease.condition == "Acne"))
length(which(manual_df$disease.condition == "Dandruff"))
length(which(manual_df$disease.condition == "Leprosy"))
length(which(manual_df$disease.condition == "Hidradenitis Suppurativa"))
length(which(manual_df$disease.condition == "Squamous Cell Carcinoma"))
length(which(manual_df$disease.condition == "Filaggrin-deficient human skin"))
length(which(manual_df$disease.condition == "Autoimmune bullous disease"))
length(which(manual_df$disease.condition == "Dystrophic epidermolysis bullosa"))
length(which(manual_df$disease.condition == "Vitiligo"))

# skin injuries:
length(which(manual_df$disease.condition == "Burn Injury"))
length(which(manual_df$disease.condition == "Diabetic Wound"))
length(which(manual_df$disease.condition == "Damaged skin"))
length(which(manual_df$disease.condition == "Wound"))
length(which(manual_df$disease.condition == "Chronic wounds"))
length(which(manual_df$disease.condition == "Open fractures"))

# allegic traits and atopic individuals:
length(which(manual_df$disease.condition == "Allergic traits"))
length(which(manual_df$disease.condition == "Atopic individuals"))

# skin pathogenic infections (bacterial, fungal and parasitic infection)
length(which(manual_df$disease.condition == "Tinea Pedis"))
length(which(manual_df$disease.condition == "Cutaneous leishmaniasis lesions"))
length(which(manual_df$disease.condition == "Buruli ulcer skin lesions"))

# other conditions:
length(which(manual_df$disease.condition == "Low bithweight"))
length(which(manual_df$disease.condition == "Obesity"))


# Location:
unique(manual_df$Location)
length(which(manual_df$Location == "China"))
length(which(manual_df$Location == "USA"))

# Europe:
length(which(manual_df$Location == "Ireland"))
length(which(manual_df$Location == "Finland"))
length(which(manual_df$Location == "Denmark"))
length(which(manual_df$Location == "Belgium"))
length(which(manual_df$Location == "Netherlands"))
length(which(manual_df$Location == "Germany"))
length(which(manual_df$Location == "France"))
length(which(manual_df$Location == "Austria"))
length(which(manual_df$Location == "UK"))
length(which(manual_df$Location == "Europe"))

length(which(manual_df$Location == "Canada and UK"))
length(which(manual_df$Location == "Canada"))
length(which(manual_df$Location == "Japan"))
length(which(manual_df$Location == "Australia"))
length(which(manual_df$Location == "Korea"))
length(which(manual_df$Location == "India"))
length(which(manual_df$Location == "Brazil"))
length(which(manual_df$Location == "Israel"))
length(which(manual_df$Location == "Egypt"))
length(which(manual_df$Location == "Benin (Africa)"))

# Taxon ID and scientific name
unique(manual_dataframe$TaxID)
unique(manual_dataframe$Scientific_Name)

length(unique(curated_df_2$TaxID))
unique(curated_df_2$Scientific_Name)
```

Contacts
-----------------------------------------------------
* Name: __Giulia Agostinetto__
* e-mail: __g.agostinetto@campus.unimib.it__
* Name: __Davide Bozzi__
* e-mail: __d.bozzi1@campus.unimib.it__

Please refer to both contacts for further information or issues related to the framework.
