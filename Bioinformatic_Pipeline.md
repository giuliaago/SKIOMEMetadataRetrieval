# SKIOME Datasets and Metadata Retrieval - Pipeline

## Download SRAdb

Install SRAdb package if not present

Use SRAdb package to download the SRAdb database containing all the metadata associated with publicly available datasets.

## Download the database (>36 Gb):
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
write_tsv(auto_list, col_names = F, path = "EDirect_metadata/SRA_studies_list.txt")
# Note: it is advisable to store this file in a sub-directory where to run EDirect and store all the metadata files (e.g. "EDirect_metadata" directory)
```

Use the generated file to download run_info file with EDirect.

## Download metadata with EDirect (from a Unix terminal window)
See https://www.ncbi.nlm.nih.gov/books/NBK179288/ for EDirect installation instructions.

```
sample_list=$(cat  SRA_studies_list.txt)
for a in $sample_list
   do
   esearch -db sra -query "$a" | efetch -format runinfo > "$a"_runinfo # download metadata for the runs
done
```

## Import metadata obtained from EDirect tool into a single dataframe
```{r Direct metadata import from automatically searched studies}
# set the directory to where the EDirect runifo files are located
setwd("FULL_PATH/EDirect_metadata")

# Create a list with the names of the EDirect files:
# convert to vector of characters
auto_list <- as.character(auto_list$studies_list)
# append the _runinfo to each character to create a vector with the names of the file you need to import
auto_list <- paste(auto_list, "_runinfo", sep = "")

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

setwd("FULL_PATH_TO_WORKING_DIRECTORY")
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

## Use SRAdb to specifically search for metadata associated with the manually retrieved list of studies
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
write_tsv(list_for_edirect, col_names = F, path = "EDirect_metadata/list_for_edirect.txt")
# Note: it is advisable to store this file in a sub-directory where to run EDirect and store all the metadata files
```

## Download metadata with EDirect (from a Unix terminal window)
```
sample_list=$(cat  list_for_edirect.txt)
for a in $sample_list
   do
   esearch -db sra -query "$a" | efetch -format runinfo > "$a"_runinfo # download metadata for the runs
done
```


## EDirect metadata import for the manual dataset list

Import metadata obtained from EDirect tool into a sigle dataframe
```{r EDirect metadata import from manually searched studies}
# set the directory to where the EDirect runifo files are located
setwd("/FULL_PATH/EDirect_metadata")

# Create another list with the names of the EDirect files:
# convert to vector of characters
edirect_list <- as.character(edirect_list$study)
# append the _runinfo to each character to create a vector with the names of the file you need to import
edirect_list <- paste(edirect_list, "_runinfo", sep = "")

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

setwd("FULL_PATH_TO_WORKING_DIRECTORY")
```

## Merging the dataframes:
```{r Merging the dataframes}
library(dplyr)

# join SRAdb and E-Direct information in a single dataset:
colnames(merged_df)[2] <- "Run"
merged_df_2 <- merge(merged_df, edirect_df, by = c("Run"), all.x = T, all.y = T)
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
colnames(merged_df_2)[146] <- "Year_of_release"
```

## Reorganize the columns:
```{r Columns removal and reorganization -part 2}
curated_df_2 <- merged_df_2 %>% select(Run, Year_of_release, ReleaseDate, LoadDate, Submission, spots.y, bases.y, spots_with_mates, avgLength, size_MB, download_path, Experiment, LibraryStrategy, LibraryLayout, InsertSize, Platform, Model, experiment_title, design_description, library_name, experiment_attribute, Sample, BioSample, TaxID, ScientificName, Sex, Body_Site, sample_alias, description, sample_attribute, CenterName, SRAStudy, BioProject, center_project_name, sradb_updated, collection_method, vRegion, clustering_used, rSequences.OTUs.ASVs., taxon_database.used_in_the_study., taxon_database_version, coocurrence_study, disease.condition, Location, MGnify.analysis, DOI, Year_Of_Publication, Journal, WOS_Research_Areas, WOS_Categories, Scopus_Research_Subject_1, Scopus_Research_Subject_2, Scopus_Research_Subject_3, Scopus_Research_Subject_4, Scopus_Research_Subject_5, Scopus_Research_Subject_6, Medicine_journal, study_description, study_abstract, Notes., Manual_Validation)

# rename columns:
colnames(curated_df_2) <- c("Run", "Year_of_release", "Release_Date", "Load_Date", "Submission", "Spots", "Bases", "Spots_with_mates", "AvgLength", "Size_MB", "Download_path", "Experiment", "Library_Strategy", "Library_Layout", "Insert_Size", "Platform", "Model", "Experiment_title", "Design_description", "Library_name", "Experiment_attribute", "Sample_ID", "BioSample", "TaxID", "Scientific_Name", "Sex", "Body_Site", "Sample_alias", "Description", "Sample_attribute", "Center_Name", "Study_ID", "BioProject", "Center_project_name", "Sradb_updated", "Collection_method", "Region_16S", "Clustering_used", "nSequences_OTUs_ASVs", "Taxon_database_used_in_the_study", "Taxon_database_version", "Coocurrence_study", "Disease_condition", "Location", "MGnify_analysis", "DOI", "Year_Of_Publication", "Journal", "WOS_Research_Areas", "WOS_Categories", "Scopus_Research_Subject_1", "Scopus_Research_Subject_2", "Scopus_Research_Subject_3", "Scopus_Research_Subject_4", "Scopus_Research_Subject_5", "Scopus_Research_Subject_6", "Medicine_journal", "Study_description", "Study_abstract", "Notes", "Manual_Validation")
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

## Collapse the Scopus research subjects into a single column:
```{r}
manual_dataframe$Scopus_Research_Subject <- paste(manual_dataframe$Scopus_Research_Subject_1, manual_dataframe$Scopus_Research_Subject_2, manual_dataframe$Scopus_Research_Subject_3, manual_dataframe$Scopus_Research_Subject_4, manual_dataframe$Scopus_Research_Subject_5, manual_dataframe$Scopus_Research_Subject_6, sep = "; ")

curated_df_2$Scopus_Research_Subject <- paste(curated_df_2$Scopus_Research_Subject_1, curated_df_2$Scopus_Research_Subject_2, curated_df_2$Scopus_Research_Subject_3, curated_df_2$Scopus_Research_Subject_4, curated_df_2$Scopus_Research_Subject_5, curated_df_2$Scopus_Research_Subject_6, sep = "; ")
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

length(which(curated_df$Year_of_release == "2015"))
length(which(curated_df_2$Year_of_release == "2017"))
length(which(manual_dataframe$Year_of_release == "2017"))
length(which(manual_dataframe$Year_of_release == "2013"))
length(which(manual_df$Year_Of_Publication == "2019"))

length(which(curated_df$Year_of_release == "2008"))
length(which(curated_df$Year_of_release == "2009"))
length(which(curated_df$Year_of_release == "2010"))
length(which(curated_df$Year_of_release == "2011"))

# collection method:
unique(manual_dataframe$Collection_method)
length(which(manual_df$collection_method == "swab"))
length(which(manual_df$collection_method == "biopsy"))
length(which(manual_df$collection_method == "scrubs buffer washes"))
length(which(manual_df$collection_method == "swab and tape-strip"))
length(which(manual_df$collection_method == "swab and biopsy"))
length(which(manual_df$collection_method == "swab, pore strip, glue strip"))

length(which(manual_dataframe$Collection_method == "swab"))
length(which(manual_dataframe$Collection_method == "biopsy"))
length(which(manual_dataframe$Collection_method == "scrubs buffer washes"))
length(which(manual_dataframe$Collection_method == "swab and tape-strip"))
length(which(manual_dataframe$Collection_method == "swab and biopsy"))
length(which(manual_dataframe$Collection_method == "swab, pore strip, glue strip"))

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
length(which(manual_df$disease.condition == "Low birthweight"))
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

# Other:
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

# Spots: dataframe1
sum(is.na(curated_df$Spots))
min(as.numeric(curated_df$Spots))
length(which(curated_df$Spots == 0))
max(as.numeric(curated_df$Spots))
median(as.numeric(curated_df$Spots))
tmp <- as.numeric(curated_df$Spots)
mean(tmp[tmp != 0])
median(tmp[tmp != 0])

# spots: dataframe2
median(as.numeric(curated_df_2$Spots))
tmp <- as.numeric(curated_df_2$Spots)
mean(tmp[tmp != 0])
median(tmp[tmp != 0])

# spots: dataframe3
median(as.numeric(manual_dataframe$Spots))
tmp <- as.numeric(manual_dataframe$Spots)
mean(tmp[tmp != 0])
median(tmp[tmp != 0])

# Bases: dataframe1
sum(is.na(curated_df$Bases))
min(as.numeric(curated_df$Bases))
length(which(curated_df$Bases == 0))
max(as.numeric(curated_df$Bases))
mean(as.numeric(curated_df$Bases))
sd(as.numeric(curated_df$Bases))
median(as.numeric(curated_df$Bases))

tmp <- as.numeric(curated_df$Bases)
mean(tmp[tmp != 0])
median(tmp[tmp != 0])

# Bases: dataframe2
median(as.numeric(curated_df_2$Bases))
tmp <- as.numeric(curated_df_2$Bases)
mean(tmp[tmp != 0])
median(tmp[tmp != 0])

# Bases: dataframe3
median(as.numeric(manual_dataframe$Bases))
tmp <- as.numeric(manual_dataframe$Bases)
mean(tmp[tmp != 0])
median(tmp[tmp != 0])

# Check AvgLength: dataframe1
length(which(curated_df$AvgLength == "0"))
tmp <- as.numeric(curated_df$AvgLength)
mean(tmp[tmp != 0])
median(tmp[tmp != 0])

# Check AvgLength: dataframe2
length(which(curated_df_2$AvgLength == "0"))
tmp <- as.numeric(curated_df_2$AvgLength)
mean(tmp[tmp != 0])
median(tmp[tmp != 0])

# Check AvgLength: dataframe3
length(which(manual_dataframe$AvgLength == "0"))
mean(as.numeric(manual_dataframe$AvgLength))
sd(as.numeric(manual_dataframe$AvgLength))
median(as.numeric(manual_dataframe$AvgLength))

# insert size: dataframe1:
length(which(curated_df$Insert_Size == "0"))
mean(as.numeric(curated_df$Insert_Size))
sd(as.numeric(curated_df$Insert_Size))
median(as.numeric(curated_df$Insert_Size))
tmp <- as.numeric(curated_df$Insert_Size)
mean(tmp[tmp != 0])
median(tmp[tmp != 0])
# insert size: dataframe2:
length(which(curated_df_2$Insert_Size == "0"))
mean(as.numeric(curated_df_2$Insert_Size))
sd(as.numeric(curated_df_2$Insert_Size))
median(as.numeric(curated_df_2$Insert_Size))
tmp <- as.numeric(curated_df_2$Insert_Size)
mean(tmp[tmp != 0])
median(tmp[tmp != 0])
# insert size: dataframe3:
length(which(manual_dataframe$Insert_Size == "0"))
mean(as.numeric(manual_dataframe$Insert_Size))
sd(as.numeric(manual_dataframe$Insert_Size))
median(as.numeric(manual_dataframe$Insert_Size))
tmp <- as.numeric(manual_dataframe$Insert_Size)
mean(tmp[tmp != 0])
median(tmp[tmp != 0])

# Sex:
unique(curated_df_2$Sex)

#dataframe1
length(which(curated_df$Sex == "male"))
length(which(curated_df$Sex == "female"))

#dataframe2
length(which(curated_df_2$Sex == "male"))
length(which(curated_df_2$Sex == "female"))

dataframe3
length(which(manual_dataframe$Sex == "male"))
length(which(manual_dataframe$Sex == "female"))

# Body site:
unique(curated_df_2$Body_Site)

# dataframe1
length(which(curated_df$Body_Site == ""))
length(which(curated_df$Body_Site == "Unknown"))
length(which(curated_df$Body_Site == "not collected"))
length(which(curated_df$Body_Site == "BODY_SITE"))
length(which(curated_df$Body_Site == "Unspecified"))
length(which(curated_df$Body_Site == "Not applicable"))
length(which(curated_df$Body_Site == "not applicable"))

# dataframe2
sum(is.na(curated_df_2$Body_Site))
length(which(curated_df_2$Body_Site == ""))
length(which(curated_df_2$Body_Site == "Unknown"))
length(which(curated_df_2$Body_Site == "not collected"))
length(which(curated_df_2$Body_Site == "BODY_SITE"))
length(which(curated_df_2$Body_Site == "Unspecified"))
length(which(curated_df_2$Body_Site == "Not applicable"))
length(which(curated_df_2$Body_Site == "not applicable"))
length(which(curated_df_2$Body_Site == "control"))
length(which(curated_df_2$Body_Site == "ctrl"))
length(which(curated_df_2$Body_Site == "epidermis"))

# dataframe3
length(which(manual_dataframe$Body_Site == ""))
length(which(manual_dataframe$Body_Site == "Unknown"))
length(which(manual_dataframe$Body_Site == "not collected"))
length(which(manual_dataframe$Body_Site == "BODY_SITE"))
length(which(manual_dataframe$Body_Site == "Unspecified"))
length(which(manual_dataframe$Body_Site == "Not applicable"))
length(which(manual_dataframe$Body_Site == "not applicable"))

# Journal
unique(manual_df$Journal)
unique(manual_df$WOS_Research_Areas)
unique(manual_df$WOS_Categories)
length(which(manual_df$Medicine_journal == "Yes"))
```
## Code used to generate plots
```{r Plots generation}
# Plots:
library(ggplot2)

# Year of release:

# data frame 2:
vec <- curated_df_2$Year_of_release
vec <- as.numeric(vec)
vec2 <- unique(vec)
count <- c()

for (j in 1:length(vec2)){
  count <- append(count, length(which(vec==vec2[j])))
}

#create a dataframe
plotdf <- data.frame(vec2, count)

ggplot(data=plotdf, aes(x=vec2, y=count)) +
  geom_line(stat="identity") +
  scale_x_continuous("vec2", labels = as.character(vec2), breaks = vec2) +
  xlab("Year of release") + ylab("Number of samples") +
  ggtitle("Year of release: samples (dataframe2)") +
  theme_bw()

ggplot(data=plotdf, aes(x=vec2, y=count)) +
  geom_bar(stat="identity") +
  scale_x_continuous("vec2", labels = as.character(vec2), breaks = vec2) +
  xlab("Year of release") + ylab("Number of samples") +
  ggtitle("Year of release: samples (dataframe2")

# data frame 3:

vec <- manual_dataframe$Year_of_release
vec <- as.numeric(vec)
vec2 <- unique(vec)
count <- c()

for (j in 1:length(vec2)){
  count <- append(count, length(which(vec==vec2[j])))
}

#create a dataframe
plotdf <- data.frame(vec2, count)

ggplot(data=plotdf, aes(x=vec2, y=count)) +
  geom_line(stat="identity") +
  scale_x_continuous("vec2", labels = as.character(vec2), breaks = vec2) +
  xlab("Year of release") + ylab("Number of samples") +
  ggtitle("Year of release: samples (dataframe3)") +
  theme_bw()

ggplot(data=plotdf, aes(x=vec2, y=count)) +
  geom_bar(stat="identity") +
  scale_x_continuous("vec2", labels = as.character(vec2), breaks = vec2) +
  xlab("Year of release") + ylab("Number of samples") +
  ggtitle("Year of release: samples (dataframe3)")


# data frame 3: (studies)

vec <- manual_df$Year_Of_Publication
add <- unlist(strsplit(vec[68], "[\n]"))
vec <- head(vec, -1)
vec <- append(vec, add)
rm(add)
vec <- as.numeric(vec)
vec2 <- unique(vec)
count <- c()

for (j in 1:length(vec2)){
  count <- append(count, length(which(vec==vec2[j])))
}

#create a dataframe
plotdf <- data.frame(vec2, count)

ggplot(data=plotdf, aes(x=vec2, y=count)) +
  geom_line(stat="identity") +
  scale_x_continuous("vec2", labels = as.character(vec2), breaks = vec2) +
  xlab("Year of release") + ylab("Number of studies") +  
  ggtitle("Year of release: studies (dataframe3)") +
  theme_bw()

ggplot(data=plotdf, aes(x=vec2, y=count)) +
  geom_bar(stat="identity") +
  scale_x_continuous("vec2", labels = as.character(vec2), breaks = vec2) +
  xlab("Year of release") + ylab("Number of studies") +
  ggtitle("Year of release: studies (dataframe3)")

# data frame 1:

vec <- curated_df$Year_of_release
vec <- as.numeric(vec)
vec2 <- unique(vec)
count <- c()

for (j in 1:length(vec2)){
  count <- append(count, length(which(vec==vec2[j])))
}

#create a dataframe
plotdf <- data.frame(vec2, count)

ggplot(data=plotdf, aes(x=vec2, y=count)) +
  geom_line(stat="identity") +
  scale_x_continuous("vec2", labels = as.character(vec2), breaks = vec2) +
  xlab("Year of release") + ylab("Number of samples") + 
  ggtitle("Year of release: samples (dataframe1)") +
  theme_bw()

ggplot(data=plotdf, aes(x=vec2, y=count)) +
  geom_bar(stat="identity") +
  scale_x_continuous("vec2", labels = as.character(vec2), breaks = vec2) +
  xlab("Year of release") + ylab("Number of samples") +
  ggtitle("Year of release: samples (dataframe1)")


# Seq Platform info - daraframe2:
vec <- curated_df_2$Platform
vec2 <- unique(vec)
count <- c()

for (j in 1:length(vec2)){
  count <- append(count, length(which(vec==vec2[j])))
}

#create a dataframe
plotdf <- data.frame(vec2, count)

ggplot(data=plotdf, aes(x=vec2, y=count)) +
  geom_bar(stat="identity") +
  labs(y= "Samples", x = "Sequencing Platform") +
  ggtitle("Seq Platform - daraframe2: (samples)")


# Seq Platform info - daraframe3:
vec <- manual_dataframe$Platform
vec2 <- unique(vec)
count <- c()

for (j in 1:length(vec2)){
  count <- append(count, length(which(vec==vec2[j])))
}

#create a dataframe
plotdf <- data.frame(vec2, count)

ggplot(data=plotdf, aes(x=vec2, y=count)) +
  geom_bar(stat="identity") +
  labs(y= "Samples", x = "Sequencing Platform") +
  ggtitle("Seq Platform - daraframe3: (samples)")


# Seq Platform Model info - daraframe2:
vec <- curated_df_2$Model
vec2 <- unique(vec)
count <- c()

for (j in 1:length(vec2)){
  count <- append(count, length(which(vec==vec2[j])))
}

#create a dataframe
plotdf <- data.frame(vec2, count)

ggplot(data=plotdf, aes(x=vec2, y=count)) +
  geom_bar(stat="identity") +
  labs(y= "Samples", x = "Sequencing Platform Model") +
  ggtitle("Seq Platform Model info - daraframe2: (samples)")


# Seq Platform Model info - daraframe3:
vec <- manual_dataframe$Model
vec2 <- unique(vec)
count <- c()

for (j in 1:length(vec2)){
  count <- append(count, length(which(vec==vec2[j])))
}

#create a dataframe
plotdf <- data.frame(vec2, count)

ggplot(data=plotdf, aes(x=vec2, y=count)) +
  geom_bar(stat="identity") +
  labs(y= "Samples", x = "Sequencing Platform Model")  +
  ggtitle("Seq Platform Model info - daraframe3: (samples)")


# Seq Platform Model info - daraframe3: (studies)
vec <- manual_df$Seq_Platform
vec2 <- unique(vec)
count <- c()

for (j in 1:length(vec2)){
  count <- append(count, length(which(vec==vec2[j])))
}

#create a dataframe
plotdf <- data.frame(vec2, count)

ggplot(data=plotdf, aes(x=vec2, y=count)) +
  geom_bar(stat="identity") +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  labs(y= "Studies", x = "Sequencing Platform Model") +
  ggtitle("Seq Platform Model info - daraframe3: (studies)")


# Info on the 16S region: (studies) - dataset 3
vec <- manual_df$vRegion
vec2 <- unique(vec)
count <- c()

for (j in 1:length(vec2)){
  count <- append(count, length(which(vec==vec2[j])))
}

#create a dataframe
plotdf <- data.frame(vec2, count)

ggplot(data=plotdf, aes(x=vec2, y=count)) +
  geom_bar(stat="identity") +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  labs(y= "Studies", x = "16S rRNA hypervariable region sequenced") +
  ggtitle("16S rRNA region - daraframe3: (studies)")

# Info on the 16S region: (samples) - dataset 3
vec <- manual_dataframe$Region_16S
vec2 <- unique(vec)
count <- c()

for (j in 1:length(vec2)){
  count <- append(count, length(which(vec==vec2[j])))
}

#create a dataframe
plotdf <- data.frame(vec2, count)

ggplot(data=plotdf, aes(x=vec2, y=count)) +
  geom_bar(stat="identity") +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  labs(y= "Samples", x = "16S rRNA hypervariable region sequenced") +
  ggtitle("16S rRNA region - daraframe3: (samples)")


# Collection method (for the samples): -dataframe3
vec <- manual_dataframe$Collection_method
vec2 <- unique(vec)
count <- c()

for (j in 1:length(vec2)){
  count <- append(count, length(which(vec==vec2[j])))
}

#create a dataframe
plotdf <- data.frame(vec2, count)

ggplot(data=plotdf, aes(x=vec2, y=count)) +
  geom_bar(stat="identity") +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  labs(y= "Samples", x = "Collection method") +
  ggtitle("Collection method - daraframe3: (samples)")

# Collection method (for the studies):-dataframe3
vec <- manual_df$collection_method
vec2 <- unique(vec)
count <- c()

for (j in 1:length(vec2)){
  count <- append(count, length(which(vec==vec2[j])))
}

#create a dataframe
plotdf <- data.frame(vec2, count)

ggplot(data=plotdf, aes(x=vec2, y=count)) +
  geom_bar(stat="identity") +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  labs(y= "Studies", x = "Collection method") +
  ggtitle("Collection method - daraframe3: (studies)")

# Clustering method used in the studies:
vec <- manual_df$clustering_used
vec2 <- unique(vec)
count <- c()

for (j in 1:length(vec2)){
  count <- append(count, length(which(vec==vec2[j])))
}

#create a dataframe
plotdf <- data.frame(vec2, count)

ggplot(data=plotdf, aes(x=vec2, y=count)) +
  geom_bar(stat="identity") +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  labs(y= "Studies", x = "Clustering method") +
  ggtitle("Clustering method - daraframe3: (studies)")

# Disease and conditions: (studies) - dataframe3
vec <- manual_df$disease.condition
vec2 <- unique(vec)
count <- c()

for (j in 1:length(vec2)){
  count <- append(count, length(which(vec==vec2[j])))
}

#create a dataframe
plotdf <- data.frame(vec2, count)

ggplot(data=plotdf, aes(x=vec2, y=count)) +
  geom_bar(stat="identity") +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  labs(y= "Studies", x = "Disease/condition") +
  ggtitle("Diseases and conditions - daraframe3: (studies)")


# Taxonomic database:
vec <- manual_df$taxon_database.used_in_the_study.
vec2 <- unique(vec)
count <- c()

for (j in 1:length(vec2)){
  count <- append(count, length(which(vec==vec2[j])))
}

#create a dataframe
plotdf <- data.frame(vec2, count)

ggplot(data=plotdf, aes(x=vec2, y=count)) +
  geom_bar(stat="identity") +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  labs(y= "Studies", x = "Taxonomic database used in the study") +
  ggtitle("Taxonomic database used in the study - daraframe3")


# Location of the study:
vec <- manual_df$Location
vec2 <- unique(vec)
count <- c()

for (j in 1:length(vec2)){
  count <- append(count, length(which(vec==vec2[j])))
}

#create a dataframe
plotdf <- data.frame(vec2, count)

ggplot(data=plotdf, aes(x=vec2, y=count)) +
  geom_bar(stat="identity") +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  labs(y= "Studies", x = "Location of the study") +
  ggtitle("Location of the study - daraframe3")

names_location <- plotdf$vec2
number <- c()
for(k in 1:length(names_location)){
  number <- append(number, length(which(manual_df$Location == names_location[k])))
}
plotdf <- data.frame(names_location, number)

library(viridis)
ggplot(plotdf, aes(x = names_location, y = names_location, size = number, fill = number)) +
  geom_point(shape = 21) +
  scale_fill_viridis_c(option = "plasma")


# Sex:

vec2 <- c("female_df1", "male_df1", "NA_df1", "female_df2", "male_df2", "NA_df2", "female_df3", "male_df3", "NA_df3")
count <- c()

count <- append(count, length(which(curated_df$Sex == "male")))
n <- length(which(curated_df$Sex == "male"))
count <- append(count, length(which(curated_df$Sex == "female")))
n <- n + length(which(curated_df$Sex == "female"))
l <- nrow(curated_df)
diff <- l-n
count <- append(count, diff)

count <- append(count, length(which(curated_df_2$Sex == "male")))
n <- length(which(curated_df_2$Sex == "male"))
count <- append(count, length(which(curated_df_2$Sex == "female")))
n <- n + length(which(curated_df_2$Sex == "female"))
l <- nrow(curated_df_2)
diff <- l-n
count <- append(count, diff)

count <- append(count, length(which(manual_dataframe$Sex == "male")))
n <- length(which(manual_dataframe$Sex == "male"))
count <- append(count, length(which(manual_dataframe$Sex == "female")))
n <- n + length(which(manual_dataframe$Sex == "female"))
l <- nrow(manual_dataframe)
diff <- l-n
count <- append(count, diff)

#create a dataframe
plotdf <- data.frame(vec2, count)

ggplot(data=plotdf, aes(x=vec2, y=count)) +
  geom_bar(stat="identity") +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  labs(y= "Samples", x = "Sex distribution among datasets") +
  ggtitle("Sex distribution among datasets")


# Taxon ID:

vec2 <- c("df1", "df2", "df3")
count <- c()
count <- append(count, length(unique(curated_df$TaxID)))
count <- append(count, length(unique(curated_df_2$TaxID)))
count <- append(count, length(unique(manual_dataframe$TaxID)))

#create a dataframe
plotdf <- data.frame(vec2, count)

ggplot(data=plotdf, aes(x=vec2, y=count)) +
  geom_bar(stat="identity") +
  xlab("Data Frames") + ylab("Number of Taxon ID") +
  ggtitle("TaxID") +
  theme_bw()


# taxID

vec <- manual_dataframe$Scientific_Name
vec2 <- unique(vec)
count <- c()

for (j in 1:length(vec2)){
  count <- append(count, length(which(vec==vec2[j])))
}

# use log 
count <- log(count)

#create a dataframe
plotdf <- data.frame(vec2, count)


ggplot(plotdf, aes(x="", y=count, fill=vec2)) +
  geom_bar(stat="identity", width=1) + 
  scale_fill_manual(values = c("dodgerblue2", 
                               "#E31A1C", 
                               "green4",
                               "#6A3D9A",
                               "#FF7F00",
                               "gold1",
                               "skyblue2", 
                               "#FB9A99",
                               "palegreen2",
                               "#CAB2D6",
                               "#FDBF6F",
                               "gray70", 
                               "khaki2",
                               "maroon", 
                               "orchid1", 
                               "deeppink1", 
                               "blue1", 
                               "steelblue4",
                               "darkturquoise", 
                               "green1", 
                               "yellow4", 
                               "yellow3",
                               "darkorange4", 
                               "brown")) +
 coord_polar("y", start=0)


# Spots Plot

vec2 <- c("df1", "df2", "df3")
count <- c()
# Spots: dataframe1
tmp <- as.numeric(curated_df$Spots)
count <- append(count, median(tmp[tmp != 0])) 
# spots: dataframe2
tmp <- as.numeric(curated_df_2$Spots)
count <- append(count, median(tmp[tmp != 0])) 
# spots: dataframe3
tmp <- as.numeric(manual_dataframe$Spots)
count <- append(count, median(tmp[tmp != 0])) 

#create a dataframe
plotdf <- data.frame(vec2, count)

ggplot(data=plotdf, aes(x=vec2, y=count)) +
  geom_bar(stat="identity") +
  xlab("Data Frames") + ylab("Spots") +
  ggtitle("Spots") +
  theme_bw()


# Bases Plot

vec2 <- c("df1", "df2", "df3")
count <- c()
# Spots: dataframe1
tmp <- as.numeric(curated_df$Bases)
count <- append(count, median(tmp[tmp != 0])) 
# spots: dataframe2
tmp <- as.numeric(curated_df_2$Bases)
count <- append(count, median(tmp[tmp != 0])) 
# spots: dataframe3
tmp <- as.numeric(manual_dataframe$Bases)
count <- append(count, median(tmp[tmp != 0])) 

#create a dataframe
plotdf <- data.frame(vec2, count)

ggplot(data=plotdf, aes(x=vec2, y=count)) +
  geom_bar(stat="identity") +
  xlab("Data Frames") + ylab("Bases") +
  ggtitle("Bases") +
  theme_bw()


# Avg Length Plot

vec2 <- c("df1", "df2", "df3")
count <- c()
# Spots: dataframe1
tmp <- as.numeric(curated_df$AvgLength)
count <- append(count, median(tmp[tmp != 0])) 
# spots: dataframe2
tmp <- as.numeric(curated_df_2$AvgLength)
count <- append(count, median(tmp[tmp != 0])) 
# spots: dataframe3
tmp <- as.numeric(manual_dataframe$AvgLength)
count <- append(count, median(tmp[tmp != 0])) 

#create a dataframe
plotdf <- data.frame(vec2, count)

ggplot(data=plotdf, aes(x=vec2, y=count)) +
  geom_bar(stat="identity") +
  xlab("Data Frames") + ylab("Average Length") +
  ggtitle("Average Length") +
  theme_bw()

# Avg Length - df1
vec <- curated_df$AvgLength
vec <- as.numeric(vec)
vec2 <- unique(vec)
count <- c()

for (j in 1:length(vec2)){
  count <- append(count, length(which(vec==vec2[j])))
}

#create a dataframe
plotdf <- data.frame(vec2, count)

ggplot(data=plotdf, aes(x=vec2, y=count)) +
  geom_bar(stat="identity") +
  scale_x_binned() +
  xlab("Average Length") + ylab("Number of samples") +
  ggtitle("Average Length (dataframe1)") +
  theme(axis.text.x = element_text(angle = 45, hjust=1))

# Avg Length - df2
vec <- curated_df_2$AvgLength
vec <- as.numeric(vec)
vec2 <- unique(vec)
count <- c()

for (j in 1:length(vec2)){
  count <- append(count, length(which(vec==vec2[j])))
}

#create a dataframe
plotdf <- data.frame(vec2, count)

ggplot(data=plotdf, aes(x=vec2, y=count)) +
  geom_bar(stat="identity") +
  scale_x_binned() +
  xlab("Average Length") + ylab("Number of samples") +
  ggtitle("Average Length (dataframe2)") +
  theme(axis.text.x = element_text(angle = 45, hjust=1))

# Avg Length - df3
vec <- manual_dataframe$AvgLength
vec <- as.numeric(vec)
vec2 <- unique(vec)
count <- c()

for (j in 1:length(vec2)){
  count <- append(count, length(which(vec==vec2[j])))
}

#create a dataframe
plotdf <- data.frame(vec2, count)

ggplot(data=plotdf, aes(x=vec2, y=count)) +
  geom_bar(stat="identity") +
  scale_x_binned() +
  xlab("Average Length") + ylab("Number of samples") +
  ggtitle("Average Length (dataframe3)") +
  theme(axis.text.x = element_text(angle = 45, hjust=1))


# Insert size Plot

vec2 <- c("df1", "df2", "df3")
count <- c()
# Spots: dataframe1
tmp <- as.numeric(curated_df$Insert_Size)
count <- append(count, median(tmp[tmp != 0])) 
# spots: dataframe2
tmp <- as.numeric(curated_df_2$Insert_Size)
count <- append(count, median(tmp[tmp != 0])) 
# spots: dataframe3
tmp <- as.numeric(manual_dataframe$Insert_Size)
count <- append(count, median(tmp[tmp != 0])) 

#create a dataframe
plotdf <- data.frame(vec2, count)

ggplot(data=plotdf, aes(x=vec2, y=count)) +
  geom_bar(stat="identity") +
  xlab("Data Frames") + ylab("Insert Size") +
  ggtitle("Insert Size") +
  theme_bw()

# Insert Size - df1
vec <- curated_df$Insert_Size
vec <- as.numeric(vec)
vec2 <- unique(vec)
vec2 <- vec2[2:length(vec2)] # remove 0
count <- c()

for (j in 1:length(vec2)){
  count <- append(count, length(which(vec==vec2[j])))
}

#create a dataframe
plotdf <- data.frame(vec2, count)
plotdf <- plotdf %>% filter(vec2 < 1000)

ggplot(data=plotdf, aes(x=vec2, y=count)) +
  geom_bar(stat="identity") +
  scale_x_binned(n.breaks = 10) +
  xlab("Insert Size") + ylab("Number of samples") +
  ggtitle("Insert Size (dataframe1)") +
  theme(axis.text.x = element_text(angle = 45, hjust=1))

# Insert Size - df2
vec <- curated_df_2$Insert_Size
vec <- as.numeric(vec)
vec2 <- unique(vec)
vec2 <- vec2[2:length(vec2)] # remove 0
count <- c()

for (j in 1:length(vec2)){
  count <- append(count, length(which(vec==vec2[j])))
}

#create a dataframe
plotdf <- data.frame(vec2, count)
plotdf <- plotdf %>% filter(vec2 < 1000)

ggplot(data=plotdf, aes(x=vec2, y=count)) +
  geom_bar(stat="identity") +
  scale_x_binned(n.breaks = 10) +
  xlab("Insert Size") + ylab("Number of samples") +
  ggtitle("Insert Size (dataframe2)") +
  theme(axis.text.x = element_text(angle = 45, hjust=1))

# Insert Size - df3
vec <- manual_dataframe$Insert_Size
vec <- as.numeric(vec)
vec2 <- unique(vec)
vec2 <- vec2[2:length(vec2)] # remove 0
count <- c()

for (j in 1:length(vec2)){
  count <- append(count, length(which(vec==vec2[j])))
}

#create a dataframe
plotdf <- data.frame(vec2, count)
plotdf <- plotdf %>% filter(vec2 < 1000)

ggplot(data=plotdf, aes(x=vec2, y=count)) +
  geom_bar(stat="identity") +
  scale_x_binned(n.breaks = 10) +
  xlab("Insert Size") + ylab("Number of samples") +
  ggtitle("Insert Size (dataframe3)") +
  theme(axis.text.x = element_text(angle = 45, hjust=1))
```

Contacts
-----------------------------------------------------
* Name: __Giulia Agostinetto__
* e-mail: __g.agostinetto@campus.unimib.it__
* Name: __Davide Bozzi__
* e-mail: __d.bozzi1@campus.unimib.it__

Please refer to both contacts for further information or issues related to the framework.
