data_dbcan<-"inst/extdata/test_data/Prueba1_Bin_154_2_1.fna.faaoverview.txt"

dbcan_table<- data_dbcan
read_dbcan<-function(dbcan_table){
library("readr")
library(dplyr)
library(tidyr)
library(stringr)
#reading data
dbcan_df<-read.delim(dbcan_table, header = T, check.names = F)
dbcan_df_format<- suppressWarnings(dbcan_df %>%
  filter( `#ofTools` >1) %>%
  separate(.data$`Gene ID`, c("Bin_name", "Scaffold_name"),
           sep = "[_|-][s|S]caffold") %>%
  mutate(Scaffold_name = paste0( "scaffold", .data$Scaffold_name),
         .data$Scaffold_name) %>%
  unite("Scaffold_name", c("Bin_name", "Scaffold_name"), remove=FALSE) %>%
  mutate(hmmer2=str_replace_all(HMMER, "[[:punct:]]", "\t")) %>%
  separate(hmmer2, c("dbNamesHMM"), sep="\t") %>%
  mutate(hotpep2=str_replace_all(Hotpep, "[[:punct:]]", "\t")) %>%
  separate(hotpep2, c("dbNameshotpep"), sep="\t") %>%
  mutate(diamond2=str_replace_all(DIAMOND, "[[:punct:]]", "\t")) %>%
  separate(diamond2, c("dbNamesdiamond"), sep="\t") %>%
  unite("dbCAN_names", dbNamesHMM, dbNameshotpep, dbNamesdiamond, sep="_", remove = F) %>%
  mutate(dbCAN_names=str_replace_all(dbCAN_names, "^_", "")) %>%
  separate(dbCAN_names, c("dbCAN_names"), sep="_") %>%
  mutate(dbCAN = case_when(
    str_detect(dbCAN_names, "CBM") ~ "carbohydrate-binding module [CBM]",
    str_detect(dbCAN_names, "CE") ~ "carbohydrate esterases [CEs]",
    str_detect(dbCAN_names, "GH") ~ "glycoside hydrolases [GHs]",
    str_detect(dbCAN_names, "GT") ~ "glycosyltransferases [GTs]",
    str_detect(dbCAN_names, "PL") ~ "polysaccharide lyases [PLs]"
  ))   %>%
  mutate(across(where(is.character), str_trim)))

initial<-dim(dbcan_df)
final<-dim(dbcan_df %>%
      filter( `#ofTools` >1))
signals<- dbcan_df %>% group_by(Signalp) %>% count()

print(paste0("Input scaffolds = " , initial[1]))
print(paste0("Output filtered = " , final[1]))
print(paste0("Percentage remained = " , round(final[1]/initial[1]*100), "%"))
print(paste0("Number of scaffolds with no signals = " , signals[2]))
return(dbcan_df)
}

table_test<-read_dbcan(dbcan_table = data_dbcan)
