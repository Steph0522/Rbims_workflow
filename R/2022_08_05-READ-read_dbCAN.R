#' @title Read the output of KofamScan/KofamKoala or KAAS.
#' @description read_ko calculates the abundance of each KO within the 
#' bins based on the KofamScan or KofamKoala output.
#' @usage read_ko(data_kofam=NULL, data_kaas=NULL, data_interpro=NULL)
#' @param data_kofam a path where KofamScan/KofamKoala output data are. They 
#' should have the extension .txt and all files in the path are the ones that
#' need to be read. Output data should have 5 columns with the bin names 
#' followed by the scaffold name divided by a '-' or '_': bin_scaffoldXX.
#' @param data_kaas a data frame with 2 columns. Contigs are expected to 
#' indicate in their names the bin name followed by the scaffold name 
#' divided by a '-' or '_': bin_scaffoldXX. 
#' @param data_interpro a data frame output of read_interpro. This
#' argument is used within mapping_KO.
#' @details This function is part of a package used for the analysis 
#' of bins metabolism.
#' @import dplyr tidyr readr stringr rlang
#' @importFrom utils read.table
#' @importFrom purrr map_dfr 
#' @examples
#' \dontrun{
#' read_ko("C:/Users/bins/")
#' }
#' @export

read_dbcan<-function(dbcan_path){
  ruta_dbcan<-dbcan_path
  # Load all the data tables results ---------------------------------------####
  lapply_read_delim_bind_rows <- function(path, pattern = "*overview.txt"){
    files = list.files(path, pattern, full.names = TRUE)
    lapply(files, read_delim, delim="\t") %>% bind_rows()
  }
  dbcan_df<-lapply_read_delim_bind_rows(ruta_dbcan)
  # Reading data ----------------------------------------------------------####
  dbcan_df_format<- suppressWarnings(
    dbcan_df %>%
      filter( `#ofTools` >1) %>%
      rename(Bin_name=`Gene ID`) %>%
      mutate( hmmer2=str_replace_all(HMMER, "[[:punct:]]",  "\t")) %>%
      separate(hmmer2, c("dbNamesHMM"),  sep="\t") %>%
      mutate(hotpep2=str_replace_all(Hotpep, "[[:punct:]]", "\t")) %>%
      separate(hotpep2, c("dbNameshotpep"), sep="\t") %>%
      mutate(diamond2=str_replace_all(DIAMOND, "[[:punct:]]", "\t")) %>%
      separate(diamond2, c("dbNamesdiamond"),msep="\t") %>%
      nite("dbCAN_names", dbNamesHMM, dbNameshotpep,  dbNamesdiamond, 
           sep="_", remove = F) %>%  
      mutate(dbCAN_names=str_replace_all(dbCAN_names, "^_", "")) %>%
      separate(dbCAN_names, c("dbCAN_names"), sep="_") %>%
      utate(dbCAN = case_when(str_detect(dbCAN_names, "CBM") ~ 
                                "carbohydrate-binding module [CBM]",
                              str_detect(dbCAN_names, "CE") ~ 
                                "carbohydrate esterases [CEs]",
                              str_detect(dbCAN_names, "GH") ~ 
                                "glycoside hydrolases [GHs]",
                              str_detect(dbCAN_names, "GT") ~ 
                                "glycosyltransferases [GTs]",
                              str_detect(dbCAN_names, "PL") ~ 
                                "polysaccharide lyases [PLs]"))  %>%
      mutate(across(where(is.character), str_trim))) %>%
    dplyr::select(Bin_name,dbCAN_names) %>%
    calc_abundance(analysis = "dbCAN") %>% 
    dplyr::select(-Scaffold_name)
  # Menssage --------------------------------------------------------------####
  initial<-dim(dbcan_df)
  final<-dim(dbcan_df %>%
               filter( `#ofTools` >1))
  signals<- dbcan_df %>% group_by(Signalp) %>% count()

  print(paste0("Input Genes = " , initial[1]))
  print(paste0("Output filtered Genes = " , final[1]))
  print(paste0("Percentage remained = " , round(final[1]/initial[1]*100), "%"))
  print(paste0("Number of Genes with no signals = " , signals[2]))

  return(dbcan_df_format)
}

