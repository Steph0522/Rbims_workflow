#' @title Read the output of dbCAN
#' @description read_dbcan3 calculates the abundance of each Gene within the 
#' bins based on the dbCAN output from run_dbcan3 V3.0.6. 
#' @usage read_dbcan(dbcan_path)
#' @param dbcan_path a path where dbCAN output data are. They 
#' should have the extension overview.txt and all files in the path are the ones that
#' need to be read. Output data should have 6 columns with the bin names 
#' followed by the Genes obtained in every algorithm (HMMER,Hotpep,DIAMOND), 
#' column 'Signalp' indcating if a Peptide signal is found and a column 
#' '#ofTools" indicating the number of algorithms that found this Gene. 
#' @details This function is part of a package used for the analysis 
#' of bins metabolism.
#' @import dplyr tidyr readr stringr rlang tidyselect
#' @examples
#' \dontrun{
#' read_dbcan3("C:/Users/bins/")
#' }
#' @export

read_dbcan3<-function(dbcan_path){
  ruta_dbcan<-dbcan_path
  # Load all the data tables results ---------------------------------------####
  lapply_read_delim_bind_rows <- function(path, pattern = "*overview.txt"){
    files = list.files(path, pattern, full.names = TRUE)
    lapply(files, read.delim, check.names=F) %>% bind_rows()
  }
  dbcan_df<-suppressWarnings(lapply_read_delim_bind_rows(ruta_dbcan))
  # Reading data ----------------------------------------------------------####
  dbcan_df_format<- suppressWarnings(
    dbcan_df %>%
      filter( .data$`#ofTools` >1) %>%
      rename(Bin_name=.data$`Gene ID`) %>%
      mutate( hmmer2=str_replace_all(.data$HMMER, "[[:punct:]]",  "\t")) %>%
      separate(.data$hmmer2, c("dbNamesHMM"),  sep="\t") %>%
      mutate(ecami2=str_replace_all(.data$eCAMI, "[[:punct:]]", "\t")) %>%
      separate(.data$ecami2, c("dbNameseCAMI"), sep="\t") %>%
      mutate(diamond2=str_replace_all(.data$DIAMOND, "[[:punct:]]", "\t")) %>%
      separate(.data$diamond2, c("dbNamesdiamond"), sep="\t") %>%
      unite("dbCAN_names", .data$dbNamesHMM, .data$dbNameseCAMI,  .data$dbNamesdiamond, 
            sep="_", remove = F) %>%  
      mutate(dbCAN_names=str_replace_all(.data$dbCAN_names, "^_", "")) %>%
      separate(.data$dbCAN_names, c("dbCAN_names"), sep="_") %>%
      mutate(dbCAN = case_when(str_detect(.data$dbCAN_names, "CBM") ~ 
                                 "carbohydrate-binding module [CBM]",
                               str_detect(.data$dbCAN_names, "CE") ~ 
                                 "carbohydrate esterases [CEs]",
                               str_detect(.data$dbCAN_names, "GH") ~ 
                                 "glycoside hydrolases [GHs]",
                               str_detect(.data$dbCAN_names, "GT") ~ 
                                 "glycosyltransferases [GTs]",
                               str_detect(.data$dbCAN_names, "PL") ~ 
                                 "polysaccharide lyases [PLs]")) %>%
      mutate_if(is.character, str_trim) %>%
      dplyr::select(.data$Bin_name, .data$dbCAN_names,.data$Signalp ) %>%
      calc_abundance(analysis = "dbCAN") %>% 
      dplyr::select(-.data$Scaffold_name))
  # Menssage --------------------------------------------------------------####
  initial<-dim(dbcan_df)
  final<-dim(dbcan_df %>%
               filter( .data$`#ofTools` >1))
  signals<- dbcan_df %>% group_by(.data$Signalp) %>% count() 
  signals2<- dbcan_df_format %>% group_by(.data$Signalp) %>% count() 
  
  
  print(paste0("Input Genes = " , initial[1]))
  print(paste0("Remained Genes after filtering = " , final[1]))
  print(paste0("Percentage of genes remained = " , round(final[1]/initial[1]*100), "%"))
  print(paste0("Number of genes with signals = " , sum(signals[-1,]$n)))
  print(paste0("Number of genes with signals that passed filtering = " , sum(signals2[-1,]$n)))
  
  
  return(dbcan_df_format)
}



