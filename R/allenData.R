### Interface to all Allen Data sets

#' Loading Allen Data into R
#'
#' This function converts an unnormalized countmatrix of genes x samples into counts-per-million normalization
#' @param cortical_area The cortical area to load e.g. 'ALM'
#' @param species The species for which to load the data e.g. human or mouse
#' @param normalization Normalization method to use out of 'exon_rpkm' and 'exon+intron_cpm'
#' @param directory Directory where the data is located
#' @keywords Allen 
#' @export
#' @examples
#' loadAllenData()

loadAllenData = function(cortical_area, normalization = 'exon+intron_cpm', directory = '/nfs/team205/aa16/AllenData/'){
  
  if (!normalization %in% c('exon_rpkm', 'exon+intron_cpm'))
  {
    stop('No valid normalization argument supplied.')
  }
  
  dataList = list()
  coldataList = list()
  for (i in 1:length(cortical_area)){
    dataList[[i]] = read.delim(paste(directory, 'mouse_', cortical_area[i], '_2018-06-14_', normalization, '-matrix.csv', sep = ''), sep = ',')
    rownames(dataList[[i]]) = dataList[[i]][,1]
    dataList[[i]] = dataList[[i]][,2:dim(dataList[[i]])[2]]
    coldataList[[i]] = read.delim(paste(directory, 'mouse_', cortical_area[i], '_2018-06-14_samples-columns.csv', sep = ''), sep = ',')
    rowdata = read.delim(paste(directory, 'mouse_', cortical_area[i], '_2018-06-14_genes-rows.csv', sep = ''), sep = ',')
  }
  
  data = dataList[[1]]
  coldata = coldataList[[1]]
  if (length(dataList) > 1){
    for (i in 1:(length(dataList)-1)){
      data = cbind(data, dataList[[i+1]])
      coldata = rbind(coldata, coldataList[[i+1]])
    }
  }
  
  rm(dataList)
  rm(coldataList)
  
  celltypes = coldata[,'cluster']
  splitted = lapply(as.character(celltypes), function(x) strsplit(x, " "))
  firstEntry = unlist(lapply(splitted, function(x) x[[1]][1]))
  
  remove = (celltypes == 'High Intronic VISp L5 Endou' | firstEntry == 'Low' | firstEntry == 'Batch' | firstEntry == 'Doublet') # Remove low quality cells
  data = data[,!remove]
  coldata = coldata[!remove,]
  
  remove = (rowSums(data != 0) == 0) # Remove unexpressed genes
  data = data[!remove,]
  rowdata = rowdata[!remove,]
  
  return(list(data, coldata, rowdata))
}
