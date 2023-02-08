# for every gene on Chr1 in the yeast Saccharomyces Cerevisiae design a 90 bp repair template (45bp ib each end flanking the gene stitched together) that would lead to a clean deletion
#output as a table with the actual strings of DNA sequences
#keep upstream and downstream sequences separate, read 'C' or 'W' for direction 
library(BSgenome.Scerevisiae.UCSC.sacCer3)
library(RFLPtools)
library(seqinr)
library(intervals)
library(VariantAnnotation)
library(GenomicFeatures)
library(GenomicRanges)
library(Biostrings)
library(stringi)

#get the entire gene sequence as a BSgenome object 
genome <- getBSgenome("BSgenome.Scerevisiae.UCSC.sacCer3")
class(genome) # BSgenome object
providerVersion(genome)

#creates a DNAString object containing the whole chrI sequence 
chrI_sequence <- genome$chrI

#get information about chrI from the gff file
annotations <- file.choose()
read.delim(annotations, header=F, comment.char="#") -> gff

#generate table with just the genes for chrI
gff.genes <- gff[gff[,3]=="gene",]

#define a function that creates a repair template given a gene sequence
get_upstream_template <- function(chr, gene_start){

  #gets 45bp upstream of the gene
  
  upstream = substr(chr,(gene_start-45),(gene_start-1))
  
  d <- DNAString(upstream)
    
  return(upstream)
}

#define a function that creates a repair template given a gene sequence
get_downstream_template <- function(chr,gene_end){
  
  #gets 45bp downstream of the gene 
  
  downstream = substr(chr,gene_end+1,(gene_end+45))
  
  #joins the two to output a 90bp DNAstring *keep separate 
  
  d <- DNAString(downstream)
  
  return(d)
}
#test getting the deletion template
#get_deletion_template(chrI_sequence, 538, 792)


#get functions to make code within the loop easier 
get_gene_start <- function(gff_file,row){
  
  gene_start = (gff_file[row,4])
  
  return(gene_start)
}


get_gene_end <- function(gff_file,row)
{
  gene_end = (gff_file[row,5])
  
  return(gene_end)
}

get_gene_info <- function(gff_file,row)
{
  gene_info = (gff_file[row,9])
  
  return(gene_info)
}

#function to determine which direction sequence goes in, true if W and false if C
get_gene_dir <- function(gene_info)
{
  dir <- substring(gene_info,10,10)
  
  if(dir == 'W')
  {
    return (TRUE) 
  }
  else
  {
    return (FALSE)
  }
}

print(get_gene_dir(get_gene_info(gff.genes,3)))


#loop through every gene in and produce the DNAstrings, store them in a character vector
#NOTE: this function is currently working for W genes but not C genes 

c_up <- character()
c_down <- character()
info <- character()
index <- numeric()
i <- 1
while(gff.genes[i,1]=="chrI"){
  #get the deletion template for row i 
  gene_start <- get_gene_start(gff.genes, i)
  gene_end <- get_gene_end(gff.genes, i)
  upstream_template <- get_upstream_template(chrI_sequence,gene_start)
  downstream_template <- get_downstream_template(chrI_sequence,gene_end)
  
  #add gene information to vector 
  gene_info = get_gene_info(gff.genes,i)
  info = append(info, paste(gene_info),i)
  
  
   #add repair template to vector in correct direction
  if (get_gene_dir(gene_info)==TRUE){
    c_up = append(c_up, paste(upstream_template))
    c_down = append(c_down, paste(downstream_template))
  }
  
  else{

    #reverse sequences for the C direction genes 
    final_upstream <- DNAString(stri_reverse(upstream_template))
    
    final_downstream <- DNAString(stri_reverse(downstream_template))
    
    c_down = append(c_down, paste(final_upstream))
    
    
    c_up = append(c_up, paste(final_downstream))
  }
  
 
  index = append(index,i)
  
  i = i + 1 
}



#store vectors in a table 
df <- data.frame(index,c_up, c_down , info)
colnames(df) <- c('Gene Number', 'Repair Template Upstream', 'Repair Template Downstream', 'Gene Information')
                  



