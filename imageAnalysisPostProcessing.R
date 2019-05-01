#!/usr/local/bin/Rscript

#  Created by Joseph Negri, original version 7/7/17.
#  copyright is maintained and must be perserved.
#  The work is provided as is; no warranty is provided, and users accept all liability.

# script for post processing of Image.csv files produced by Cell Profiler
# pipeline punctaCounting.cpproj.

# this script accepts 2 argument: 
# *an experiment or screen id 
# *path to the dir containing the *Image.csv CellProfiler output files
# *path to output directory
# if these are not provided, the basename and
# path of the current directory are used

# this script is designed to analyze data from 384well plates
# in which negative controls are positioned in columns: 1, 23
# and positive controls are positioned in columns: 2, 24

# the script will return:
# *.csv of well level, nuclei and puncta count data
# *.csv of Z-prime per assay plate
# *figures for evaluting assay performance:
#    *puncta count row-wise
#    *puncta count column-wise
#    *puncta count plate view
#    *normalized puncta count plate view
#    *nuclei count plate view
#    *control histogram

# Load essential libraries----
library(plyr)
library(ggplot2)
library(data.table)

# Set constants----
experiment <- basename(getwd())

args = commandArgs(trailingOnly=TRUE)

if(is.null(args[1])){
  experiment <- basename(getwd())
  imageDirPath <- getwd()
  outputDirPath <- getwd()
}

experiment <- args[1]

imageDirPath <- args[2]

outputDirPath <- args[3]

#  input per image data from CellProfiler output----
image.files <- 
  dir(
    path = imageDirPath,
    pattern = "*Image.csv",
    full.names = TRUE
    )

image.data <- rbindlist(lapply(image.files,function(x){
  fread(dir(pattern = x))
}))

cols <- c("ImageNumber","AreaOccupied_AreaOccupied_Puncta","Count_Nuclei",
          "Count_Puncta","Metadata_Plate","Metadata_Replicate",
          "Metadata_Well","Metadata_Row","Metadata_Column","Metadata_Field")

image.data <- image.data[,.SD,,.SDcols=cols]
setnames(image.data,
         c("AreaOccupied_AreaOccupied_Puncta",
           "Count_Nuclei","Count_Puncta","ImageNumber","Metadata_Column",
           "Metadata_Field","Metadata_Plate","Metadata_Row","Metadata_Well",
           "Metadata_Replicate"),
         c("area_occupied_puncta", "count_nuclei","count_puncta",
           "image_number","well_column","field","plate","well_row","well","replicate"))

#  Aggregate per image data to per well data----
well.data <- 
  image.data[
    ,lapply(.SD, sum, na.rm = TRUE),
    by = .(plate,replicate,well,well_row,well_column),
    .SDcols = c("area_occupied_puncta","count_nuclei","count_puncta")
  ]

well.data[,`:=`(
  plate=as.numeric(plate),
  replicate=as.factor(replicate),
  plate.rep=as.factor(paste0(plate,"_",replicate)),
  well_row=factor(well_row, levels = LETTERS[16:1]),
  well_column=as.numeric(well_column)
)]

# Define Neg and Pos Controls----
well.data[well_column %in% c(1,23),type:="N"]
well.data[well_column %in% c(2,24),type:="P"]
well.data[well_column %in% c(3:22),type:="C"]

# Normalize count_puncta to Neg and Pos Controls----
sapply(unique(well.data$plate.rep),function(x){
         negcon.median <- well.data[plate.rep==x][
           type=="N",.(median.negcon=median(count_puncta))]
         poscon.median <- well.data[plate.rep==x][
           type=="P",.(median.negcon=median(count_puncta))]
         well.data[plate.rep==x,
                   norm_count_puncta:=(count_puncta-negcon.median)/(negcon.median-poscon.median),
                   by=well]
         return(NULL)
         }
       )

#  Reassign data types----
well.data[,well_type:=as.factor(type)]
well.data[,well_row:=factor(well_row, levels = LETTERS[16:1])]

#  Write well.data to csv----
write.csv(
  well.data,
  file=paste0(
    outputDirPath,
    Sys.Date(),"_",
    experiment,
    "_well-stats.csv")
  )

#  Generate ICCB submission files----
lapply(unique(well.data$plate),function(x){
  cnt.nuclei.wide <- dcast.data.table(
    well.data[plate==x], plate + well + type ~ plate.rep,
    value.var = "count_nuclei"
  )
  colnames(cnt.nuclei.wide)[4:5] <- c("count_nuclei_A","count_nuclei_B")
  
  cnt.puncta.wide <- dcast.data.table(
    well.data[plate==x], plate + well + type ~ plate.rep,
    value.var = "count_puncta"
  )
  colnames(cnt.puncta.wide)[4:5] <- c("count_puncta_A", "count_puncta_B")
  
  norm.cnt.puncta.wide <- dcast.data.table(
    well.data[plate==x], plate + well + type ~ plate.rep,
    value.var = "norm_count_puncta"
  )
  colnames(norm.cnt.puncta.wide)[4:5] <- c("norm_count_puncta_A", "norm_count_puncta_B")
  
  plate.well.wide <- cbind(
    cnt.nuclei.wide,
    cnt.puncta.wide[,c(4:5),with=FALSE],
    norm.cnt.puncta.wide[,c(4:5),with=FALSE]
  )
  
  write.csv(
    plate.well.wide,
    file = paste0(
      outputDirPath,
      experiment,"_",x,".csv")
    )
  return(NULL)
})

#  Calculate stats per treatment, conc----
treatment.mean <- well.data[
  ,lapply(.SD,mean),by=.(plate.rep,type),
  .SDcols=c("count_nuclei","count_puncta","norm_count_puncta")]
setnames(treatment.mean,
         c("count_nuclei","count_puncta","norm_count_puncta"),
         c("mean_count_nuclei","mean_count_puncta","mean_norm_count_puncta")
          )

treatment.sd <- well.data[
  ,lapply(.SD,sd),by=.(plate.rep,type),
  .SDcols=c("count_nuclei","count_puncta","norm_count_puncta")]
setnames(treatment.sd,
         c("count_nuclei","count_puncta","norm_count_puncta"),
         c("sd_count_nuclei","sd_count_puncta","sd_norm_count_puncta")
          )

treatment.stats <- merge(treatment.mean,treatment.sd,by=c("plate.rep","type"))

z.prime.puncta.count <- sapply(unique(treatment.stats$plate.rep),function(x){
NegCon.mean <- treatment.stats[plate.rep==x][type=="N",mean_count_puncta]
PosCon.mean <- treatment.stats[plate.rep==x][type=="P",mean_count_puncta]
NegCon.sd <- treatment.stats[plate.rep==x][type=="N",sd_count_puncta]
PosCon.sd <- treatment.stats[plate.rep==x][type=="P",sd_count_puncta]
z.prime <- 1 - (3*(NegCon.sd+PosCon.sd)/(NegCon.mean-PosCon.mean))
}
)

z.prime.plate <- data.table(
  plate.rep=unique(treatment.stats$plate.rep),
  z.prime.puncta.count=z.prime.puncta.count
  )

write.csv(
  z.prime.plate, 
  paste0(
    outputDirPath,
    Sys.Date(),
    "_z-prime-per-plate.csv"
    ),
  row.names = FALSE
  )

#  Figures-----

#  Plate Effect: Puncta Count across Rows
plate.effect.well.row <- ggplot(
  data=well.data,
  aes(x=well_row,y=count_puncta,
      group=well_row,color=type)
)

plate.effect.well.row+
  geom_point()+
  scale_x_discrete(limits=LETTERS[1:16])+
  facet_wrap(~plate.rep,ncol=2,labeller = label_both)+
  scale_color_brewer(palette = "Set1")+
  theme_bw()+
  ggtitle("Plate Effect: Well Row")

ggsave(
  filename=paste0(outputDirPath,Sys.Date(),"_plate.effect.well-row.pdf"),
  width = 40,
  height = 40,
  units = "cm",
  device = pdf
  )

#  Plate Effect: Puncta Count across Columns
plate.effect.well.col <- ggplot(
  data=well.data,
  aes(x=well_column,y=count_puncta,
      group=well_column,color=type)
)

plate.effect.well.col+
  geom_boxplot()+
  scale_x_discrete(limits=c(1:24))+
  facet_wrap(~plate.rep,ncol=2,labeller = label_both)+
  scale_color_brewer(palette = "Set1")+
  theme_bw()+
  ggtitle("Plate Effect: Well Column")

ggsave(
  filename=paste0(outputDirPath,Sys.Date(),"_plate.effect.well-column.pdf"),
  width = 40,
  height = 40,
  units = "cm",
  device = pdf
)

#  Plate Effect: Puncta Count Plateview
puncta.count.plateview <- ggplot(
  well.data,
  aes(x=well_column, y=well_row, fill=count_puncta))

puncta.count.plateview+
  geom_tile()+
  scale_fill_gradient2(low = "blue", mid ="white", high = "red",
                       midpoint = well.data[type=="N",median(count_puncta)])+
  scale_x_continuous(breaks = c(1:24))+
  theme_bw()+
  facet_wrap(~plate.rep,ncol=2,labeller = label_both)+
  ggtitle("Puncta Count - Plate View")

ggsave(
  filename=paste0(outputDirPath,Sys.Date(),"_puncta.count.plateview.pdf"),
  width = 40,
  height = 40,
  units = "cm",
  device = pdf
)

#  Plate Effect: Puncta Count Norm Plateview
norm.puncta.count.plateview <- ggplot(
  well.data,
  aes(x=well_column, y=well_row, fill=norm_count_puncta))

norm.puncta.count.plateview+
  geom_tile()+
  scale_fill_gradient2(low = "blue", mid ="white", high = "red",
                       midpoint = well.data[type=="N",median(norm_count_puncta)])+
  scale_x_continuous(breaks = c(1:24))+
  theme_bw()+
  facet_wrap(~plate.rep,ncol=2,labeller = label_both)+
  ggtitle("Normalized Puncta Count - Plate View")

ggsave(
  filename=paste0(outputDirPath,Sys.Date(),"_norm.puncta.count.plateview.pdf"),
  width = 40,
  height = 40,
  units = "cm",
  device = pdf
)

#  Plate Effect: Nuclei Count Plateview
nuclei.count.plateview <- ggplot(
  well.data,
  aes(x=well_column, y=well_row, fill=count_nuclei))

nuclei.count.plateview+
  geom_tile()+
  scale_fill_gradient2(low = "blue", mid ="white", high = "red",
                       midpoint = median(well.data$count_nuclei))+
  scale_x_continuous(breaks = c(1:24))+
  theme_bw()+
  facet_wrap(~plate.rep,ncol=2,labeller = label_both)+
  ggtitle("Nuclei Count - Plate View")

ggsave(
  filename=paste0(outputDirPath,Sys.Date(),"_nuclei.count.plateview.pdf"),
  width = 40,
  height = 40,
  units = "cm",
  device = pdf
)

#  Pos/Neg Con Histogram
con.histogram <- 
  ggplot(
  well.data[type!="C"],
  aes(x=count_puncta, 
      #color=well_type,
      fill=type)
  )

con.histogram+
  geom_histogram(binwidth = 50)+
  facet_wrap(~plate.rep,ncol=2,labeller = label_both)+
  theme_bw()+
  scale_fill_manual(values = c("#377eb8","#4daf4a"))+
  ggtitle(paste0(experiment, "\n", "Puncta Count - Distribution of Control Wells"))

ggsave(
  filename=paste0(outputDirPath,Sys.Date(),"_puncta.count.cntrl.dist.pdf"),
  width = 40,
  height = 40,
  units = "cm",
  device = pdf
)