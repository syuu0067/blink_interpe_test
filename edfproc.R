library(eyelinkReader)
library(tidyverse)
library(patchwork)

dat<-eyelinkReader::read_edf('blinktest.edf', import_samples=TRUE)

gorgx<-ggplot(data=dat$samples)+geom_line( aes(x=time_rel, y=gxR), color='red')
gorgx<-gorgx+ylim(dat$display_coords[1], dat$display_coords[3])
gorgy<-ggplot(data=dat$samples)+geom_line( aes(x=time_rel, y=gyR), color='blue')
gorgy<-gorgy+ylim(dat$display_coords[2], dat$display_coords[4])

gppl<-ggplot(data=dat$samples)+geom_line(aes(x=time_rel, y=paR), color='green')

gorgx / gorgy / gppl

tmpwid<-100
tblink<-1
prd<-(dat$blinks$sttime_rel[tblink]/(1000/dat$headers$rec_sample_rate)-tmpwid):(dat$blinks$entime_rel[tblink]/(1000/dat$headers$rec_sample_rate)+tmpwid)
plot(prd,dat$samples$paR[prd])

