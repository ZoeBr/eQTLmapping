# eQTLmapping

Calculating molecular population statistics
```
salloc --nodes=5 --ntasks-per-node=5 --cpus-per-task=4 --mem=150G --job-name=copy --time=24:00:00 --partition=general --account=a_ortiz_barrientos_coe srun --export=PATH,TERM,HOME,LANG --pty /bin/bash -l
module load anaconda3
conda activate pixy
pixy --stats pi dxy --vcf /scratch/user/s4480088/wgs_dataset_sf1.vcf.gz --populations /scratch/user/s4480088/Populations1.txt --window_size 10000 --n_cores 10 --chromosomes scaffold_1 --output_prefix 'wgs_test1' --output_folder '/scratch/user/s4480088'
```
Plotting pixy output
```
#In R
library(ggplot2)
library(readr)
library(dplyr)
library(tidyr)
install.packages('zoo')
library(zoo)

inp<-read.table("wgs_test1_pi.txt",sep="\t",header=T)

chroms <- unique(inp$chromosome)
chrOrder <- sort(chroms)
inp$chrOrder <- factor(inp$chromosome,levels=chrOrder)

if("avg_pi" %in% colnames(inp)){
  pops <- unique(inp$pop)
  for (p in pops){
    thisPop <- subset(inp, pop == p)
    popPlot <- ggplot(thisPop, aes(window_pos_1, avg_pi, color=chrOrder)) +
      geom_point()+
      facet_grid(. ~ chrOrder)+
      labs(title=paste("Pi for pop", p))+
      labs(x="position of window start", y="Pi")+
      scale_colour_manual(values=rep(c("black","gray"),ceiling((length(chrOrder)/2))))+
      theme_classic()+
      theme(legend.position = "none")
    ggsave(paste("piplot_",p, ".png", sep=""), plot = popPlot, device = "png", dpi = 300)
  }
} else {
  print("Pi not found in this file")
}

#Just WAT1 region

if("avg_pi" %in% colnames(inp)){
  pops <- unique(inp$pop)
  for (p in pops){
    thisPop <- subset(inp, pop == p & 
                        chromosome == "scaffold_1" & 
                        window_pos_1 >= 35848322 & 
                        window_pos_2 <= 35891354)
    
    popPlot <- ggplot(thisPop, aes(window_pos_1, avg_pi, color = chromosome)) +
      geom_point() +
      facet_grid(. ~ chromosome) +
      labs(title = paste("Pi for pop", p),
           x = "Position of window start",
           y = "Pi") +
      scale_colour_manual(values = rep(c("black", "gray"), ceiling((length(chrOrder) / 2)))) +
      theme_classic() +
      theme(legend.position = "none")
    
    ggsave(paste("piplot_", p, ".png", sep = ""), plot = popPlot, device = "png", dpi = 300)
  }
} else {
  print("Pi not found in this file")
}

#whole chromosome line graph - this is good code

pixy_to_long <- function(pixy_files){
  pixy_df <- list()
  for(i in 1:length(pixy_files)){
    stat_file_type <- gsub(".*_|.txt", "", pixy_files[i])
    if(stat_file_type == "pi"){
      df <- read_delim(pixy_files[i], delim ="\t")
      df <- df %>%
        gather(-pop, -window_pos_1, -window_pos_2, -chromosome, key = "statistic", value = "value") %>%
        rename(pop1 = pop) %>%
        mutate(pop2 = NA)
      pixy_df[[i]] <- df
      
    } else{
      df <- read_delim(pixy_files[i], delim = "\t")
      df <- df %>%
        gather(-pop1, -pop2, -window_pos_1, -window_pos_2, -chromosome, key = "statistic", value = "value")
      pixy_df[[i]] <- df
    }
  }
  bind_rows(pixy_df) %>%
    arrange(pop1, pop2, chromosome, window_pos_1, statistic)
}

pixy_folder <- "wgs_test1"
pixy_files <- "wgs_test1_pi.txt"
#pixy_files <- list.files(pixy_folder, pattern = "\\.txt$", full.names = TRUE)
pixy_df <- pixy_to_long(pixy_files)

pixy_labeller <- as_labeller(c(avg_pi = "pi"),
                               default = label_parsed)
pixy_df %>%
  filter(chromosome == "scaffold_1") %>%
  filter(statistic %in% c("avg_pi")) %>%
  mutate(chr_position = ((window_pos_1 + window_pos_2)/2)/100000) %>%
  ggplot(aes(x = chr_position, y = value, colour = statistic))+
  geom_line(size = 0.25)+
  facet_grid(statistic ~ .,
             scales = "free_y", switch = "x", space = "free_x",
             labeller = labeller(statistix = pixy_labeller,
                                 value = label_value))+
  xlab("Position on Chromosome (Mb)")+
  ylab("statistic value")+
  theme_bw()+
  theme(panel.spacing = unit(0.1, "cm"),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.position = "none")+
  scale_x_continuous(expand = c(0, 0))+
  scale_y_continuous(expand = c(0, 0))+
  scale_colour_brewer(palette = "Set1")

#Just WAt1 region
pixy_df %>%
  filter(chromosome == "scaffold_1") %>%
  filter(statistic %in% c("avg_pi")) %>%
  filter(window_pos_1 >= 35600000, window_pos_2 <= 35961354) %>%  # 
  mutate(chr_position = ((window_pos_1 + window_pos_2)/2)/1e6) %>%  # Mb for x-axis
  ggplot(aes(x = chr_position, y = value, colour = statistic)) +
  geom_line(size = 0.25) +
  facet_grid(statistic ~ .,
             scales = "free_y", switch = "x", space = "free_x",
             labeller = labeller(statistic = pixy_labeller, value = label_value)) +
  xlab("Position on Chromosome (Mb)") +
  ylab("Statistic Value") +
  theme_bw() +
  theme(panel.spacing = unit(0.1, "cm"),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.position = "none") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_colour_brewer(palette = "Set1")


#Attempt to smooth
smoothed_df <- pixy_df %>%
  filter(chromosome == "1", statistic == "avg_pi") %>%  # adjust as needed
  arrange(chr_position) %>%  # make sure data is sorted
  mutate(
    rolling_value = rollmean(value, k = 10, fill = NA, align = "center"),  # window of 10
    chr_position = ((window_pos_1 + window_pos_2)/2)/100000
  )

smoothed_df <- pixy_df %>%
  filter(chromosome == "scaffold_1", statistic == "avg_pi") %>%  # adjust as needed
  mutate(chr_position = ((window_pos_1 + window_pos_2)/2)/100000) %>%
  arrange(chr_position) %>%  # now this works!
  mutate(
    rolling_value = rollmean(value, k = 10, fill = NA, align = "center")
  )

ggplot(smoothed_df, aes(x = chr_position, y = rolling_value, colour = statistic)) +
  geom_line(size = 0.5) +
  facet_grid(statistic ~ .,
             scales = "free_y", switch = "x", space = "free_x") +
  xlab("Position on Chromosome (Mb)") +
  ylab("Rolling Mean of Statistic") +
  theme_bw() +
  theme(panel.spacing = unit(0.1, "cm"),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.position = "none") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_colour_brewer(palette = "Set1")
```
