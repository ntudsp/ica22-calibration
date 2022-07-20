Replication code for: Assessment Of A Cost-Effective Headphone
Calibration Procedure For Soundscape Evaluations
================

## Data Preparation

The survey data was collected via a matlab GUI
(<https://github.com/kenowr/satp-gui>). The data from each participant
is stored in a single csv file. All participants’ data would be merged
into a single dataframe.

### Acoustic parameters

Extract single value

``` r
#load computed single values
data.sv<-read_excel("./Acoustic Params/acousticParams.xlsx",skip = 1) %>%
  `colnames<-`(sv.colnames)

#clean data to wide format
data.sv<-data.sv %>%  
  select(c("UCLFilename","Channel",params.udrTest)) %>% #select relevant
  mutate(Channel=ifelse(Channel=="Ch 3","Left", #rename channel
                        ifelse(Channel=="Ch 4","Right",Channel))) %>%
  mutate(method=ifelse(grepl('UCL',UCLFilename),"UCL", #add method
                       ifelse(grepl("NTU",UCLFilename),"NTU","Insitu"))) %>%
  mutate(UCLFilename=str_remove_all(UCLFilename, #remove all unnecessary
                                    "\\s*\\([^\\)]+\\)|(_UCL)|(_NTU)"))
```

### Listening test data

``` r
#merge all csv files containing subjective data based on UCL calibration
allcsvnames.UCL <- list.files(path = "./UCL Result",
                        pattern = "*.csv", full.names = TRUE)
subj.UCL <- ldply(allcsvnames.UCL,read.csv,header=FALSE)

#merge all csv files containing subjective data based on NTU calibration
allcsvnames.NTU <- list.files(path = "./NTU Result",
                        pattern = "*.csv", full.names = TRUE)
subj.NTU <- ldply(allcsvnames.NTU,read.csv,header=FALSE)

#extract participant ID
participantID.UCL <- as.numeric(gsub(".*?([0-9]+).*", 
                                 "\\1", allcsvnames.UCL))
participantID.NTU <- as.numeric(gsub(".*?([0-9]+).*", 
                                 "\\1", allcsvnames.NTU))

noOfStimuli <- 27 #27 stimuli from SATP project
participantIDvec.UCL <- rep(1:length(participantID.UCL), each=noOfStimuli)
participantIDvec.NTU <- rep(1:length(participantID.NTU), each=noOfStimuli)

stimuli.Name <- read.csv2(file = "stimuliIDkey.csv",sep = ",") %>%
  `colnames<-`(c("stimuliID","UCLFilename"))

#add participant ID col & stimuliName
subj.UCL<-cbind(participantIDvec.UCL,subj.UCL) %>% mutate(Calibration="UCL")
subj.NTU<-cbind(participantIDvec.NTU,subj.NTU) %>% mutate(Calibration="NTU") 

colnames(subj.UCL)<-c("participantID","stimuliID","pleasant", "chaotic",
                      "vibrant","uneventful","calm", "annoying",
                      "eventful","monotonous","check","duration","Calibration")
colnames(subj.NTU)<-c("participantID","stimuliID","pleasant", "chaotic",
                      "vibrant","uneventful","calm", "annoying",
                      "eventful","monotonous","check","duration","Calibration")

#remove outliers
subj.UCL<- subj.UCL %>%
  #merge(.,stimuli.Name) %>%
  filter(participantID!=4) %>%
  filter(participantID!=16)

subj.NTU<- subj.NTU %>%
  #merge(.,stimuli.Name) %>%
  filter(participantID!=1)

#Paired data
pairing.key<-read_delim("key_correspondence.txt",delim = "/",col_names = F) %>%
  `colnames<-`(c("NTU.ID","UCL.ID")) %>%
  mutate(UCL.ID=ifelse(UCL.ID=="No corresponding UCL",NA,UCL.ID)) %>% #remove unpaired
  filter(!NTU.ID==1) %>% #remove outlier
  merge(.,data.frame(UCL.ID=unique(subj.UCL$participantID)),all = TRUE)

#merge matched ID data from both NTU and UCL sets with corresponding key pairs

subj.comb<-rbind(subj.UCL %>% 
  filter(participantID %in% pairing.key$UCL.ID) %>%
  mutate(UCL.ID=participantID) %>%
  merge(.,pairing.key, by="UCL.ID"),subj.NTU %>% 
  filter(participantID %in% pairing.key$NTU.ID) %>%
  mutate(NTU.ID=participantID) %>%
  merge(.,pairing.key, by="NTU.ID"))  %>%
  select(!participantID)
```

## Analysis

### Single value acoustic parameters

``` r
defvar<-c("UCLFilename","method","Channel")
n.stimuliID<-27

#t-distribution at 95% Conf Lvl
td<-qt(0.975,df=n.stimuliID-1)

#generate stats for BA plots
data.sv.ba2<-data.sv %>%
  pivot_longer(!defvar,names_to = "Params",values_to = "Values") %>%
  pivot_wider(names_from = Channel,values_from = "Values") %>% #expand chn to mean 
  mutate(Values=rowMeans(select(.,c("Left","Right")),na.rm = T)) %>% #mean left & right
  select(!c("Left","Right")) %>% #remove channels
  pivot_wider(names_from = method, 
              values_from = "Values") %>% #expand method
  #compute diff in calib for each stimuli for all params
  dplyr::mutate(d_Insitu.UCL=Insitu-UCL,
                d_Insitu.NTU=Insitu-NTU,
                d_NTU.UCL=NTU-UCL,
                ave_Insitu.UCL=rowMeans(select(.,c("Insitu","UCL")),na.rm=TRUE),
                ave_Insitu.NTU=rowMeans(select(.,c("Insitu","NTU")),na.rm=TRUE),
                ave_NTU.UCL=rowMeans(select(.,c("NTU","UCL")),na.rm=TRUE)) %>%
  select(!c("Insitu","UCL","NTU")) %>% #remove channels
  dplyr::group_by(Params) %>% 
  #mean of differences
  dplyr::mutate(m_Insitu.UCL=mean(d_Insitu.UCL,na.rm=T),#mean of diff
                m.lwr_Insitu.UCL=m_Insitu.UCL-  #lwr 95% CI mean
                  td*sd(d_Insitu.UCL,na.rm = T)/sqrt(n.stimuliID), 
                m.upr_Insitu.UCL=m_Insitu.UCL+  #upr 95% CI mean
                  td*sd(d_Insitu.UCL,na.rm = T)/sqrt(n.stimuliID),
                m_Insitu.NTU=mean(d_Insitu.NTU,na.rm=T), #mean of diff
                m.lwr_Insitu.NTU=m_Insitu.NTU-  #lwr 95% CI mean
                  td*sd(d_Insitu.NTU,na.rm = T)/sqrt(n.stimuliID), 
                m.upr_Insitu.NTU=m_Insitu.NTU+  #upr 95% CI mean
                  td*sd(d_Insitu.NTU,na.rm = T)/sqrt(n.stimuliID),
                m_NTU.UCL=mean(d_NTU.UCL,na.rm=T), #mean of diff
                m.lwr_NTU.UCL=m_NTU.UCL-  #lwr 95% CI mean
                  td*sd(d_NTU.UCL,na.rm = T)/sqrt(n.stimuliID), 
                m.upr_NTU.UCL=m_NTU.UCL+  #upr 95% CI mean
                  td*sd(d_NTU.UCL,na.rm = T)/sqrt(n.stimuliID)) %>%
  #lower limits of agreement (LoA)
  dplyr::mutate(lwr_Insitu.UCL=m_Insitu.UCL- #lwr LoA
                  1.96*sd(d_Insitu.UCL,na.rm = T),
                lwr.lwr_Insitu.UCL=lwr_Insitu.UCL-  #lwr 95% CI lwr LoA
                  td*sqrt((1/n.stimuliID + 3.8416/(2*(n.stimuliID-1))))
                  *sd(d_Insitu.UCL,na.rm=T),
                lwr.upr_Insitu.UCL=lwr_Insitu.UCL+  #upr 95% CI lwr LoA
                  td*sqrt((1/n.stimuliID + 3.8416/(2*(n.stimuliID-1))))
                  *sd(d_Insitu.UCL,na.rm=T),
                lwr_Insitu.NTU=m_Insitu.NTU-   #lwr LoA
                  1.96*sd(d_Insitu.NTU,na.rm = T), 
                lwr.lwr_Insitu.NTU=lwr_Insitu.NTU-  #lwr 95% CI lwr LoA
                  td*sqrt((1/n.stimuliID + 3.8416/(2*(n.stimuliID-1))))
                  *sd(d_Insitu.NTU,na.rm=T),
                lwr.upr_Insitu.NTU=lwr_Insitu.NTU+  #upr 95% CI lwr LoA
                  td*sqrt((1/n.stimuliID + 3.8416/(2*(n.stimuliID-1))))
                  *sd(d_Insitu.NTU,na.rm=T),
                lwr_NTU.UCL=m_NTU.UCL - #lwr LoA
                  1.96*sd(d_NTU.UCL,na.rm = T),
                lwr.lwr_NTU.UCL=lwr_NTU.UCL-  #lwr 95% CI lwr LoA
                  td*sqrt((1/n.stimuliID + 3.8416/(2*(n.stimuliID-1))))
                  *sd(d_NTU.UCL,na.rm=T),
                lwr.upr_NTU.UCL=lwr_NTU.UCL+  #upr 95% CI lwr LoA
                  td*sqrt((1/n.stimuliID + 3.8416/(2*(n.stimuliID-1))))
                  *sd(d_NTU.UCL,na.rm=T)) %>%
  #upper limits of agreement (LoA)
  dplyr::mutate(upr_Insitu.UCL=m_Insitu.UCL+ #upr LoA
                  1.96*sd(d_Insitu.UCL,na.rm = T),
                upr.lwr_Insitu.UCL=upr_Insitu.UCL-  #lwr 95% CI upr LoA
                  td*sqrt((1/n.stimuliID + 3.8416/(2*(n.stimuliID-1))))
                  *sd(d_Insitu.UCL,na.rm=T),
                upr.upr_Insitu.UCL=upr_Insitu.UCL+  #upr 95% CI upr LoA
                  td*sqrt((1/n.stimuliID + 3.8416/(2*(n.stimuliID-1))))
                  *sd(d_Insitu.UCL,na.rm=T),
                upr_Insitu.NTU=m_Insitu.NTU+   #upr LoA
                  1.96*sd(d_Insitu.NTU,na.rm = T), 
                upr.lwr_Insitu.NTU=upr_Insitu.NTU-  #lwr 95% CI upr LoA
                  td*sqrt((1/n.stimuliID + 3.8416/(2*(n.stimuliID-1))))
                  *sd(d_Insitu.NTU,na.rm=T),
                upr.upr_Insitu.NTU=upr_Insitu.NTU+  #upr 95% CI upr LoA
                  td*sqrt((1/n.stimuliID + 3.8416/(2*(n.stimuliID-1))))
                  *sd(d_Insitu.NTU,na.rm=T),
                upr_NTU.UCL=m_NTU.UCL+ #upr LoA
                  1.96*sd(d_NTU.UCL,na.rm = T),
                upr.lwr_NTU.UCL=upr_NTU.UCL-  #lwr 95% CI upr LoA
                  td*sqrt((1/n.stimuliID + 3.8416/(2*(n.stimuliID-1))))
                  *sd(d_NTU.UCL,na.rm=T),
                upr.upr_NTU.UCL=upr_NTU.UCL+  #upr 95% CI upr LoA
                  td*sqrt((1/n.stimuliID + 3.8416/(2*(n.stimuliID-1))))
                  *sd(d_NTU.UCL,na.rm=T)) %>%
  ungroup() %>%
  pivot_longer(!c("UCLFilename","Params"),
               names_to = c(".value","set"),
               names_pattern = "(.+)_(.+)") %>% #match pattern to multiple var
  mutate(set=ifelse(set=="Insitu.UCL","In-situ~OCV",
                    ifelse(set=="Insitu.NTU","In-situ~HATS","HATS~OCV")))

#Bland-Altman (parametric)

params.LC<-c("L(C)","L5(C)","L10(C)","L50(C)","L90(C)","L95(C)")
data.LC=data.sv.ba2 %>% 
         filter(Params %in% params.LC)

params.LA<-c("L(A)","L5(A)","L10(A)","L50(A)","L90(A)","L95(A)")
data.LA=data.sv.ba2 %>% 
         filter(Params %in% params.LA)

params.N<-c("N5","N10","N50","N90","N95")
data.N=data.sv.ba2 %>% 
         filter(Params %in% params.N)

params.S<-c("S","S5","S10","S50","S90","S95")
data.S=data.sv.ba2 %>% 
         filter(Params %in% params.S)

params.R<-c("R","R5","R10","R50","R90","R95")
data.R=data.sv.ba2 %>% 
         filter(Params %in% params.R)

params.T<-c("Tone5","Tone10","Tone50","Tone90","Tone95")
data.T=data.sv.ba2 %>% 
         filter(Params %in% params.T)

params.F<-c("Fluc5","Fluc10","Fluc50","Fluc90","Fluc95")
data.F=data.sv.ba2 %>% 
         filter(Params %in% params.F)

df.Names <- list(data.LC=data.LC, data.LA=data.LA,
                 data.S=data.S, data.N=data.N,
                 data.R=data.R, data.T=data.T,
                 data.F=data.F)

params.list<- list(params.LC=params.LC,params.LA=params.LA,
                   params.S=params.S,params.N=params.N,
                   params.R=params.R,params.T=params.T,
                   params.F=params.F)

set1clr<-brewer.pal(n = 9,"Set1")
xlim.upp<-list(LC=120,LA=120,S=3,N=120,R=0.2,Tone=5,Fluc=0.4)
xlim.low<-list(LC=30,LA=30,S=0,N=0,R=-0.05,Tone=-.05,Fluc=-0.1)
ylim.upp<-list(LC=20,LA=20,S=1.2,N=35,R=0.06,Tone=2,Fluc=0.08)
geom.text.size<-4
theme.size<-16
geom.point.size<-4

#BA Plots
for (i in 1:length(df.Names)) {
  df<-df.Names[[i]]  
  g<-baplotOpts(df,X=ave,Y=d,color=set, #ggplot aes
     fg.XY=as.formula(set~Params), #facet grid formula
     np=length(params.list[[i]]),#no. of params
     #hline for mean diff, upper LoA, lower LoA
     hl.mean=m, hl.upper=upr, hl.lower=lwr, 
     #lower/upper 95% limits of lwr/upr LoA & mean
     lwr.ymin=lwr.lwr,lwr.ymax=lwr.upr,
     upr.ymin=upr.lwr,upr.ymax=upr.upr,
     m.ymin=m.lwr,m.ymax=m.upr,
     geom.text.size,theme.size,geom.point.size,colorp,
     xlim.low[[i]],xlim.upp[[i]],ylim.low=-ylim.upp[[i]],ylim.upp[[i]],
     ylabel="Difference between methods",
     xlabel="Arithemic average between methods")
  g
  ggsave(paste(i,"BAPlot.pdf"),plot = g, width = 1500, height = 800, units = "px",scale = 2.5)
}
```

    ## Warning: Removed 12 rows containing missing values (geom_point).
    ## Removed 12 rows containing missing values (geom_point).
    ## Removed 12 rows containing missing values (geom_point).

    ## Warning: Removed 10 rows containing missing values (geom_point).

    ## Warning: Removed 1 rows containing missing values (geom_rect).

    ## Warning: Removed 12 rows containing missing values (geom_point).

    ## Warning: Removed 10 rows containing missing values (geom_point).
    ## Removed 10 rows containing missing values (geom_point).

### Optimal pooled t-test

``` r
paq<-c("eventful","vibrant","pleasant","calm",
       "uneventful","monotonous","annoying","chaotic")

#initialise dataframe to store p.values of optt
optt.colnames<-c("stimuliID", "PAQ", "p.value","diff")
optt.pval.df<-data.frame(
  matrix(ncol=4,nrow=0, 
         dimnames=list(NULL, optt.colnames)))

#optimal pooled t-test for each stimuli and PAQ combination
for(st.ID in 1:n.stimuliID){
  for(paq.ID in paq){
    
    t<-subj.comb %>%
      filter(stimuliID==st.ID) %>%
      select(c("NTU.ID","UCL.ID",paq.ID,"Calibration")) %>%
      pivot_wider(names_from = Calibration,values_from = paq.ID)
    
    optt<-t.test.partial(as.data.frame(t[c("UCL","NTU")]))
    optt.pval.df<-rbind(optt.pval.df,
                      data.frame(st.ID,paq.ID,
                                 optt$p.value,optt$estimate,
                                 row.names = NULL) %>%
                        `colnames<-`(optt.colnames))
  }
}
optt.pval.df<-merge(optt.pval.df,stimuli.Name) %>%
  #Bonferroni correction
  dplyr::mutate(adj.pval=p.adjust(p = p.value,
                                  method = "bonferroni",
                                  n = n.stimuliID*length(paq))) %>%
  #compute ave abs diff for ranking
  dplyr::group_by(UCLFilename) %>% 
  dplyr::mutate(ave=mean(abs(diff))) %>% 
  dplyr::ungroup() %>%
  dplyr::mutate(UCLFilename=as.factor(UCLFilename)) %>%
  dplyr::mutate(UCLFilename=fct_reorder(UCLFilename,-ave))

optt.pval.df$stars <- cut(optt.pval.df$adj.pval, 
                     breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), 
                     label=c("***", "**", "*", ""))

#Plot heatmap of estimated mean of differences + signif stars
g<-ggplot(optt.pval.df,aes(UCLFilename,PAQ,fill=diff)) + 
  geom_tile() +
  scale_fill_distiller(palette = "PuOr",limits=c(-60,60)) +
  #scale_y_discrete(expand = c(0, 0)) +
  geom_text(aes(label=stars), color="black", 
            size=5,vjust = 0.75,angle = 30) +
  labs(fill='Mean of \ndifferences',
       y="ISO/TS 12913-2 Perceived affective quality attributes",
       x="Stimuli from SATP") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 30,vjust = 0.75, hjust=0.5),
        axis.text.y = element_text(angle = 30,vjust = 0.5, hjust=0.75)) 
g
```

<img src="main_files/figure-gfm/optt-1.png" height="150%" />

``` r
ggsave("optt.pdf",plot = g, width = 1800, height = 1000, units = "px",scale = 1.5)
```
