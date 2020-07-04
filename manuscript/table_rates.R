library("xtable")
library("tidyverse")

#CD4 progression
load("C:/Users/ahauser/Documents/Step1/Rfiles/used scripts/survivalanalysis/mat1.RData")
mat=mat1 %>%
  filter(treat %in% c("T1","F1","T2","F2")) %>%
  mutate(treat_1=substr(treat,1,1),
         treat_2=substr(treat,2,2)) %>%
  transmute(parameters = paste0("$1/\\nu_{CD4}^{",treat_1,"_",treat_2,"}$"),
            order=factor(treat,levels=unique(treat),labels=1:length(unique(treat))),
            cd4=cd4_in,
            description=paste0("Average time to progress from one to another CD4 class, at ", "$",treat_1,"_",treat_2,"$"),
            space=" ",
            est = paste0(round(as.numeric(median)/4.35,0),
                            " [", round(as.numeric(lower_b)/4.35,0),",",round(as.numeric(upper_b)/4.35,0),"]")) %>%
  spread(cd4,est) %>%
  arrange(order) %>%
  select(-order)
print(xtable(mat), sanitize.text.function = function(x) {x},include.rownames=FALSE)

#CD4 progression
load("C:/Users/ahauser/Documents/Step1/Rfiles/used scripts/survivalanalysis/mat2.RData")
mat=mat2 %>%
  filter(treat %in% c("T1","S1","T2","S2")) %>%
  mutate(treat_1=substr(treat,1,1),
         treat_2=substr(treat,2,2)) %>%
  transmute(parameters = paste0("$1/\\tilde{\\nu}_{CD4}^{",treat_1,"_",treat_2,"}$"),
            order=factor(treat,levels=unique(treat),labels=1:length(unique(treat))),
            cd4=cd4_in,
            description=paste0("Average time to progress from one to another CD4 class, at ", "$",treat_1,"_",treat_2,"$"),
            space=" ",
            est = paste0(round(as.numeric(median)/4.35,0),
                         " [", round(as.numeric(lower_b)/4.35,0),",",round(as.numeric(upper_b)/4.35,0),"]")) %>%
  spread(cd4,est) %>%
  arrange(order) %>%
  select(-order)
print(xtable(mat), sanitize.text.function = function(x) {x},include.rownames=FALSE)

#Care stages
load("C:/Users/ahauser/Documents/Step1/Rfiles/used scripts/survivalanalysis/mat5.RData")
mat=mat5 %>%
  filter(name %in% c("T1ToS1","T1ToF1","S1ToF1","S1ToF1","F1ToS1","F1ToT2","T2ToS2","T2ToF2","S2ToF2","F2ToS2")) %>%
  mutate(treat_1=substr(name,1,1),
         treat_2=substr(name,2,2),
         treat_3=substr(name,5,5),
         treat_4=substr(name,6,6)) %>%
  transmute(parameters = paste0("$1/\\gamma_{",treat_1,"_",treat_2,"\\rightarrow ",treat_3,"_",treat_4,"}$"),
            order=factor(name,levels=unique(name),labels=1:length(unique(name))),
            cd4=cd4,
            description=paste0("Time from ", "$",treat_1,"_",treat_2,"$", " to ", "$",treat_3,"_",treat_4,"$"),
            est = paste0("\\makecell{",round(as.numeric(lambda)/4.35,1),"\\\\[0cm] ",
                         " [", round(as.numeric(lambda_min)/4.35,1),",",round(as.numeric(lambda_max)/4.35,1),"]","}")) %>%
  spread(cd4,est) %>%
  arrange(order) %>%
  select(-order)
print(xtable(mat), sanitize.text.function = function(x) {x},include.rownames=FALSE)
