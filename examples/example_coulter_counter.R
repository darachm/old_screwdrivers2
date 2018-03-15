
# Load that function!
source("../scripts/coulterZ2reader.R")

# Read in all the Z2 files from a dir. Extension can be changed.
some_datar <- read_directory_of_Z2("data/example_coulter_counter/")

# To make things tidier later
library(tidyverse)

#This next part must be customized for your application! 
#Look at the comments!
datar_to_filter <- some_datar %>% as_tibble %>% 
  mutate(Filename=sub(".=#Z2","",basename(as.character(Path)))) %>% 
  select(-Path) %>% 
  separate(Filename,into=c("Date","Time","Treatment","FlaskOrNot"),
    sep="_")

# Now, how about filtering out the low end bubbles? Look at the blank
# controls. Lots of bubbles. Threshold is up to you.
datar_to_filter %>% 
  ggplot()+aes(x=BinCenter,y=BinCount,col=factor(Time))+
  geom_point()+
  facet_wrap(~FlaskOrNot)+
  scale_y_log10()

datar_to_filter %>% filter(BinCenter > 3) %>%
  ggplot()+aes(x=BinCenter,y=BinCount,col=factor(Time))+
  geom_point()+
  facet_wrap(~FlaskOrNot)+
  scale_y_log10()

datar_to_filter %>% filter(BinCenter > 2.5) %>%
  ggplot()+aes(x=BinCenter,y=BinCount,col=factor(Time))+
  geom_point()+
  facet_wrap(~FlaskOrNot)+
  scale_y_log10()

# Filter with that threshold
filtered_datar <- datar_to_filter %>% 
  filter(BinCenter > 2.5,FlaskOrNot=="flask") 

# Here, I expand the observations into a vector. Inefficient 
# computationally, easy to use/explain with other tools.
filtered_datar_expanded <- filtered_datar %>% nest(BinCenter,BinCount) %>%
  mutate(Observations=map(data,function(x){ 
      rep(x[["BinCenter"]],x[["BinCount"]]) 
    } ))

summary_stats <- filtered_datar_expanded %>%
  mutate(
    TotalCounts=unlist(map(Observations,length)),
    MeanDiameter=unlist(map(Observations,mean,na.rm=T)),
    MedianDiameter=unlist(map(Observations,mean,na.rm=T)),
    Minutes=as.numeric( 
      as.numeric(sub("(\\d\\d?)(\\d\\d)","\\1",Time))*60+
      as.numeric(sub("(\\d\\d?)(\\d\\d)","\\2",Time))
      )
    ) %>%
  select(-Observations,-data,-FlaskOrNot)
summary_stats

summary_stats %>% 
  gather(Variable,Value,TotalCounts,MeanDiameter,MedianDiameter)%>%
  ggplot()+aes(x=Minutes,y=Value,col=Treatment)+
  facet_wrap(~Variable,scales="free")+
  geom_point()

summary_stats %>% 
  ggplot()+aes(x=Minutes,y=TotalCounts,col=Treatment)+
  scale_y_log10()+
  geom_point()

