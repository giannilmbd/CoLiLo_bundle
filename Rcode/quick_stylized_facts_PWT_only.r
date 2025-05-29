{ 
  rm(list = setdiff(ls(), lsf.str()))

if(!require('IRkernel'))install.packages('IRkernel')
IRkernel::installspec()
options(repos=c("https://stat.ethz.ch/CRAN/","https://cran.r-project.org"))
# setRepositories( addURLs = c(ETH="https://stat.ethz.ch/CRAN/",
# DE="https://ftp.gwdg.de/pub/misc/cran/",
# FR="https://mirror.ibcp.fr/pub/CRAN/",
# IT="https://cran.stat.unipd.it/"),ind=c(1,2,3,4,5))
# IF problems with graphics API reinstall ragg
if(.Platform$OS.type=="unix"){
  setwd('~/Dropbox/projects/ANNAGIANCARLOGIANNI/CoLiLoGlobal/models/CoLiLoCodes4Circulation/Rcode/')
}else{
  setwd('C:/Users/gi003182/Dropbox (BIS)/Apps/ANNAGIANCARLO/CoLiLoGlobal/models/CoLiLoCodes4Circulation/Rcode/')
}
if(!require("jpeg",character.only=TRUE)) install.packages("jpeg", dependencies = TRUE)
library(openxlsx)

if(!require('multidplyr'))install.packages('multidplyr')
library(multidplyr)
library(dplyr, warn.conflicts = FALSE)
library(ggplot2)
library(stargazer)
library(plotly)
library(httpgd)
library(languageserver)
library(gridExtra)
if(!require("e1071")) install.packages("e1071")
library(e1071)
if(!require("plm")) install.packages("plm")
library(plm)
if(!require("e1071")) install.packages("e1071")
library(e1071)
if(!require('visdat'))install.packages("visdat")
library(visdat)
if(!require('UpSetR'))install.packages("UpSetR")
library(UpSetR)
if(!require('naniar'))install.packages("naniar")
library(naniar)
if(!require('tidyverse'))install.packages("tidyverse")
library(tidyverse)
if(!require('colorDF'))install.packages("colorDF")
library(colorDF)

options(vsc.use_httpgd = TRUE)
if(.Platform$OS.type=="unix"){
  source('~/Dropbox/Econometrics/R/Functions/PlotThemeScale/PlotThemeScale.r',echo=T)
  source("~/Dropbox/Econometrics/R/Functions/Filters/myfilter.r")
  source("~/Dropbox/Econometrics/R/Functions/Coverage/show_coverage.r")
}else{
  source('~/Dropbox/Econometrics/R/Functions/PlotThemeScale/PlotThemeScale.r',echo=T)
  source('~/Dropbox/Econometrics/R/Functions/Filters/myfilter.r')
  source("~/Dropbox/Econometrics/R/Functions/Coverage/show_coverage.r")
}
}
system('mkdir figures')
# remotes::install_github("nx10/unigd")
# library(unigd)
# X11()
# ugd()
# dev.off()
 # system('soffice --headless --convert-to xlsx:"Calc MS Excel 2007 XML" kappa_data_s.csv')
# system('soffice --headless --convert-to xlsx:"Calc MS Excel 2007 XML" kappa_data_v.csv')


if(!require("readstata13",character.only = T))install.packages("readstata13")


#----- Filter -----
filt<-function(data_=NULL,type_="BK",fr_=4,polyord=1){
    names(data_)<-c("year","variable");
                tmp<-try(myfilter(data_=data_,series_ = "variable",
                    lg_=0,fr_ = fr_,type_ = type_),silent = T);
                if(class(tmp)[1]!="try-error"){x<-tmp$cycle}
                else{
                  x<-NA};
                return(x)}

# estimate an AR(1) on the cycle

# Define a function to estimate AR(1) equation for each cc unit
ar1_estimator <- function(cc_data=NULL,var1='rgdp_cycle') {
  # Estimate AR(1) model
  form=formula(paste(var1,'~dplyr::lag(',var1,',1)',sep=''))
  ar1_model <- lm(form, data = cc_data)
  # Extract autocorrelation coefficient
  autocorr_coef <- coef(ar1_model)[2]
  # Extract standard deviation of innovation
  innovation_sd<- summary(ar1_model)$sigma

  # Calculate residuals
  residuals <- ar1_model$residuals
   # Extract mean
 mean_inno<-mean(residuals,na.rm=T)
  # Calculate skewness of residuals
  skewness <- skewness(residuals,na.rm=T)
  # Calculate kurtosis of residuals
  kurtosis <- kurtosis(residuals,na.rm=T)+3
  # SIMULATE MOMENTS OF AR PROCESS
  Time_10000;
  eps1=sample(residuals,T,replace=T)
  y=eps1[1];
  for (cnt in 2:T){
    y[cnt]<-autocorr_coef*y[cnt-1]+eps1[cnt]
  }
  y=y[ceiling(.1*T):T]
mean_ar<-mean(y)
var_ar<-var(y)
  # Calculate skewness of process
  skewness_ar <- skewness(y)
  # Calculate kurtosis of process
  kurtosis_ar <- kurtosis(y)+3
  # Return output as a data frame
  return(data.frame(cc = unique(cc_data$cc), autocorr_coef = autocorr_coef,mean=mean_inno,
   stdev = innovation_sd, skewness = skewness, kurtosis = kurtosis,
   mean_ar=mean_ar,var_ar=var_ar,skewness_ar=skewness_ar,kurtosis_ar=kurtosis_ar))
} # END ar_estimation function


# USING PEN WORLD TABLE -----------

if (!require('pwt10',character.only=TRUE))install.packages("pwt10")
library(pwt10)
pwt=pwt10.01
names(pwt)
if (!require('countrycode',character.only = T))install.packages('countrycode')
pwt$cc=pwt$isocode
pwt$date=pwt$year
# pwt<-pdata.frame(pwt,index=c("cc","year"),drop.NA.series = T)

# CHECK IF BALANCED -------------------------------------------------------

is.pbalanced(pwt)



# DEFINE GEOGRAPHICAL REGIONS ------------------------------------------
cntrcds<-openxlsx::read.xlsx(xlsxFile = '~/Dropbox/Econometrics/R/IMF Country Class/FMEconGroup.xlsx',sheet=4)
AEnms<-countrycode::countrycode(sourcevar=cntrcds$Advanced.Economies,origin='country.name',destination = 'iso3c')
AEnms<-AEnms[!is.na(AEnms)]
pwt<-pwt %>% group_by(cc) %>%
  mutate(AE=ifelse(cc[1]%in%AEnms,'AE','Other'))
EMnms<-countrycode::countrycode(sourcevar=cntrcds$Emerging.Market.Economies,origin='country.name',
                                destination = 'iso3c',custom_match = c('Micronesia'="FSM",'Kosovo'='XXK'))
EMnms<-EMnms[!is.na(EMnms)]
pwt<-pwt %>% group_by(cc) %>%
  mutate(AE=ifelse(cc[1]%in%EMnms,'EM',AE))

# on the basis of very irregular rgdpNA we remove the following
cnt_excl=c('YEM','UZB','UKR','TJK','TKM','SYR',
           'SXM','SRB','SLE','RUS','MSR','MNE','MKD','MDA','LVA','LTU','LBR',
           'LBN','KWT','KAZ','IRQ','HRV','GUY','GNQ','GEO','CUW','COD',
           'ARM','AZE','BLR','VEN','VGB')
# Ad-hoc subsample of countries

cntry_sel<-c("ABW" ,"AGO" ,"ALB" ,"ARE" ,"ARG" ,"ATG" ,"AUS" ,"AUT" ,"BEL" ,"BGR" ,"BHR" ,"BHS" ,"BIH" ,"BLZ" ,"BOL" ,"BRA" ,"BRB" ,"BRN" ,"BWA" ,"CAN" ,"CHE" ,"CHL" ,"CHN" ,"COL" ,"CPV" ,"CRI" ,"CYP",
 "CZE" ,"DEU" ,"DMA" ,"DNK" ,"DOM" ,"DZA" ,"ECU" ,"EGY" ,"ESP" ,"EST" ,"FIN" ,"FJI" ,"FRA" ,"GAB" ,"GBR" ,"GRC" ,"GRD" ,"GTM" ,"HKG" ,"HUN" ,"IDN" ,"IND" ,"IRL" ,"IRN" ,"ISL" ,"ISR" ,"ITA",
 "JAM" ,"JOR" ,"JPN" ,"KNA" ,"KOR" ,"LCA" ,"LKA" ,"LUX" ,"MAC" ,"MAR" ,"MDV" ,"MEX" ,"MLT" ,"MNG" ,"MUS" ,"MYS" ,"NAM" ,"NLD" ,"NOR" ,"NZL" ,"OMN" ,"PAK" ,"PAN" ,"PER" ,"PHL" ,"POL" ,"PRT",
 "PRY" ,"QAT" ,"ROU" ,"SAU" ,"SGP" ,"SLV" ,"SUR" ,"SVK" ,"SVN" ,"SWE" ,"SWZ" ,"SYC" ,"THA" ,"TTO" ,"TUN" ,"TUR" ,"TWN" ,"URY" ,"USA" ,"VCT" ,"ZAF")
inenglish<-countrycode::countrycode(sourcevar=cntry_sel,destination='country.name', origin= 'iso3c')
View(inenglish)
cntry_sel0<-c("Angola", "Argentina",
 "Australia", "Austria", "Belgium", "Bulgaria",
 "Bolivia",
"Brazil", "Canada", "Switzerland",
"Chile", "China", "Colombia", "Cyprus",
"Czechia", "Germany", "Denmark",
"Egypt", "Spain", "Estonia", "Finland",
"France",  "United Kingdom", "Greece",
"Hong Kong SAR China", "Hungary", "Indonesia", "India",
"Ireland", "Iran", "Iceland", "Israel", "Italy",
"Japan", "South Korea",
"Luxembourg", "Mexico",
"Malta", "Malaysia", "Namibia", "Netherlands",
"Norway", "New Zealand","Peru",
"Philippines", "Poland", "Portugal", "Paraguay", "Qatar", "Romania",
"Saudi Arabia", "Singapore", "Slovakia",
"Slovenia", "Sweden", "Thailand",
"Tunisia", "Turkey", "Taiwan", "Uruguay", "United States",
"South Africa")
cntry_sel<-countrycode::countrycode(sourcevar=cntry_sel0, origin='country.name', destination = 'iso3c')
# DEFINE WORLD AS Subset of AE+EM+couple of others
World<-cntry_sel

# TRUNCATE INITIAL DATE ALSO FOR FILTERS

pwt_sub <- pwt[!pwt$cc %in%cnt_excl&pwt$date>1959,]


pwt_sub <- pwt_sub[pwt_sub$cc%in% World,]

pwt_sub$year<-as.numeric(pwt_sub$date)
summary(pwt_sub)
# MISSING DATA -----
nas<-show_coverage(pwt_sub,cols_=c("rgdpna","pop"),id='cc',time_='year')
fob<-nas[[1]]
fob[,2]<-as.numeric(fob[,2])
fob[,3]<-as.numeric(fob[,3])
fob[,'dif']<-fob[,2]<fob[,3]
knitr::kable(fob)
knitr::kable(cbind(cntry_sel0,fob[,1:3]),format='latex',col.names=c('Country','ISO',"GDP",'population'),caption='Starting Observation')
# plot missing data
ggplot(pwt_sub,aes(year,rgdpna))+geom_point()+geom_miss_point()+facet_wrap(~cc,scales='free')


# export import share ----------
# transform share of import in GDP into share of inport in domesti absorption C+I+G
pwt_sub %<>% group_by(cc) %>%
  mutate(imp_sh=csh_m/(csh_c+csh_i+csh_g))

ggplot(pwt_sub)+
  geom_point(aes(year,imp_sh),shape=1)+
  geom_point(aes(year,csh_m),shape=2)+facet_wrap(~cc,scales='free_y')

# OECD MEMBERS -------
oecd_countries_df <- data.frame(
  country = c(
    "Australia", "Austria", "Belgium", "Canada", "Chile", "Colombia",
    "Costa Rica", "Czech Republic", "Denmark", "Estonia", "Finland",
    "France", "Germany", "Greece", "Hungary", "Iceland", "Ireland",
    "Israel", "Italy", "Japan", "Korea", "Latvia", "Lithuania",
    "Luxembourg", "Mexico", "Netherlands", "New Zealand", "Norway",
    "Poland", "Portugal", "Slovak Republic", "Slovenia", "Spain",
    "Sweden", "Switzerland", "Türkiye", "United Kingdom",
    "United States"
  ),
  iso3 = c(
    "AUS", "AUT", "BEL", "CAN", "CHL", "COL",
    "CRI", "CZE", "DNK", "EST", "FIN",
    "FRA", "DEU", "GRC", "HUN", "ISL", "IRL",
    "ISR", "ITA", "JPN", "KOR", "LVA", "LTU",
    "LUX", "MEX", "NLD", "NZL", "NOR",
    "POL", "PRT", "SVK", "SVN", "ESP",
    "SWE", "CHE", "TUR", "GBR",
    "USA"
  ),
  stringsAsFactors = FALSE
)

oecd_countries_df


# CREATE REGIONAL AGGREGATES ----------------
# per capita GDP
# Modified on 17 August 2023

pwt_sub<-pwt_sub %>% group_by(cc,year) %>%
  mutate(rgdp_pc=rgdpna/pop,
        pop_RoW=sum(pwt_sub[pwt_sub$year==year[1],'pop'],na.rm=T)-pop,
        pop_AE=sum(pwt_sub[pwt_sub$year==year[1]&pwt_sub$AE=='AE','pop'],na.rm=T)-ifelse(AE=='AE',pop,0),
        pop_EM=sum(pwt_sub[pwt_sub$year==year[1]&pwt_sub$AE=='EM','pop'],na.rm=T)-ifelse(AE=='EM',pop,0),

         rgdp_RoW=sum(pwt_sub[pwt_sub$year==year[1],'rgdpna'],na.rm=T)/pop_RoW,
         rgdp_AE=(sum(pwt_sub[pwt_sub$year==year[1]&pwt_sub$AE=='AE','rgdpna'],na.rm=T)-ifelse(AE=='AE',rgdpna,0))
        /pop_AE,
         rgdp_EM=(sum(pwt_sub[pwt_sub$year==year[1]&pwt_sub$AE=='EM','rgdpna'],na.rm=T)-ifelse(AE=='EM',rgdpna,0))
        /pop_EM,
        # TFP
        rtfp_AE=(sum(pwt_sub[pwt_sub$year==year[1]&pwt_sub$AE=='AE','rtfpna'],na.rm=T)-
                   ifelse(AE=='AE',rtfpna,0)),
        rtfp_EM=(sum(pwt_sub[pwt_sub$year==year[1]&pwt_sub$AE=='EM','rtfpna'],na.rm=T)-
                   ifelse(AE=='EM',rtfpna,0)),
        rtfp_RoW=(sum(pwt_sub[pwt_sub$year==year[1]&pwt_sub$AE=='RoW','rtfpna'],na.rm=T)-
                   ifelse(AE=='RoW',rtfpna,0))
  )



# pwt_sub<-pwt_sub %>% group_by(cc,year) %>%
#   mutate(rgdp_pc=rgdpna/pop,
#         rgdp_RoW=sum(pwt_sub[pwt_sub$cc!=cc[1]&pwt_sub$year==year[1],'rgdpna'],na.rm=T)/sum(pwt_sub[pwt_sub$cc!=cc[1]&pwt_sub$year==year[1],'pop'],na.rm=T),
#
#
#          rgdp_AE=sum(pwt_sub[pwt_sub$cc!=cc[1]&pwt_sub$year==year[1]&pwt_sub$AE=='AE','rgdpna'],na.rm=T)/sum(pwt_sub[pwt_sub$cc!=cc[1]&pwt_sub$year==year[1]&pwt_sub$AE=='AE','pop'],na.rm=T),
#
#          rgdp_EM=sum(pwt_sub[pwt_sub$cc!=cc[1]&pwt_sub$year==year[1]&pwt_sub$AE=='EM','rgdpna'],na.rm=T)/sum(pwt_sub[pwt_sub$cc!=cc[1]&pwt_sub$year==year[1]&pwt_sub$AE=='EM','pop'],na.rm=T)
#
#         )

pwt_sub<-pwt_sub %>% group_by(cc,year) %>%
  mutate(rgdp_DEU=pwt_sub[pwt_sub$cc=='DEU'&pwt_sub$year==year[1],'rgdpna']$rgdpna/pwt_sub[pwt_sub$cc=='DEU'&pwt_sub$year==year[1],'pop']$pop,
         rgdp_USA=pwt_sub[pwt_sub$cc=='USA'&pwt_sub$year==year[1],'rgdpna']$rgdpna/pwt_sub[pwt_sub$cc=='USA'&pwt_sub$year==year[1],'pop']$pop,
         rtfp_DEU=pwt_sub[pwt_sub$cc=='DEU'&pwt_sub$year==year[1],'rtfpna']$rtfpna,
         rtfp_USA=pwt_sub[pwt_sub$cc=='USA'&pwt_sub$year==year[1],'rtfpna']$rtfpna,

         pop_DEU=pwt_sub[pwt_sub$cc=='DEU'&pwt_sub$year==year[1],'pop']$pop,
         pop_USA=pwt_sub[pwt_sub$cc=='USA'&pwt_sub$year==year[1],'pop']$pop
  )




pwt_sub<-as.data.frame(pwt_sub)

summary(pwt_sub)
# NOTE THAT rgdp_OTHER might be NA if the OTHERS are too few-small
# ggplot(pwt_sub)+geom_line(aes(year,rgdp_pc))+facet_wrap(~cc,scales='free_y')
#
# ggplot(pwt_sub)+geom_line(aes(year,rgdp_AE))+facet_wrap(~cc,scales='free_y')
# ggplot(pwt_sub)+geom_line(aes(year,rgdp_EM))+facet_wrap(~cc,scales='free_y')

# ggplot(pwt_sub)+geom_line(aes(year,rgdp_Other))+facet_wrap(~cc,scales='free_y')

cntrs<-as.character(pwt_sub$cc%>%unique())



# detrend -------------------------
filtr="PS" # "Qt"
pwt_sub<-pwt_sub %>% group_by(cc) %>%
           mutate(
             rgdp_cycle=filt(data_=data.frame(date,L=log(rgdp_pc),fr_=1),type_=filtr),
             rtfp_cycle=filt(data_=data.frame(date,L=log(rtfpna),fr_=1),type_=filtr),

             rgdp_cycle_AE=filt(data_=data.frame(date,L=log(rgdp_AE),fr_=1),type_=filtr),
             rtfp_cycle_AE=filt(data_=data.frame(date,L=log(rtfp_AE),fr_=1),type_=filtr),
             rgdp_cycle_EM=filt(data_=data.frame(date,L=log(rgdp_EM),fr_=1),type_=filtr),
             rtfp_cycle_EM=filt(data_=data.frame(date,L=log(rtfp_EM),fr_=1),type_=filtr),

             # rgdp_cycle_Other=filt(data_=data.frame(date,L=log(rgdp_Other),fr_=1),type_=filtr),
             rgdp_cycle_DEU=filt(data_=data.frame(date,L=log(rgdp_DEU),fr_=1),type_=filtr),
             rtfp_cycle_DEU=filt(data_=data.frame(date,L=log(rtfp_DEU),fr_=1),type_=filtr),
             rgdp_cycle_USA=filt(data_=data.frame(date,L=log(rgdp_USA),fr_=1),type_=filtr),
             rtfp_cycle_USA=filt(data_=data.frame(date,L=log(rtfp_USA),fr_=1),type_=filtr),

             rgdp_cycle_RoW=filt(data_=data.frame(date,L=log(rgdp_RoW),fr_=1),type_=filtr),
             rtfp_cycle_RoW=filt(data_=data.frame(date,L=log(rtfp_RoW),fr_=1),type_=filtr)
           )






# ggplot(pwt_subt)+geom_line(aes(year,rgdp_cycle))+facet_wrap(~cc,scales='free_y')
# ggplot(pwt_subt)+geom_line(aes(year,rgdp_cycle_AE))+facet_wrap(~cc,scales='free_y')
# ggplot(pwt_subt)+geom_line(aes(year,rgdp_cycle_EM))+facet_wrap(~cc,scales='free_y')
# # ggplot(pwt_subt)+geom_line(aes(year,rgdp_cycle_Other))+facet_wrap(~cc,scales='free_y')
# ggplot(pwt_subt)+geom_line(aes(year,rgdp_cycle_RoW))+facet_wrap(~cc,scales='free_y')


ncntr<-pwt_sub$cc %>% unique() %>% length()
nblks=ceiling(ncntr/30)
# {i=1
# while(i+30*(i-1)<=ncntr){
#     a=i+30*(i-1)
#     b=min(c(30*i,ncntr))
#
#   p<-ggplot(pwt_sub[pwt_sub$cc%in%cntrs[a:b],])+
#   geom_line(aes(x=date,y=log(rgdpna),color='rgdp'),linetype='solid',linewidth=2)+
#   geom_line(aes(x=date,y=log(rgdpna)-rgdp_cycle,color='trend'))+
#   facet_wrap(~cc,scale='free_y')+
#   labs(title=paste(c('from ',cntrs[a],' to ',cntrs[b]),collapse=''))
#   print(p)
#   i=i+1
#
# }
# }

# Compute the moments of the cycle
{
  # STDEV
mgdp_sd<-aggregate(data=pwt_sub,c(rgdp_cycle)~cc,sd)
names(mgdp_sd)<-c('cc','stdev')

mgdp_sd_AE<-aggregate(data=pwt_sub,c(rgdp_cycle_AE)~cc,sd)
names(mgdp_sd_AE)<-c('cc','stdev_AE')

mgdp_sd_EM<-aggregate(data=pwt_sub,c(rgdp_cycle_EM)~cc,sd)
names(mgdp_sd_EM)<-c('cc','stdev_EM')

# mgdp_sd_Other<-aggregate(data=pwt_sub,c(rgdp_cycle_Other)~cc,sd)
# names(mgdp_sd_Other)<-c('cc','stdev_Other')

mgdp_sd_RoW<-aggregate(data=pwt_sub,c(rgdp_cycle_RoW)~cc,sd)
names(mgdp_sd_RoW)<-c('cc','stdev_RoW')

# SKEWNESS

mgdp_s<-aggregate(data=pwt_sub,c(rgdp_cycle)~cc,skewness)
names(mgdp_s)<-c('cc','skewness')
mgdp_s_AE<-aggregate(data=pwt_sub,c(rgdp_cycle_AE)~cc,skewness)
names(mgdp_s_AE)<-c('cc','skewness_AE')

mgdp_s_EM<-aggregate(data=pwt_sub,c(rgdp_cycle_EM)~cc,skewness)
names(mgdp_s_EM)<-c('cc','skewness_EM')

# mgdp_s_Other<-aggregate(data=pwt_sub,c(rgdp_cycle_Other)~cc,skewness)
# names(mgdp_s_Other)<-c('cc','skewness_Other')

mgdp_s_RoW<-aggregate(data=pwt_sub,c(rgdp_cycle_RoW)~cc,skewness)
names(mgdp_s_RoW)<-c('cc','skewness_RoW')

# kurtosis
mgdp_k<-aggregate(data=pwt_sub,c(rgdp_cycle)~cc,kurtosis,type=1)
names(mgdp_k)<-c('cc','kurtosis')
mgdp_k_AE<-aggregate(data=pwt_sub,c(rgdp_cycle_AE)~cc,kurtosis,type=1)
names(mgdp_k_AE)<-c('cc','kurtosis_AE')

mgdp_k_EM<-aggregate(data=pwt_sub,c(rgdp_cycle_EM)~cc,kurtosis,type=1)
names(mgdp_k_EM)<-c('cc','kurtosis_EM')

# mgdp_k_Other<-aggregate(data=pwt_sub,c(rgdp_cycle_Other)~cc,kurtosis,type=1)
# names(mgdp_k_Other)<-c('cc','kurtosis_Other')

mgdp_k_RoW<-aggregate(data=pwt_sub,c(rgdp_cycle_RoW)~cc,kurtosis,type=1)
names(mgdp_k_RoW)<-c('cc','kurtosis_RoW')
}


{
allmoms<-full_join(mgdp_sd,mgdp_sd_AE,by='cc')
allmoms<-full_join(allmoms,mgdp_sd_EM,by='cc')
# allmoms<-full_join(allmoms,mgdp_sd_Other,by='cc')
allmoms<-full_join(allmoms,mgdp_sd_RoW,by='cc')
allmoms<-full_join(allmoms,mgdp_s,by='cc')
allmoms<-full_join(allmoms,mgdp_s_AE,by='cc')
allmoms<-full_join(allmoms,mgdp_s_EM,by='cc')
# allmoms<-full_join(allmoms,mgdp_s_Other,by='cc')
allmoms<-full_join(allmoms,mgdp_s_RoW,by='cc')
allmoms<-full_join(allmoms,mgdp_k,by='cc')
allmoms<-full_join(allmoms,mgdp_k_AE,by='cc')
allmoms<-full_join(allmoms,mgdp_k_EM,by='cc')
# allmoms<-full_join(allmoms,mgdp_k_Other,by='cc')
allmoms<-full_join(allmoms,mgdp_k_RoW,by='cc')
}
allmoms$cc<-as.character(allmoms$cc)

head(allmoms)
# ADD SINGLE BASE COUNTRIES DEU AND USA
allmoms$stdev_DEU=allmoms[allmoms$cc=='DEU','stdev']
allmoms$skewness_DEU=allmoms[allmoms$cc=='DEU','skewness']
allmoms$kurtosis_DEU=allmoms[allmoms$cc=='DEU','kurtosis']

allmoms$stdev_USA=allmoms[allmoms$cc=='USA','stdev']
allmoms$skewness_USA=allmoms[allmoms$cc=='USA','skewness']
allmoms$kurtosis_USA=allmoms[allmoms$cc=='USA','kurtosis']

# AR actual regressions -----------
# Use sapply to estimate AR(1) equation for each cc unit and store output in a data frame
library(parallel)
library(multidplyr)
ggplot(pwt_sub)+geom_line(aes(year,rgdp_cycle_EM))+facet_wrap(~cc,scales='free_y')
if(T==T){ #parallel
  cl <- detectCores()-1
cluster <-new_cluster(cl)
groups<-pwt_sub %>% group_by(cc)  %>% partition(cluster=cluster) # UNCOMMENT FOR PARRALLEL
groups
cluster_copy(cluster,"ar1_estimator")
cluster_library(cluster,c('dplyr','e1071','stats','base'))
ar_output_gdp <- groups %>%
  do(ar1_estimator(.,var1='rgdp_cycle')) %>%
  ungroup() %>% collect()
ar_output_AE <- groups %>%
  do(ar1_estimator(.,var1='rgdp_cycle_AE')) %>%
  ungroup() %>% collect()

ar_output_EM <- groups %>%
  do(ar1_estimator(.,var1='rgdp_cycle_EM')) %>%
  ungroup() %>% collect()

# ar_output_Other <- groups %>%
#   do(ar1_estimator(.,var1='rgdp_cycle_Other')) %>%
#   ungroup() %>% collect()

ar_output_RoW <- groups %>%
  do(ar1_estimator(.,var1='rgdp_cycle_RoW')) %>%
  ungroup() %>% collect()

ar_output_DEU <- groups %>%
  do(ar1_estimator(.,var1='rgdp_cycle_DEU')) %>%
  ungroup() %>% collect()

ar_output_USA <- groups %>%
  do(ar1_estimator(.,var1='rgdp_cycle_USA')) %>%
  ungroup() %>% collect()

rm(cluster)
}else{


ar_output_gdp <-pwt_sub %>%  group_by(cc) %>%
  do(ar1_estimator(.,var1='rgdp_cycle')) %>%
  ungroup()


ar_output_AE <-pwt_sub %>%  group_by(cc) %>%
  do(ar1_estimator(.,var1='rgdp_cycle_AE')) %>%
  ungroup()

ar_output_EM <-pwt_sub %>%  group_by(cc) %>%
  do(ar1_estimator(.,var1='rgdp_cycle_EM')) %>%
  ungroup()

# ar_output_Other <-pwt_sub %>%  group_by(cc) %>%
#   do(ar1_estimator(.,var1='rgdp_cycle_Other')) %>%
#   ungroup()

ar_output_RoW <-pwt_sub %>%  group_by(cc) %>%
  do(ar1_estimator(.,var1='rgdp_cycle_RoW')) %>%
  ungroup()

ar_output_DEU <- groups %>%
  do(ar1_estimator(.,var1='rgdp_cycle_DEU')) %>%
  ungroup() %>% collect()

ar_output_USA <- groups %>%
  do(ar1_estimator(.,var1='rgdp_cycle_USA')) %>%
  ungroup() %>% collect()
}

# exclude lower and upper 5% do this on the full data to ensure same list of countries
trm_ser<-function(data=NULL,var=NULL){
  data[,'nm_str_ser']<-data[,var]
  lwr=quantile(data$nm_str_ser,p=c(0.05,0.95))
  out<-data[data$nm_str_ser>=lwr[1]&data$nm_str_ser<=lwr[2],]
  out<-out[,names(out)!='nm_str_ser']
  return(out)
}


if(T==F){#MUTED this created extra columns by country. Instead what is needed is extra rows (lots)
# ar_output_gdp_trim<-trm_ser(data=ar_output_gdp,var='kurtosis_ar')
 ncolumns<-ncol(ar_output_gdp)
# nmscntr<-unique(ar_output_gdp_trim$cc)
# ar_output_AE_trim<-trm_ser(data=ar_output_AE,var='kurtosis_ar')
names(ar_output_AE)[2:ncolumns]<-paste(names(ar_output_AE)[2:ncolumns],'_AE',sep='')
# ar_output_EM_trim<-trm_ser(data=ar_output_EM,var='kurtosis_ar')
names(ar_output_EM)[2:ncolumns]<-paste(names(ar_output_EM)[2:ncolumns],'_EM',sep='')
# ar_output_Other_trim<-trm_ser(data=ar_output_Other,var='kurtosis_ar')
# names(ar_output_Other)[2:ncolumns]<-paste(names(ar_output_Other)[2:ncolumns],'_Other',sep='')
# ar_output_RoW_trim<-trm_ser(data=ar_output_RoW,var='kurtosis_ar')
names(ar_output_RoW)[2:ncolumns]<-paste(names(ar_output_RoW)[2:ncolumns],'_RoW',sep='')
nmscntr<-unique(ar_output_gdp$cc)
ar_output<-full_join(ar_output_gdp,ar_output_AE[ar_output_AE$cc%in%nmscntr,],by='cc')
ar_output<-full_join(ar_output,ar_output_EM[ar_output_EM$cc%in%nmscntr,],by='cc')
# ar_output<-full_join(ar_output,ar_output_Other[ar_output_Other$cc%in%nmscntr,],by='cc')
ar_output<-full_join(ar_output,ar_output_RoW[ar_output_RoW$cc%in%nmscntr,],by='cc')
}

# ADD POPULATION
{
pwt_pop <- aggregate(data=pwt_sub,cbind(pop,pop_AE,pop_EM,pop_RoW,pop_DEU,pop_USA)~cc,FUN=function(x)last(x))

ar_output_gdp<-full_join(ar_output_gdp,pwt_pop[,c('cc','pop')],by='cc')
summary(ar_output_gdp)
ar_output_AE<-full_join(ar_output_AE,pwt_pop[,c('cc','pop_AE')],by='cc')
names(ar_output_AE)[ncol(ar_output_AE)]<-c("pop")
ar_output_EM<-full_join(ar_output_EM,pwt_pop[,c('cc','pop_EM')],by='cc')
names(ar_output_EM)[ncol(ar_output_EM)]<-c("pop")
ar_output_RoW<-full_join(ar_output_RoW,pwt_pop[,c('cc','pop_RoW')],by='cc')
names(ar_output_RoW)[ncol(ar_output_RoW)]<-c("pop")
ar_output_DEU<-full_join(ar_output_DEU,pwt_pop[,c('cc','pop_DEU')],by='cc')
names(ar_output_DEU)[ncol(ar_output_DEU)]<-c("pop")
ar_output_USA<-full_join(ar_output_USA,pwt_pop[,c('cc','pop_USA')],by='cc')
names(ar_output_USA)[ncol(ar_output_USA)]<-c("pop")
}


# Add import shares -----------------------------------------------------------------------------------------------

# use average 2000-2007 & 2010-2019

pwt_sub_m<-pwt_sub[(pwt_sub$year>=2000&pwt_sub$year<=2007)|(pwt_sub$year>=2010&pwt_sub$year<=2019),]
unique(pwt_sub_m$year)
pwt_m <- aggregate(data=pwt_sub_m,cbind(imp_sh)~cc,FUN=function(x)mean(x,na.rm=T))
# create cc extension for regions
{
  ar_output_AE_cc<-ar_output_AE
  ar_output_AE_cc$cc<-paste(ar_output_AE$cc,'_A',sep='')
  ar_output_EM_cc<-ar_output_EM
  ar_output_EM_cc$cc<-paste(ar_output_EM$cc,'_E',sep='')
  ar_output_RoW_cc<-ar_output_RoW
  ar_output_RoW_cc$cc<-paste(ar_output_RoW$cc,'_R',sep='')
  ar_output_DEU_cc<-ar_output_DEU
  ar_output_DEU_cc$cc<-paste(ar_output_DEU$cc,'_D',sep='')
  ar_output_USA_cc<-ar_output_USA
  ar_output_USA_cc$cc<-paste(ar_output_USA$cc,'_U',sep='')

}
ar_output<-rbind(ar_output_gdp,ar_output_AE_cc,ar_output_EM_cc,ar_output_RoW_cc,ar_output_DEU_cc,ar_output_USA_cc)


# ar_output_trim<-trm_ser(data=ar_output,var='kurtosis_ar')
ar_output_trim<-ar_output
summary(ar_output_trim)



knitr::kable(ar_output_trim,format='pipe')



write.csv(ar_output_trim,file='ar_data_pwt.csv',row.names=F)
ncol(ar_output_trim)
names(ar_output_trim)
nrow(ar_output_trim)
