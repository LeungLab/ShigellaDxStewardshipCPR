#######################################################
# Code for predicting who benefits from microbiologic testing
#######################################################
rm(list=ls())
graphics.off()


library(tidyverse)
library(lubridate)
library(gridExtra)
library(pROC)
library(data.table)
library(fields)
library(zoo)
library(ks)
library(KernSmooth)
library(ranger)
library(viridis)
library(purrr)
library(broom)
library(profvis)
library(furrr)
library(mice)
library(glmnet)
library(glmnetUtils)
library(cvAUC)
library(table1)
library(data.table)


setwd('')

################### FUNCTion for var screening, regression fitting, AUC calc ####
CPR.funct <- function(data,outcome,iter,nvars_opts){
  out=ranger(as.formula(paste(outcome,'~',paste(names,collapse="+"),sep="")),data=data,num.trees=1000,importance="impurity")
  imps=importance(out)
  df_imps_full=data.frame(names=names(imps),var_red=as.numeric(imps)) %>% arrange(desc(var_red))
  
  
  result=data.frame(iter=NA,nvar=NA,true=NA,pred_glm=NA,pred_RF=NA)
  test_record <- NA
  train_record <- NA
  
  for (each in 1:iter){
    print(each)
    train=data %>% sample_frac(.80,replace=F)
    
    test=data[-which(data$index %in% train$index),]
    
    train_record <- c(train_record,table(train[,outcome])[["1"]])
    test_record <- c(test_record,table(test[,outcome])[["1"]])
    
    
    out=ranger(as.formula(paste(outcome,'~',paste(names,collapse="+"),sep="")),data=train,num.trees=1000,importance="impurity")
    df_imps=data.frame(names=names(ranger::importance(out)),imps=ranger::importance(out)) %>% arrange(desc(imps))
    for (nvars in nvars_opts){
      
      print(nvars)
      out1=glm(as.formula(paste(outcome,'~',paste(df_imps$names[1:nvars],collapse="+"),sep="")),data=train,family="binomial",control=glm.control(maxit=50))
      out2=ranger(as.formula(paste(outcome,'~',paste(df_imps$names[1:nvars],collapse="+"),sep="")),data=train,num.trees=1000)
      
      df=data.frame(iter=each,nvar=nvars,true=test[[outcome]],pred_glm=as.numeric(predict(out1,newdata=test,type="response")),pred_RF=as.numeric(predict(out2,data=test,type="response")$predictions))
      result=rbind(result,df)
    }
  }
  result<-result[-1,]
  train_record<-train_record[-1]
  test_record<-test_record[-1]
  
  AUCs<-result %>% split(.$nvar) %>% purrr::map(~ci.cvAUC(.$pred_glm,.$true,folds=.$iter))
  AUCs2<-result %>% split(.$nvar) %>% purrr::map(~ci.cvAUC(.$pred_RF,.$true,folds=.$iter))
  
  AUC_df<-rbind(bind_rows(AUCs %>% purrr::map(~data.frame(AUC=.$cvAUC,SE=.$se,lower=.$ci[1],upper=.$ci[2],level=.$confidence,Model="LR"))),
                bind_rows(AUCs2 %>% purrr::map(~data.frame(AUC=.$cvAUC,SE=.$se,lower=.$ci[1],upper=.$ci[2],level=.$confidence,Model="RF"))))
  AUC_df$nvar<-rep(nvars_opts,2)
  AUC_df
  
  
  calib_fits=data.frame(nvar=NA,iter=NA,intc=NA,intc_LCI=NA,intc_UCI=NA,slope=NA,slope_LCI=NA,slope_UCI=NA)
  for (nvars in nvars_opts){
    for (each in 1:iter){
      data.temp <- result %>% filter(nvar==nvars & iter==each)
      
      intercept <- glm(true~1,offset=log(pred_glm/(1-pred_glm)),family="binomial",data=data.temp)
      slope <- glm(true~log(pred_glm/(1-pred_glm)),family="binomial",data=data.temp)
      
      df=data.frame(nvar=nvars,iter=each,
                    intc=coef(intercept),intc_LCI=confint(intercept)[1],intc_UCI=confint(intercept)[2],
                    slope=coef(slope)[2],slope_LCI=confint(slope)[2,1],slope_UCI=confint(slope)[2,2])
      calib_fits=rbind(calib_fits,df)
      
    }
  }
  calib_fits<-calib_fits[-1,]
  calib <- calib_fits %>% group_by(nvar) %>% summarize(mean(intc),mean(intc_LCI),mean(intc_UCI),
                                                       mean(slope),mean(slope_LCI),mean(slope_UCI))
  names(calib) <- c("nvar","intc","intc_LCI","intc_UCI","slope","slope_LCI","slope_UCI")  #renaming
  
  decilesCC <- result %>% split(.,list(.$iter,.$nvar),drop=TRUE) %>% purrr::map(. %>% arrange(pred_glm)) %>% #now have a df for each iteration of each nvar
    purrr::map(~mutate(.x, decile_glm=ntile(pred_glm,10))) %>% #create predicted glm decile groups; equivalent to: purrr::map(list_resB, ~mutate(.x, decile_glm=ntile(pred_glm,10))); str(temp3)
    bind_rows(.) %>% split(.,f=.$nvar) %>% #a list of df for each nvar which contains all iter for that nvar. "nest" might be better for this
    purrr::map(., . %>% group_by(decile_glm) %>% summarize(mean(true),mean(pred_glm))) #for each decile in each nvar, have an avg true and avg predicted
  
  output<-list(df_imps=df_imps_full,result=result,train_record=train_record,test_record=test_record,AUC_df=AUC_df,decilesCC=decilesCC,calib=calib,iter=iter,nvars_opts=nvars_opts)
  
}

################### FUNCTion for Se/Sp loop to generate data for plotting ####
SeSp.funct <- function(data,outcome,predictors,new.data,POC.Se,POC.Sp,min,max,increment){
  twovar_fit <- glm(as.formula(paste(outcome,'~',paste(predictors,collapse="+"),sep="")),data=data,family="binomial")
  val_data <- new.data
  val_data$predicted <- predict(twovar_fit, newdata=val_data,type="response")
  val_data <- val_data %>% select(all_of(append(x=predictors, values=c(outcome,"predicted"))))
  
  results=data.frame(lower_cutoff=NA,upper_cutoff=NA,
                     low.pred.prob.truly.neg=NA,low.pred.prob.truly.pos=NA,
                     mid.pred.prob.truly.neg=NA,mid.pred.prob.truly.pos=NA,
                     high.pred.prob.truly.neg=NA, high.pred.prob.truly.pos=NA)
  
  for (i in seq(min,max-increment,increment)){ 
    for (j in seq(i+increment,max,increment)){ 
      
      low.pred.prob.truly.neg <- length(val_data[[outcome]][val_data$predicted<i & val_data[[outcome]]==0])
      low.pred.prob.truly.pos <- length(val_data[[outcome]][val_data$predicted<i & val_data[[outcome]]==1])
      mid.pred.prob.truly.neg <- length(val_data[[outcome]][val_data$predicted>=i & val_data$predicted<=j & val_data[[outcome]]==0])
      mid.pred.prob.truly.pos <- length(val_data[[outcome]][val_data$predicted>=i & val_data$predicted<=j & val_data[[outcome]]==1])
      high.pred.prob.truly.neg <- length(val_data[[outcome]][val_data$predicted>j & val_data[[outcome]]==0])
      high.pred.prob.truly.pos <- length(val_data[[outcome]][val_data$predicted>j & val_data[[outcome]]==1])
      
      df=data.frame(lower_cutoff=i,upper_cutoff=j,
                    low.pred.prob.truly.neg=low.pred.prob.truly.neg,low.pred.prob.truly.pos=low.pred.prob.truly.pos,
                    mid.pred.prob.truly.neg=mid.pred.prob.truly.neg,mid.pred.prob.truly.pos=mid.pred.prob.truly.pos,
                    high.pred.prob.truly.neg=high.pred.prob.truly.neg, high.pred.prob.truly.pos=high.pred.prob.truly.pos)
      results=rbind(results,df)
      
      
    }
  }
  results=results[-1,]
  
  results <- results %>% mutate(TN = low.pred.prob.truly.neg + (mid.pred.prob.truly.neg * POC.Sp),
                                TP = high.pred.prob.truly.pos + (mid.pred.prob.truly.pos * POC.Se),
                                FN = low.pred.prob.truly.pos + (mid.pred.prob.truly.pos * (1-POC.Se)),
                                FP = high.pred.prob.truly.neg + (mid.pred.prob.truly.neg * (1-POC.Sp)),
                                
                                prop.tested=round(((mid.pred.prob.truly.neg+mid.pred.prob.truly.pos)/dim(val_data)[1])*100,2),
                                
                                #of regimen as a whole, incorporating Se/Sp of the POC test itself
                                Se = round(TP/(TP+FN),2), #test positive | truly positive
                                Sp = TN/(TN+FP), #test negative | truly negative
                                PPV = TP/(TP+FP), #truly positive | test positive
                                NPV = TN/(TN+FN), #truly negative | test negative
                                FPR = FP/(FP+TN), #test positive | truly negative
                                FNR = FN/(TP+FN) #test negative | truly positive
                                
  )
  
  #Se/Sp of bloody diarrhea as only test
  TN.bld = length(val_data[[outcome]][val_data$f4a_drh_blood==0 & val_data[[outcome]]==0])
  TP.bld = length(val_data[[outcome]][val_data$f4a_drh_blood==1 & val_data[[outcome]]==1])
  FN.bld = length(val_data[[outcome]][val_data$f4a_drh_blood==0 & val_data[[outcome]]==1])
  FP.bld = length(val_data[[outcome]][val_data$f4a_drh_blood==1 & val_data[[outcome]]==0])
  
  Se.bld = round(TP.bld/(TP.bld+FN.bld),2) #test positive | truly positive
  Sp.bld = round(TN.bld/(TN.bld+FP.bld),2) #test negative | truly negative
  PPV.bld = round(TP.bld/(TP.bld+FP.bld),2) #truly positive | test positive
  NPV.bld = round(TN.bld/(TN.bld+FN.bld),2) #truly negative | test negative
  FPR.bld = round(FP.bld/(FP.bld+TN.bld),2) #test positive | truly negative
  FNR.bld = round(FN.bld/(TP.bld+FN.bld),2) #test negative | truly positive
  
  results <- results %>% mutate(TN.pct.chng = ((TN-TN.bld)/TN.bld)*100, #percent change
                                TP.pct.chng = ((TP-TP.bld)/TP.bld)*100,
                                FN.pct.chng = ((FN-FN.bld)/FN.bld)*100,
                                FP.pct.chng = ((FP-FP.bld)/FP.bld)*100,
                                Se.pct.chng = ((Se-Se.bld)/Se.bld)*100,
                                Sp.pct.chng = ((Sp-Sp.bld)/Sp.bld)*100,
                                PPV.pct.chng = ((PPV-PPV.bld)/PPV.bld)*100,
                                NPV.pct.chng = ((NPV-NPV.bld)/NPV.bld)*100,
                                FPR.pct.chng = ((FPR-FPR.bld)/FPR.bld)*100,
                                FNR.pct.chng = ((FNR-FNR.bld)/FNR.bld)*100
  )
  
  output<-list(glm=twovar_fit,
               results=results,
               blood=c(TN.bld=TN.bld,TP.bld=TP.bld,
                       FN.bld=FN.bld,FP.bld=FP.bld,
                       Se.bld=Se.bld,Sp.bld=Sp.bld,
                       PPV.bld=PPV.bld,NPV.bld=NPV.bld,
                       FPR.bld=FPR.bld,FNR.bld=FNR.bld),
               val_data=val_data)
  
}


####################
#import data, create variables, subset by age
################### import/merge ####
#starting with full dataset - moderate/severe diarrhea
gems1_orig <- read.csv("gems1.csv", header=T)

cases <- gems1_orig %>% filter(type=="Case")

cases=cases %>% select(site,	type, caseid, f3_gender,	f3_drh_turgor	,	f3_drh_iv	,	f3_drh_hosp,
                            f4a_relationship,	f4a_dad_live	,					
                            f4a_prim_schl,	f4a_ppl_house	,	f4a_yng_children	,			
                            f4a_slp_rooms,	f4a_floor	,	f4a_house_elec	,  			
                            f4a_house_bike,	f4a_house_phone	,	f4a_house_tele	,    			
                            f4a_house_car	,	f4a_house_cart	,	f4a_house_scoot	,    			
                            f4a_house_fridge	,	f4a_house_agland	,	f4a_house_radio	,   			
                            f4a_house_boat	,	f4a_house_none	,	f4a_fuel_elec	,   			
                            f4a_fuel_biogas	,	f4a_fuel_grass	,	f4a_fuel_propane	,   			
                            f4a_fuel_coal	,	f4a_fuel_dung	,	f4a_fuel_natgas	, 			
                            f4a_fuel_charcoal	,	f4a_fuel_crop	,	f4a_fuel_kero	,   			
                            f4a_fuel_wood	,	f4a_fuel_other	,	f4a_ani_goat	,     			
                            f4a_ani_sheep	,	f4a_ani_dog	,	f4a_ani_cat	,      			
                            f4a_ani_cow	,	f4a_ani_rodents	,	f4a_ani_fowl	,       			
                            f4a_ani_other	,	f4a_ani_no	,	f4a_water_house	,    			
                            f4a_water_covwell	,	f4a_water_yard	,	f4a_water_covpwell	, 			
                            f4a_water_pubtap	,	f4a_water_prospring	,	f4a_water_well	,			
                            f4a_water_unspring	,	f4a_water_pubwell	,	f4a_water_river	,    			
                            f4a_water_pond	,	f4a_water_deepwell	,	f4a_water_rain	,   			
                            f4a_water_shallwell	,	f4a_water_bought	,	f4a_water_othr	,    			
                            f4a_water_bore	,	f4a_ms_water	,	f4a_fetch_water	,    			
                            f4a_trip_day	,	f4a_trip_week	,	f4a_water_avail	,   			
                            f4a_store_water	,	f4a_trt_water	,	f4a_trt_method	,   			
                            f4a_notrt_water	,	f4a_disp_feces	,    					
                            f4a_fac_waste	,	f4a_share_fac	,	f4a_wash_eat	,    			
                            f4a_wash_cook	,	f4a_wash_nurse	,	f4a_wash_def	,      			
                            f4a_wash_animal	,	f4a_wash_child	,	f4a_wash_othr	,      			
                            f4a_wash_use	,	f4a_breastfed	,	f4a_drh_days	,     			
                            f4a_max_stools	,	f4a_drh_blood	,	f4a_drh_vomit	,      			
                            f4a_drh_thirst	,	f4a_drh_lessdrink	,    					
                            f4a_drh_bellypain	,	f4a_drh_restless	,   					
                            f4a_drh_lethrgy	,	f4a_drh_consc	,	f4a_drh_strain	,  			
                            f4a_drh_prolapse	,	f4a_drh_cough	,    					
                            f4a_drh_conv	,	f4a_cur_thirsty	,	f4a_cur_skin	,    			
                            f4a_cur_restless	,	f4a_cur_drymouth	,   					
                            f4a_cur_fastbreath	,	f4a_hometrt_ors	,	f4a_hometrt_maize	,  			
                            f4a_hometrt_milk	,	f4a_hometrt_herb	,	f4a_hometrt_zinc	, 			
                            f4a_hometrt_none	,	f4a_hometrt_othrliq	,	f4a_hometrt_ab	,  			
                            f4a_hometrt_othr1	,	f4a_hometrt_othr2	,	f4a_offr_drink	,    			
                            f4a_seek_outside	,	f4a_seek_pharm	,	f4a_seek_friend	,    			
                            f4a_seek_healer	,	f4a_seek_doc	,	f4a_seek_privdoc	,   			
                            f4a_seek_remdy	,	f4a_seek_other	,	f4b_haz	,  			
                            f4b_muac	,	f4b_temp	,	f4b_resp	,           			
                            f4b_chest_indrw	,	f4b_eyes	,	f4b_mouth	,          			
                            f4b_skin	,	f4b_mental	,	f4b_rectal	,         			
                            f4b_bipedal	,	f4b_abn_hair	,	f4b_under_nutr	,     			
                            f4b_skin_flaky	,	f4b_observe_stool	,	f4b_nature_stool	,   			
                            f4b_recommend	,	f4b_volume	,  					
                            f4b_admit	,	f9_memory_aid	,	wealth_index	,       			
                            wiq	,	base_age, f4b_date, f4b_outcome, agegroup,
                            f5_date ,
                            f11_consistency, f11_blood, f11_pus, f11_mucus
)


#drop the extra 1 at the end of each caseid
cases$caseid <- as.numeric(substr(cases$caseid,1,9))
AFes <- read.csv("GEMS with AFes.csv", header=T)
#is no data here for controls, so just drop those
AFes <- AFes %>% filter(case.control==1)
cases <- cases %>% inner_join(AFes, by=c("caseid"="Case.ID")) #only want cases that have AFes
#gems1 is caseid, originally 10 digits long
#AFes is Case.ID, 9 digits long

################### define variables ####
#create binary etiology variables
#has to be a more elegant way of doing this
cases$astro <- cases$astrovirus_afe
cases$NoV <- cases$norovirus_gii_afe
cases$rota <- cases$rotavirus_afe
cases$sapo <- cases$sapovirus_afe
cases$adeno <- cases$adenovirus_40_41_afe
cases$aero <- cases$aeromonas_afe
cases$campy <- cases$c_jejuni_coli_afe
cases$crypto <- cases$cryptosporidium_afe
cases$ehisto <- cases$e_histolytica_afe
cases$cyclo <- cases$cyclospora_afe
cases$hpylori <- cases$h_pylori_afe
cases$isospora <- cases$isospora_afe
cases$salm <- cases$salmonella_afe
cases$shigella <- cases$shigella_eiec_afe
cases$cholera <- cases$v_cholerae_afe
cases$EAEC <- cases$EAEC_afe
cases$ST_ETEC <- cases$ST_ETEC_afe
cases$LT_ETEC <- cases$LT_ETEC_afe
cases$TEPEC <- cases$TEPEC_afe
cases$STEC <- cases$STEC_afe

func_bi <- function(x, na.rm=F) (x=ifelse(x>=0.5,1,0)) #function if AFes>=0.5, 1(present), else 0(absent)
cases <- cases %>% mutate_at(c("astro","NoV","rota","sapo","adeno",
                               "aero","campy","crypto","ehisto","cyclo",
                               "hpylori","isospora","salm","shigella","cholera",
                               "EAEC","ST_ETEC","LT_ETEC","TEPEC","STEC"), 
                             func_bi)
cases <- cases %>% rename_with(~paste0(., "_afe0.5"), astro:STEC) %>%
  mutate(shigella_afe0.3 = ifelse(shigella_eiec_afe>=0.3,1,0),
                        shigella_afe0.7 = ifelse(shigella_eiec_afe>=0.7,1,0),
                        cholera_afe0.3 = ifelse(v_cholerae_afe>=0.3,1,0),
                        cholera_afe0.7 = ifelse(v_cholerae_afe>=0.7,1,0))

cases$f4b_date_date <- as.Date(as.character(cases$f4b_date))
cases$f5_date_date <- as.Date(as.character(cases$f5_date))
cases$fup_days <- as.numeric(cases$f5_date_date - cases$f4b_date_date)

cases=cases %>% mutate(any_breast_fed=factor(case_when((f4a_breastfed==1|f4a_breastfed==2)~1,TRUE~0))) #SMA dichotomizing breastfeeding
cases=cases %>% mutate(any_breast_fed2=factor(case_when((f4a_breastfed==0|f4a_breastfed==1)~0,TRUE~1)))
cases=cases %>% mutate(cont=case_when(site %in% c(1,2,3,4) ~ 1,
                                      TRUE ~ 2)) #SMA creating "cont"inent variable
cases$index=1:dim(cases)[1] #SMA creating an ID for each observation

 
#>>> recode these dk as no
func_dk <- function(x, na.rm=F) (x=ifelse(x==9,0,x)) #function if==9 (unknown), then recodes as 0(no)
cases <- cases %>% mutate_at(c("f4a_drh_blood","f4a_drh_vomit","f4a_drh_thirst",
                               "f4a_drh_lessdrink","f4a_drh_bellypain",
                               "f4a_drh_restless","f4a_drh_lethrgy",
                               "f4a_drh_consc","f4a_drh_strain","f4a_drh_prolapse",
                               "f4a_drh_cough","f4a_drh_conv"), 
                             #f4a_drh_undrink, f4a_drh_fever, f4a_drh_breath in codebook but in data
                             func_dk)

#relabeling site
cases <- cases %>% mutate(site=case_when(site==1~"The Gambia",
                                        site==2~"Mali",
                                        site==3~"Mozambique",
                                        site==4~"Kenya",
                                        site==5~"India",
                                        site==6~"Bangladesh",
                                        site==7~"Pakistan",
                                        TRUE~"no")) 

cases$f4a_drh_lethrgy <- as.numeric(cases$f4a_drh_lethrgy)
cases <- cases %>% mutate(f4a_drh_lethrgy_miss=case_when((f4a_drh_lethrgy==9) ~ as.numeric(NA), TRUE~f4a_drh_lethrgy))
cases$f4a_drh_restless <- as.numeric(cases$f4a_drh_restless)
cases <- cases %>% mutate(f4a_drh_restless_miss=case_when((f4a_drh_restless==9) ~ as.numeric(NA), TRUE~f4a_drh_restless))

# #f4b_outcome: 1=resolved; 2=improved; 3=no better; 4=worse; 5=died at hosp; 6=unknown/LTF
# #f5_child_health: 1=appears healthy; 2=improved but not back to normal; 3=no better; 4=worse; 5=died
#create variable where NA missing have been set to value 9
cases=cases %>% mutate(f4b_outcome_miss=replace_na(f4b_outcome,9))

#very few 5 so put in other category
cases$f4a_disp_feces <- as.numeric(cases$f4a_disp_feces)
cases=cases %>% mutate(f4a_disp_feces=(case_when(f4a_disp_feces==5 ~ 6, TRUE~f4a_disp_feces)))

#combining categories. not sure what 0 is, isn't in dictionary
table(cases$f4a_trt_method)
cases$f4a_trt_method <- as.numeric(cases$f4a_trt_method)
cases=cases %>% mutate(f4a_trt_method=(case_when(f4a_trt_method==1 | f4a_trt_method==6 | f4a_trt_method==7 ~ 7, TRUE~f4a_trt_method)))
table(cases$test)
# 0    2    3    4    5    7 
# 6740  656  751 1197   73   22 

cases <- cases %>% mutate(month=month(f4b_date_date))

#combine to fewer categories: f4a_ms_water
table(cases$f4a_ms_water)
#   1    2    3    4    5    6    7    8   10   11   12   13   15   16   17   18 
# 668  933 3481   72  272  109  972  703  127   58   56  374  529  658  336   91 
cases=cases %>% mutate(f4a_ms_water=(case_when((f4a_ms_water==6 | f4a_ms_water==13 | f4a_ms_water==14) ~ 0, #surface 
                                                         (f4a_ms_water==4 | f4a_ms_water==5 | f4a_ms_water==12 | f4a_ms_water==16) ~ 1, #unimproved
                                                         (f4a_ms_water==3 | f4a_ms_water==7 | f4a_ms_water==8 | f4a_ms_water==9 |
                                                            f4a_ms_water==10 | f4a_ms_water==11 | f4a_ms_water==15 | f4a_ms_water==17) ~ 2, #other improved
                                                         (f4a_ms_water==1 | f4a_ms_water==2) ~ 3, #piped
                                                         TRUE~4))) #other
# table(cases$test)
# 0    1    2    3    4 
# 483 1058 6206 1601   91 
# #use JMP drinking water services ladder
# #surface (subset of unimproved)
# 6-pond/lake
# 13-river/stream
# 14-dam/earth pan
# #unimproved
# 4-open well in house/yard
# 5-open public well
# 12-unprotected spring
# 16-bought
# #other improved
# 3-public tap
# 7-deep tube well
# 8-shallow tube well
# 9-covered well in house/yard
# 10-covered public well
# 11-protected spring
# 15-rainwater
# 17-bore hole
# #safely managed/ piped into dwelling/plot/yard
# 1-piped into house
# 2-piped into yard
# #other
# 18-other

cases=cases %>% mutate(f4a_floor=(case_when((f4a_floor==1 | f4a_floor==2 |
                                                         f4a_floor==3 | f4a_floor==4 | f4a_floor==10) ~ 0, #natural, rudimentary, other floor
                                                      TRUE~1))) #finished floor
# table(cases$test)
# # 0    1 
# # 1933 3371

table(cases$f4a_prim_schl)
# 1    2    3    4    5    6    7 
# 1742 1081 1594  319   73  475   19 
cases$f4a_prim_schl <- as.numeric(cases$f4a_prim_schl)
cases=cases %>% mutate(f4a_prim_schl=(case_when((f4a_prim_schl==7) ~ 1,
                                                TRUE~f4a_prim_schl)))

#creating a combined category for non-father male relation OR non relation (4,6,8,9)
cases$f4a_relationship <- as.numeric(cases$f4a_relationship)
cases=cases %>% mutate(f4a_relationship=(case_when((f4a_relationship==4 | f4a_relationship==6 | f4a_relationship==8 | f4a_relationship==9) ~ 9, #non-father male relation OR 
                                                             TRUE~f4a_relationship))) #what was originally


cases$f4b_date_date <- as.Date(as.character(cases$f4b_date))
cases$f5_date_date <- as.Date(as.character(cases$f5_date))
cases$fup_days <- as.numeric(cases$f5_date_date - cases$f4b_date_date)
cases <- cases %>% mutate(month=month(f4b_date_date))
cases <- cases %>% mutate(season=(case_when((month==10 | month==11 | month==12 | month==1 | month==2 | month==3) ~ 0,
                                            TRUE~1))) #0=cold months, 1=warm months

#convert these to factors
vars <- c("f4a_ms_water","f4a_fac_waste","f4a_dad_live","f4b_recommend",
          "f4a_relationship","f4a_prim_schl","f4a_floor","f4a_disp_feces",
          "f4a_wash_use","f4a_water_avail","f4a_trt_method","f4a_drh_blood",
          "f4a_drh_vomit","f4a_drh_thirst","f4a_drh_lessdrink","f4a_drh_bellypain",
          "f4a_drh_restless_miss","f4a_drh_lethrgy_miss","f4a_drh_consc","f4a_drh_strain",
          "f4a_drh_prolapse","f4a_drh_cough","f4a_drh_conv","f4a_cur_thirsty",
          "f4a_cur_skin","f4a_cur_restless","f4a_cur_drymouth","f4a_cur_fastbreath",
          "f4b_nature_stool","month","site","season"
)
cases[vars] <- lapply(cases[vars], factor)

#collected as categorical, but ordinal so leaving as numeric for now: 
#f4a_prim_schl, f4a_offr_drink, f4a_max_stools, f4a_breastfed, f4b_mouth, f4b_skin, f4b_mental


################### inclusion/exclusion ####
#already limited to only those w/ known etiology in import code 
cases <- cases %>% filter(f4b_haz>=-7.0 & f4b_haz <=7.0)

################### drop missing since can't have missing in RF, define "names" variables interested in ####
cases_full <- cases
cases <- cases %>% filter(!is.na(f4a_dad_live)&!is.na(f4a_prim_schl)&!is.na(f4a_slp_rooms)&!is.na(f4a_water_avail)&
                            !is.na(f4a_disp_feces)&!is.na(f4a_share_fac)&!is.na(f4a_wash_use)&!is.na(f4a_drh_days)&
                            !is.na(f4a_drh_thirst)&!is.na(f4a_drh_restless)&!is.na(f4a_drh_lethrgy)&!is.na(f4a_drh_conv)&
                            !is.na(f4a_cur_skin)&!is.na(f4a_cur_restless)&!is.na(f4a_cur_drymouth)&!is.na(f4a_cur_fastbreath)&
                            !is.na(f4a_offr_drink)&!is.na(f4b_temp)&!is.na(f4b_chest_indrw)&!is.na(f4b_mouth)&
                            !is.na(f4b_skin)&!is.na(f4b_under_nutr)&!is.na(f4b_nature_stool)&!is.na(f3_drh_iv)&
                            !is.na(f4a_ppl_house)&!is.na(f4a_store_water)&!is.na(f4b_haz)&!is.na(f4b_resp)&
                            !is.na(f4a_drh_cough)&!is.na(f4b_abn_hair)
                            &!is.na(f11_consistency)&!is.na(f11_blood)&!is.na(f11_pus)&!is.na(f11_mucus)
)

# select variables we're interested in 
names <- c("f3_gender","f3_drh_turgor","f3_drh_iv","f3_drh_hosp", 
           "site",
           "f4a_relationship","f4a_dad_live",
           "f4a_prim_schl","f4a_ppl_house","f4a_yng_children",
           "f4a_slp_rooms","f4a_floor","f4a_house_elec",  
           "f4a_house_bike","f4a_house_phone","f4a_house_tele",    
           "f4a_house_car","f4a_house_cart","f4a_house_scoot",    
           "f4a_house_fridge","f4a_house_agland","f4a_house_radio",   
           "f4a_house_boat","f4a_house_none","f4a_fuel_elec",   
           "f4a_fuel_biogas","f4a_fuel_grass","f4a_fuel_propane",   
           "f4a_fuel_coal","f4a_fuel_dung","f4a_fuel_natgas", 
           "f4a_fuel_charcoal","f4a_fuel_crop","f4a_fuel_kero",   
           "f4a_fuel_wood","f4a_fuel_other","f4a_ani_goat",     
           "f4a_ani_sheep","f4a_ani_dog","f4a_ani_cat",      
           "f4a_ani_cow","f4a_ani_rodents","f4a_ani_fowl",       
           "f4a_ani_other","f4a_ani_no","f4a_water_house",    
           "f4a_water_covwell","f4a_water_yard","f4a_water_covpwell", 
           "f4a_water_pubtap","f4a_water_prospring","f4a_water_well",
           "f4a_water_unspring","f4a_water_pubwell","f4a_water_river",    
           "f4a_water_pond","f4a_water_deepwell","f4a_water_rain",   
           "f4a_water_shallwell","f4a_water_bought","f4a_water_othr",    
           "f4a_water_bore",
           "f4a_ms_water",
           "f4a_water_avail",   
           "f4a_store_water","f4a_trt_water","f4a_trt_method",   
           "f4a_disp_feces",    
           "f4a_fac_waste","f4a_share_fac","f4a_wash_eat",    
           "f4a_wash_cook","f4a_wash_nurse","f4a_wash_def",      
           "f4a_wash_animal","f4a_wash_child","f4a_wash_othr",      
           "f4a_wash_use","f4a_breastfed","f4a_drh_days",     
           "f4a_max_stools","f4a_drh_blood","f4a_drh_vomit",      
           "f4a_drh_thirst","f4a_drh_lessdrink",    
           "f4a_drh_bellypain","f4a_drh_restless",   
           "f4a_drh_lethrgy_miss","f4a_drh_consc","f4a_drh_strain",  
           "f4a_drh_prolapse","f4a_drh_cough",    
           "f4a_drh_conv","f4a_cur_thirsty","f4a_cur_skin",    
           "f4a_cur_restless","f4a_cur_drymouth",   
           "f4a_cur_fastbreath","f4a_hometrt_ors","f4a_hometrt_maize",  
           "f4a_hometrt_milk","f4a_hometrt_herb","f4a_hometrt_zinc", 
           "f4a_hometrt_none","f4a_hometrt_othrliq","f4a_hometrt_ab",  
           "f4a_hometrt_othr1","f4a_hometrt_othr2","f4a_offr_drink",    
           "f4a_seek_outside","f4a_seek_pharm","f4a_seek_friend",    
           "f4a_seek_healer","f4a_seek_doc","f4a_seek_privdoc",   
           "f4a_seek_remdy","f4a_seek_other","f4b_haz",  
           "f4b_temp","f4b_resp",           
           "f4b_chest_indrw","f4b_eyes","f4b_mouth",          
           "f4b_skin","f4b_mental","f4b_rectal",         
           "f4b_bipedal","f4b_abn_hair","f4b_under_nutr",     
           "f4b_skin_flaky",
           "f4b_recommend",  
           "f4b_admit",
           "base_age"
)

#save(cases, file = "")

################### into age groups, continents, cholera site, bloody ####
#into age groups
table(cases$agegroup)
summary(cases$agegroup)
#agegroup 1=0-11mo, 2=12-23mo, 3=24-50mo, 4=0-23mo
cases_age1 <- cases %>% filter(agegroup == 1)
cases_age2 <- cases %>% filter(agegroup == 2)
cases_age3 <- cases %>% filter(agegroup == 3)
cases_age4 <- cases %>% filter(agegroup == 1 | agegroup==2)

Gambia <- cases %>% filter(site == "The Gambia")
Mali <- cases %>% filter(site == "Mali")
Moz <- cases %>% filter(site == "Mozambique")
Kenya <- cases %>% filter(site == "Kenya")
India <- cases %>% filter(site == "India")
Bdesh <- cases %>% filter(site == "Bangladesh")
Pak <- cases %>% filter(site == "Pakistan")

Afr <- cases %>% filter(site=="The Gambia" | site=="Mali" | site=="Mozambique" | site=="Kenya")
Asia <- cases %>% filter(site=="India" | site=="Pakistan" | site=="Bangladesh")
Asia.noBdesh <- cases %>% filter(site=="India" | site=="Pakistan")

cases.cholera <- cases %>% filter(site=="Bangladesh" | site=="India" | site=="Mozambique" | site=="Pakistan")
cases.noncholera <- cases %>% filter(site=="Kenya" | site=="Mali" | site=="The Gambia")
cases.noBdesh <- cases %>% filter(site!="Bangladesh")

table(cases$site,cases$cholera)
#              0   1
# Bangladesh 845  11
# India      790  41
# Kenya      623   1
# Mali       828   4
# Mozambique 434  13
# Pakistan   692  70
# The Gambia 662   1

cases.bloody <- cases %>% filter(f4a_drh_blood == 1)
cases.nonbloody <- cases %>% filter(f4a_drh_blood == 0)

################### descriptive for pub ####
dim(cases_full)
# [1] 5287  209

table(cases_full$site)
table(cases_full$shigella_afe0.3,cases_full$site)
round((table(cases_full$shigella_afe0.3,cases_full$site)[2,])/(table(cases_full$site))*100,2)
table(cases_full$shigella_afe0.5,cases_full$site)
round((table(cases_full$shigella_afe0.5,cases_full$site)[2,])/(table(cases_full$site))*100,2)
table(cases_full$shigella_afe0.7,cases_full$site)
round((table(cases_full$shigella_afe0.7,cases_full$site)[2,])/(table(cases_full$site))*100,2)

dim(cases_full)
table(cases_full$shigella_afe0.3)
round((table(cases_full$shigella_afe0.3)[2])/(dim(cases_full)[1])*100,2)
table(cases_full$shigella_afe0.5)
round((table(cases_full$shigella_afe0.5)[2])/(dim(cases_full)[1])*100,2)
table(cases_full$shigella_afe0.7)
round((table(cases_full$shigella_afe0.7)[2])/(dim(cases_full)[1])*100,2)


table(cases_full$agegroup)
table(cases_full$shigella_afe0.3,cases_full$agegroup)
round((table(cases_full$shigella_afe0.3,cases_full$agegroup)[2,])/(table(cases_full$agegroup))*100,2)
table(cases_full$shigella_afe0.5,cases_full$agegroup)
round((table(cases_full$shigella_afe0.5,cases_full$agegroup)[2,])/(table(cases_full$agegroup))*100,2)
table(cases_full$shigella_afe0.7,cases_full$agegroup)
round((table(cases_full$shigella_afe0.7,cases_full$agegroup)[2,])/(table(cases_full$agegroup))*100,2)

#overlapping histograms comparing distributions of top predictive variables
#Shigella
p1 <- hist(cases[which(cases$shigella_afe0.5==0),]$base_age,freq=F)
p2 <- hist(cases[which(cases$shigella_afe0.5==1),]$base_age,freq=F)
p3 <- hist(data_MALED[which(data_MALED$shigella_afe0.5==0),]$base_age,freq=F)
p4 <- hist(data_MALED[which(data_MALED$shigella_afe0.5==1),]$base_age,freq=F)
#jpeg("overlap_hist_age_shigella.jpg",width=480,height=480,quality=400)
par(mfrow=c(2,1))
plot( p1, col=rgb(1,0,0,1/4), xlim=c(0,60),ylim=c(0,0.1),freq=F,xlab="Age (months)",main="Shigella etiology by age in GEMS")
plot( p2, col=rgb(0,0,1,1/4),  xlim=c(0,60),ylim=c(0,0.1),freq=F, add=T)
legend('topright',c('no Shigella','Shigella'),fill = c(rgb(1,0,0,1/4),rgb(0,0,1,1/4)), bty = 'n')
plot( p3, col=rgb(1,0,0,1/4), xlim=c(0,60),ylim=c(0,0.1),freq=F,xlab="Age (months)",main="Shigella etiology by age in MAL-ED")
plot( p4, col=rgb(0,0,1,1/4),  xlim=c(0,60),ylim=c(0,0.1),freq=F, add=T)
legend('topright',c('no Shigella','Shigella'),fill = c(rgb(1,0,0,1/4),rgb(0,0,1,1/4)), bty = 'n')
dev.off()
#non-overlaping version for R2
#jpeg("overlap_hist_age_shigella_unoverlap.jpg",width=480,height=480,quality=400)
par(mfcol=c(2,2))
plot( p1, col=rgb(1,0,0,1/4), xlim=c(0,60),ylim=c(0,1200),freq=T,xlab="Age (months)",main="Non-Shigella etiology by age in GEMS")
plot( p2, col=rgb(0,0,1,1/4),  xlim=c(0,60),ylim=c(0,400),freq=T, add=F,xlab="Age (months)",main="Shigella etiology by age in GEMS")
#legend('topright',c('no Shigella','Shigella'),fill = c(rgb(1,0,0,1/4),rgb(0,0,1,1/4)), bty = 'n')
plot( p3, col=rgb(1,0,0,1/4), xlim=c(0,60),ylim=c(0,175),freq=T,xlab="Age (months)",main="Non-Shigella etiology by age in MAL-ED")
plot( p4, col=rgb(0,0,1,1/4),  xlim=c(0,60),ylim=c(0,50),freq=T, add=F,xlab="Age (months)",main="Shigella etiology by age in MAL-ED")
#legend('topright',c('no Shigella','Shigella'),fill = c(rgb(1,0,0,1/4),rgb(0,0,1,1/4)), bty = 'n')
dev.off()


table(cases$shigella_afe0.5)
table(cases$shigella_afe0.5,cases$f4a_drh_blood)

p1 <- hist(cases[which(cases$f4a_drh_blood==0),]$base_age,freq=F)
p2 <- hist(cases[which(cases$f4a_drh_blood==1),]$base_age,freq=F)
p3 <- hist(data_MALED[which(data_MALED$f4a_drh_blood==0),]$base_age,freq=F)
p4 <- hist(data_MALED[which(data_MALED$f4a_drh_blood==1),]$base_age,freq=F)
#jpeg("overlap_hist_ageBD.jpg",width=480,height=480,quality=400)
par(mfrow=c(2,1))
plot( p1, col=rgb(1,0,0,1/4), xlim=c(0,60),ylim=c(0,0.1),freq=F,xlab="Age (months)",main="Dystentery (bloody diarrhea) by age in GEMS")
plot( p2, col=rgb(0,0,1,1/4),  xlim=c(0,60),ylim=c(0,0.1),freq=F, add=T)
legend('topright',c('no dysentery','dysentery'),fill = c(rgb(1,0,0,1/4),rgb(0,0,1,1/4)), bty = 'n')
plot( p3, col=rgb(1,0,0,1/4), xlim=c(0,60),ylim=c(0,0.1),freq=F,xlab="Age (months)",main="Dystentery (bloody diarrhea) by age in MAL-ED")
plot( p4, col=rgb(0,0,1,1/4),  xlim=c(0,60),ylim=c(0,0.1),freq=F, add=T)
legend('topright',c('no dysentery','dysentery'),fill = c(rgb(1,0,0,1/4),rgb(0,0,1,1/4)), bty = 'n')
dev.off()
#non-overlaping version for R2
#jpeg("overlap_hist_ageBD_unoverlap.jpg",width=480,height=480,quality=400)
par(mfcol=c(2,2))
plot( p1, col=rgb(1,0,0,1/4), xlim=c(0,60),ylim=c(0,1000),freq=T,xlab="Age (months)",main="No Dystentery (bloody diarrhea) \nby age in GEMS")
plot( p2, col=rgb(0,0,1,1/4),  xlim=c(0,60),ylim=c(0,250),freq=T, add=F,xlab="Age (months)",main="Dystentery (bloody diarrhea) \nby age in GEMS")
#legend('topright',c('no Shigella','Shigella'),fill = c(rgb(1,0,0,1/4),rgb(0,0,1,1/4)), bty = 'n')
plot( p3, col=rgb(1,0,0,1/4), xlim=c(0,60),ylim=c(0,175),freq=T,xlab="Age (months)",main="No Dystentery (bloody diarrhea) \nby age in MAL-ED")
plot( p4, col=rgb(0,0,1,1/4),  xlim=c(0,60),ylim=c(0,20),freq=T, add=F,xlab="Age (months)",main="Dystentery (bloody diarrhea) \nby age in MAL-ED")
#legend('topright',c('no Shigella','Shigella'),fill = c(rgb(1,0,0,1/4),rgb(0,0,1,1/4)), bty = 'n')
dev.off()


set.seed(11)
p1 <- hist(rnorm(100, mean=0.3, sd=0.1), freq=F)
p2 <- hist(rnorm(100, mean=0.6, sd=0.1), freq=F)
#jpeg("",width=600,height=480,quality=100)
plot( p1, col=rgb(1,0,0,1/4), freq=F, xlim=c(0.1,0.9),ylim=c(0,8),xlab="predicted probability of treatable bacteria",
      main="Schematic of prediction rule + \npoint-of-care diagnostic testing regimen")
plot( p2, col=rgb(0,0,1,1/4), freq=F, xlim=c(0.1,0.9),ylim=c(0,8),add=T)
legend('topright',c('Non-Treatable','Treatable'),fill = c(rgb(1,0,0,1/4),rgb(0,0,1,1/4)), bty = 'n')
dev.off()


####################
#shigella AFe>=0.5 only MUAC
################### var screening, AUC ####
names<-append(x=names, values=c("f4b_muac")) #add MUAC to list
names <- names[!names %in% c("f4b_haz")] #take HAZ off list

shigella0.5MUAC <- CPR.funct(data=cases,outcome="shigella_afe0.5",iter=10,nvars_opts=c(1:10,15,20,30,40,50))
#save(shigella0.5MUAC, file = "shigella0.5MUAC_10iter.Rdata")
shigella0.5MUAC[["df_imps"]]
shigella0.5MUAC[["AUC_df"]]
shigella0.5MUAC[["calib"]]

# names    var_red
# 1               base_age 90.2725970
# 2          f4a_drh_blood 80.6839740
# 3               f4b_muac 33.0383724
# 4               f4b_resp 31.1394297
# 5               f4b_temp 29.3389628
# 6               f4b_eyes 27.4325513
# 7          f4a_ppl_house 23.2759654
# 8                   site 20.6930675
# 9           f4a_drh_days 15.8631886
# 10         f4a_slp_rooms 14.5896492
# 11         f4a_share_fac 13.9315800
# 12        f4a_drh_strain 13.1195783
# 13      f4a_yng_children 12.4728347
# 14         f4a_drh_vomit 11.7743292
# 15         f4a_prim_schl 11.7210086
# 16         f4a_breastfed 11.0273569
# 17             f4b_mouth 10.1072213
# 18        f4a_fuel_grass 10.0887164
# 19        f4a_max_stools  9.2834566
# 20        f4a_offr_drink  9.0130183
# 21         f4b_recommend  8.8758314
# 22       f4a_store_water  8.6716567
# 23         f4a_fuel_crop  8.0536639
# 24     f4a_drh_bellypain  7.8530248
# 25         f4a_fac_waste  7.8034451
# 26       f4a_water_avail  7.2097737
# 27        f4a_disp_feces  7.1626022
# 28         f4a_drh_cough  7.1399925
# 29       f4a_cur_thirsty  7.0889230
# 30        f4a_trt_method  6.6589941
# 31        f4a_drh_thirst  6.6055259
# 32          f4a_ms_water  6.5979236
# 33          f4a_dad_live  6.5457569
# 34            f4b_mental  5.7503421
# 35      f4a_hometrt_none  5.3114327
# 36      f4a_cur_drymouth  5.1900804
# 37       f4a_hometrt_ors  5.1791854
# 38             f3_gender  5.1729355
# 39          f4a_wash_use  5.0463397
# 40         f4a_fuel_dung  4.9659207
# 41         f4a_wash_cook  4.9555953
# 42        f4a_wash_child  4.9018448
# 43        f4a_wash_nurse  4.8976888
# 44     f4a_fuel_charcoal  4.8203296
# 45  f4a_drh_lethrgy_miss  4.7924440
# 46        f4a_house_tele  4.6665433
# 47        f4a_house_bike  4.6545415
# 48          f4a_wash_def  4.6248748
# 49              f4b_skin  4.6216457
# 50      f4a_seek_outside  4.5740790
# 51          f4a_ani_fowl  4.5563203
# 52           f4a_ani_cat  4.4501613
# 53             f4a_floor  4.4239691
# 54      f4a_cur_restless  4.4083922
# 55       f4a_house_phone  4.3925041
# 56      f4a_relationship  4.3848435
# 57          f4a_wash_eat  4.3845270
# 58      f4a_drh_restless  4.3391110
# 59          f4a_cur_skin  4.2401432
# 60         f4a_fuel_wood  4.2270292
# 61        f4a_seek_pharm  4.1682246
# 62       f4a_house_radio  4.1352316
# 63             f4b_admit  4.0109713
# 64    f4a_water_deepwell  3.9997955
# 65        f4a_house_elec  3.9688527
# 66       f4a_house_scoot  3.9206740
# 67           f3_drh_hosp  3.8596025
# 68           f4a_ani_dog  3.8562681
# 69       f4a_ani_rodents  3.8321439
# 70          f4a_ani_goat  3.7471518
# 71      f4a_house_fridge  3.7433486
# 72      f4a_water_pubtap  3.7124833
# 73      f4a_house_agland  3.7028233
# 74         f4a_trt_water  3.6875504
# 75           f4a_ani_cow  3.4868889
# 76     f4a_drh_lessdrink  3.4613735
# 77     f4a_hometrt_maize  3.4525924
# 78    f4a_cur_fastbreath  3.4315573
# 79             f3_drh_iv  3.4215893
# 80        f4b_under_nutr  3.3980023
# 81         f3_drh_turgor  3.3874569
# 82     f4a_hometrt_othr1  3.3219818
# 83        f4a_hometrt_ab  3.0179101
# 84       f4a_fuel_natgas  2.9149662
# 85            f4a_ani_no  2.8033651
# 86   f4a_water_shallwell  2.7557537
# 87         f4a_ani_sheep  2.7414064
# 88       f4a_wash_animal  2.6265525
# 89         f4a_wash_othr  2.5515715
# 90      f4a_hometrt_herb  2.3797147
# 91        f4a_water_yard  2.3715187
# 92         f4a_house_car  2.3605984
# 93       f4a_water_house  2.2867022
# 94          f4a_drh_conv  2.2622012
# 95         f4a_ani_other  2.2416818
# 96        f4a_house_cart  2.2270992
# 97      f4a_water_bought  2.0725486
# 98       f4a_seek_healer  2.0054332
# 99         f4a_fuel_kero  1.7856934
# 100    f4a_water_pubwell  1.7418388
# 101       f4a_house_none  1.5351469
# 102     f4a_seek_privdoc  1.5154535
# 103     f4a_hometrt_zinc  1.5093559
# 104     f4a_fuel_propane  1.4981548
# 105       f4a_water_bore  1.4936243
# 106       f4a_seek_remdy  1.4655070
# 107         f4b_abn_hair  1.4357273
# 108       f4a_water_well  1.3691646
# 109         f4a_seek_doc  1.3007742
# 110       f4a_fuel_other  1.2971766
# 111     f4a_hometrt_milk  1.2696419
# 112        f4a_fuel_coal  1.2477525
# 113       f4a_water_rain  1.1984798
# 114    f4a_hometrt_othr2  1.1539242
# 115      f4b_chest_indrw  1.1305168
# 116        f4a_fuel_elec  1.1144067
# 117        f4a_drh_consc  1.1135748
# 118       f4b_skin_flaky  1.1050455
# 119    f4a_water_covwell  1.0208597
# 120     f4a_drh_prolapse  1.0074668
# 121   f4a_water_covpwell  0.8626334
# 122       f4a_water_othr  0.7526533
# 123      f4a_water_river  0.7389567
# 124      f4a_fuel_biogas  0.7384299
# 125       f4a_seek_other  0.7059038
# 126  f4a_hometrt_othrliq  0.6791278
# 127       f4a_house_boat  0.6492218
# 128       f4a_water_pond  0.5181040
# 129          f4b_bipedal  0.4752391
# 130           f4b_rectal  0.4374171
# 131      f4a_seek_friend  0.3502214
# 132   f4a_water_unspring  0.2974277
# 133  f4a_water_prospring  0.1550710

# AUC          SE     lower     upper level Model nvar
# 1  0.6873757 0.001743926 0.6839577 0.6907938  0.95    LR    1
# 2  0.8020408 0.001619849 0.7988660 0.8052157  0.95    LR    2
# 3  0.8005508 0.001601236 0.7974124 0.8036891  0.95    LR    3
# 4  0.7994617 0.001610332 0.7963055 0.8026179  0.95    LR    4
# 5  0.7933739 0.001675916 0.7900891 0.7966586  0.95    LR    5
# 6  0.7954947 0.001668657 0.7922242 0.7987652  0.95    LR    6
# 7  0.7952680 0.001673445 0.7919881 0.7985479  0.95    LR    7
# 8  0.7989847 0.001667338 0.7957168 0.8022526  0.95    LR    8
# 9  0.7989317 0.001668425 0.7956616 0.8022017  0.95    LR    9
# 10 0.7992743 0.001667419 0.7960062 0.8025424  0.95    LR   10
# 11 0.7994886 0.001666742 0.7962218 0.8027553  0.95    LR   15
# 12 0.8044519 0.001640187 0.8012371 0.8076666  0.95    LR   20
# 13 0.8026810 0.001655248 0.7994368 0.8059252  0.95    LR   30
# 14 0.7998712 0.001672516 0.7965931 0.8031493  0.95    LR   40
# 15 0.7973070 0.001686331 0.7940019 0.8006122  0.95    LR   50
# 16 0.6825276 0.001761671 0.6790748 0.6859804  0.95    RF    1
# 17 0.8061309 0.001624775 0.8029464 0.8093154  0.95    RF    2
# 18 0.8078385 0.001578268 0.8047451 0.8109318  0.95    RF    3
# 19 0.7969692 0.001646976 0.7937412 0.8001972  0.95    RF    4
# 20 0.8094017 0.001593652 0.8062782 0.8125252  0.95    RF    5
# 21 0.8234887 0.001544742 0.8204611 0.8265164  0.95    RF    6
# 22 0.8259862 0.001534107 0.8229794 0.8289930  0.95    RF    7
# 23 0.8322012 0.001511203 0.8292393 0.8351631  0.95    RF    8
# 24 0.8331447 0.001510321 0.8301845 0.8361049  0.95    RF    9
# 25 0.8335223 0.001508489 0.8305657 0.8364789  0.95    RF   10
# 26 0.8362177 0.001498771 0.8332801 0.8391552  0.95    RF   15
# 27 0.8442255 0.001450305 0.8413829 0.8470680  0.95    RF   20
# 28 0.8446873 0.001449458 0.8418465 0.8475282  0.95    RF   30
# 29 0.8447125 0.001448470 0.8418735 0.8475514  0.95    RF   40
# 30 0.8435809 0.001453397 0.8407323 0.8464295  0.95    RF   50

# nvar    intc intc_LCI intc_UCI slope slope_LCI slope_UCI
# <dbl>   <dbl>    <dbl>    <dbl> <dbl>     <dbl>     <dbl>
# 1     1 -0.0124   -0.158    0.131 1.00      0.723      1.29
# 2     2 -0.0125   -0.176    0.148 0.997     0.860      1.14
# 3     3 -0.0129   -0.177    0.148 0.996     0.859      1.14
# 4     4 -0.0133   -0.177    0.148 0.995     0.858      1.14
# 5     5 -0.0133   -0.178    0.148 0.993     0.856      1.13
# 6     6 -0.0116   -0.177    0.151 0.993     0.858      1.13
# 7     7 -0.0117   -0.177    0.151 0.990     0.855      1.13
# 8     8 -0.0134   -0.180    0.150 0.985     0.852      1.13
# 9     9 -0.0135   -0.180    0.150 0.985     0.851      1.12
# 10    10 -0.0133   -0.180    0.150 0.983     0.850      1.12
# 11    15 -0.0142   -0.181    0.150 0.968     0.836      1.11
# 12    20 -0.0144   -0.182    0.151 0.963     0.832      1.10
# 13    30 -0.0148   -0.184    0.151 0.936     0.808      1.07
# 14    40 -0.0153   -0.185    0.151 0.911     0.786      1.04
# 15    50 -0.0155   -0.186    0.152 0.896     0.771      1.03

#intercept <0, overestimation
#intercept >0, underestimation
#slopes <1, estimated risks too extreme in both directions
#slope >1 risk estimates to moderate
#slightly positive slope indicates slight underestimation

temp <- shigella0.5MUAC[["decilesCC"]][c("1","2","3","4","5","6","7","8","9","10")]
names(temp) <- c("1-var","2-var","3-var","4-var","5-var","6-var","7-var","8-var","9-var","10-var")  #renaming
#jpeg("",width=600,height=480,quality=400)
plot(x=seq(0,1.0,by=0.05),y=seq(0,1.0,by=0.05),type="l",
     xlab="Predicted Probability",ylab="Observed Proportion",
     main=expression(paste("Calibration Curve: Shigella in cases 0-59mo in GEMS")),
     cex.axis=1.5,cex.lab=1.5,cex.main=1.5)
points(temp$`2-var`$`mean(pred_glm)`,temp$`2-var`$`mean(true)`,col="black",pch=3,cex=2,lwd=2)
points(temp$`5-var`$`mean(pred_glm)`,temp$`5-var`$`mean(true)`,col="red",pch=1,cex=2,lwd=2)
points(temp$`10-var`$`mean(pred_glm)`,temp$`10-var`$`mean(true)`,col="blue",pch=2,cex=2,lwd=2)
legend("topleft",col=c("black","red","blue"),c("2-variable","5-variable","10-variable"),pch=c(3,1,2),cex=1.5)
dev.off()


AUC_df <- shigella0.5MUAC[["AUC_df"]]
#jpeg("",width=600,height=480,quality=100)
par(mar=c(5,5,4,2))
plot(AUC_df$nvar[1:length(shigella0.5MUAC[["nvars_opts"]])[1]],AUC_df$AUC[1:length(shigella0.5MUAC[["nvars_opts"]])[1]],
     xlab="number of variables",ylab="AUC",
     main="Shigella in cases 0-59mo in GEMS",
     ylim=c(0.5,1.0),
     pch=1,col="red",cex=2,lwd=2,cex.axis=1.5,cex.lab=1.5,cex.main=1.5)
points(AUC_df$nvar[1:length(shigella0.5MUAC[["nvars_opts"]])[1]],AUC_df$AUC[(length(shigella0.5MUAC[["nvars_opts"]])[1]+1):dim(AUC_df)[1]],
       pch=2,col="blue",cex=2,lwd=2)
legend("topleft",c("logistic reg","random forest"),col=c("red","blue"),pch=c(1,2),cex=1.5)
dev.off()



################### shigella muac ROC curve ####
#want a single ROC surve for CV in test datasets. use "result"
result<-shigella0.5MUAC[["result"]]

roc.data=result %>% split(.,list(.$nvar),drop=TRUE) %>% .$"2" %>% #swith this .$"x" for number of var; #subset to all the iter's (v-fold cross validations) for a given nvar (number of predictor variables in CPR)
  split(.,list(.$iter),drop=TRUE) #now have a list for each iter

#iter are the folds
# table($iter)
#>>> now  per iteration

predict_combo=list(NA)
true_combo=list(NA)

for (i in 1:shigella0.5MUAC[["iter"]]){ #this iter from the main RF/LR loop
  temp.list <- list(roc.data[[i]]$pred_glm)
  predict_combo <- c(predict_combo,temp.list)
  
  temp.list2 <- list(roc.data[[i]]$true)
  true_combo <- c(true_combo,temp.list2)
  
}
predict_combo=predict_combo[-1]
true_combo=true_combo[-1]
str(predict_combo)
str(true_combo)
combo <- list(predict_combo,true_combo)
names(combo) <- c("pred_glm","true")  #renaming
str(combo)

CV.roc.data <- cvAUC(predictions=combo$pred_glm, labels=combo$true) 

CV.roc.data.ShigellaMUAC.2 <- CV.roc.data
CV.roc.data.ShigellaMUAC.5 <- CV.roc.data
CV.roc.data.ShigellaMUAC.10 <- CV.roc.data
#to save
CV.roc.data.ShigellaMUAC <- list(CV.roc.data.ShigellaMUAC.2,CV.roc.data.ShigellaMUAC.5,CV.roc.data.ShigellaMUAC.10)
names(CV.roc.data.ShigellaMUAC) <- c("Shig.2","Shig.5","Shig.10")  #renaming
str(CV.roc.data.ShigellaMUAC)
#save(CV.roc.data.ShigellaMUAC, file = "")


#load(file = "ShigellaMUAC_CV.roc.data.Rdata")
#(there's a better way to do this as a list...)
CV.roc.data.Shig.2 <- CV.roc.data.ShigellaMUAC$Shig.2
CV.roc.data.Shig.5 <- CV.roc.data.ShigellaMUAC$Shig.5
CV.roc.data.Shig.10 <- CV.roc.data.ShigellaMUAC$Shig.10

#Plot CV AUC a single line of the averaged cross-validated ROC curve
#tiff("ShigellaMUAC_roc.tif",units="px",width=2000,height=2000,res=300)
plot(CV.roc.data.Shig.2$perf, avg="vertical", main="Cross-validated ROC Curves for Shigella",
     col="#1c61b6", lwd=2, lty=1) 
plot(CV.roc.data.Shig.5$perf, avg="vertical", col="#1c61b6", lwd=2, lty=2, add=TRUE) 
plot(CV.roc.data.Shig.10$perf, avg="vertical", col="#1c61b6", lwd=2, lty=3, add=TRUE) 
legend("bottomright", 
       legend = c("2-predictors","5-predictors","10-predictors"), 
       col = c("#1c61b6"),
       lty = c(1,2,3),
       lwd = 2)
segments(x0=0,y0=0.8,x1=0.42,y1=0.8,lty=2,col="gray")
segments(x0=0.42,y0=0.8,x1=0.42,y1=0.0,lty=2,col="gray")
dev.off()


################### shigella Se/Sp loop internal performance GEMS 0-59mo####
shigella<-SeSp.funct(data=cases,outcome="shigella_afe0.5",predictors=c("base_age","f4a_drh_blood"),new.data=cases_full,
                     POC.Se=0.9,POC.Sp=0.9,min=0,max=1.0,increment=0.05)
str(shigella)
summary(shigella[["glm"]])
shigella[["results"]]
shigella[["blood"]]

colors <- c("Se"="#7b3294", "Sp"="#c2a5cf", "FNR"="#008837", "FPR"="#a6dba0")
#jpeg("",width=600,height=480)
ggplot(shigella[["results"]], aes(x=prop.tested)) + 
  theme_bw() + ylim(0,1) +
  ggtitle("Accuracy of Clinical Prediction Rule (CPR)-Guided Testing for Shigella \n point-of-care Se=0.90; point-of-care Sp=0.90 \n Derived in GEMS 0-59mo, assessed in GEMS 0-59mo") +
  geom_smooth(aes(y=Se,color="Se",alpha = "CPR-derived"),method="loess",se=FALSE,size=1.2) + #se removes CI
  geom_smooth(aes(y=Sp,color="Sp",alpha = "BD-derived"),method="loess",se=FALSE,size=1.2) +
  geom_smooth(aes(y=FNR,color="FNR"),method="loess",se=FALSE,size=1.2) +
  geom_smooth(aes(y=FPR,color="FPR"),method="loess",se=FALSE,size=1.2) +
  geom_hline(yintercept=shigella[["blood"]][["Se.bld"]], linetype="dashed",color="#7b3294",size=1.1) +
  geom_hline(yintercept=shigella[["blood"]][["Sp.bld"]], linetype="dashed",color="#c2a5cf",size=1.1) +
  geom_hline(yintercept=shigella[["blood"]][["FNR.bld"]], linetype="dashed",color="#008837",size=1.1) +
  geom_hline(yintercept=shigella[["blood"]][["FPR.bld"]], linetype="dashed",color="#a6dba0",size=1.1) +
  labs(x="Proportion Tested",y="Test Accuracy")+
  scale_color_manual(name="Legend",
                     breaks=c("Se", "Sp", "FNR", "FPR"),
                     values=c("Se"="#7b3294", "Sp"="#c2a5cf", "FNR"="#008837", "FPR"="#a6dba0"),
                     guide = guide_legend(override.aes = list(size = 5))) +
  scale_alpha_manual(name = NULL,
                     values = c(1, 1),
                     breaks = c("CPR-derived", "BD-derived"),
                     guide = guide_legend(override.aes = list(linetype = c(1, 3),
                                                              shape = c(NA, NA),
                                                              color = "black")))
dev.off()


####################
#MAL-ED: dataset including etiology 
################### import data, define variables ####
load(file = "maled_tac_sadl.Rdata")
#looks like this is one observation per stool sample, including samples taken at diarrhea episode and monthly surveillance

temp1 <- maled_tac %>% filter(stooltype=="D1") %>% #now only stool samples from diarrhea episodes
  rename_with(~paste0(., "_Ct"), adenovirus_pan:trichuris) %>% #rename pathogen abbrev names to indicate are Ct values
  rename(astro = astrovirus_afe, #rename AFe's to abbreviations
         NoV = norovirus_afe,
         NoV.i = norovirus_gi_afe,
         NoV.ii = norovirus_gii_afe,
         rota = rotavirus_afe,
         sapo = sapovirus_afe,
         adeno = adenovirus_40_41_afe,
         aero = aeromonas_afe,
         campy = campylobacter_jejuni_coli_afe,
         crypto = cryptosporidium_afe,
         ehisto = e_histolytica_afe,
         cyclo = cyclospora_afe,
         hpylori = h_pylori_afe,
         isospora = isospora_afe,
         salm = salmonella_afe,
         shigella = shigella_eiec_afe,
         cholera = v_cholerae_afe,
         EAEC = EAEC_afe,
         ST_ETEC = ST_ETEC_afe,
         LT_ETEC = LT_ETEC_afe,
         TEPEC = tEPEC_afe,
         STEC = STEC_afe) %>%
  mutate(shigella_afe0.5 = ifelse(shigella>=0.5,1,0), #various AFe cutoffs
         shigella_afe0.3 = ifelse(shigella>=0.3,1,0),
         shigella_afe0.7 = ifelse(shigella>=0.7,1,0),
         cholera_afe0.5 = ifelse(cholera>=0.5,1,0),
         cholera_afe0.3 = ifelse(cholera>=0.3,1,0),
         cholera_afe0.7 = ifelse(cholera>=0.7,1,0),
         agegroup = ifelse(month_ss<12,1,
                           ifelse(month_ss>=12 & month_ss<24,2,
                                  ifelse(month_ss>=24 & month_ss<60,3,0)))) %>%
  rename(base_age=month_ss,
         f4a_drh_blood=maxb)

temp1$f4a_drh_blood <- as.factor(temp1$f4a_drh_blood)

################### into age groups ####
#want 1 observation per ID in each age group
data_MALED <- temp1 %>% group_by(pid) %>% sample_n(1) %>% ungroup 
data_age1_MALED <- temp1 %>% filter(agegroup==1) %>% 
  group_by(pid) %>% sample_n(1) %>% ungroup() 
data_age2_MALED <- temp1 %>% filter(agegroup == 2) %>% 
  group_by(pid) %>% sample_n(1) %>% ungroup()
data_age3_MALED <- temp1 %>% filter(agegroup ==3) %>% 
  group_by(pid) %>% sample_n(1) %>% ungroup()
data_age4_MALED <- temp1 %>% filter(agegroup==1 | agegroup==2) %>%
  group_by(pid) %>% sample_n(1) %>% ungroup()

length(unique(data_age1_MALED$pid)) 
length(unique(data_age2_MALED$pid)) 
length(unique(data_age3_MALED$pid)) 

################### descriptive for pub ####
table(data_MALED$site)
table(data_MALED$shigella_afe0.3,data_MALED$site)
round((table(data_MALED$shigella_afe0.3,data_MALED$site)[2,])/(table(data_MALED$site))*100,2)
table(data_MALED$shigella_afe0.5,data_MALED$site)
round((table(data_MALED$shigella_afe0.5,data_MALED$site)[2,])/(table(data_MALED$site))*100,2)
table(data_MALED$shigella_afe0.7,data_MALED$site)
round((table(data_MALED$shigella_afe0.7,data_MALED$site)[2,])/(table(data_MALED$site))*100,2)
table(data_MALED$cholera_afe0.3,data_MALED$site)
round((table(data_MALED$cholera_afe0.3,data_MALED$site)[2,])/(table(data_MALED$site))*100,2)
table(data_MALED$cholera_afe0.5,data_MALED$site)
round((table(data_MALED$cholera_afe0.5,data_MALED$site)[2,])/(table(data_MALED$site))*100,2)
table(data_MALED$cholera_afe0.7,data_MALED$site)
round((table(data_MALED$cholera_afe0.7,data_MALED$site)[2,])/(table(data_MALED$site))*100,2)

dim(data_MALED)
table(data_MALED$shigella_afe0.3)
round((table(data_MALED$shigella_afe0.3)[2])/(dim(data_MALED)[1])*100,2)
table(data_MALED$shigella_afe0.5)
round((table(data_MALED$shigella_afe0.5)[2])/(dim(data_MALED)[1])*100,2)
table(data_MALED$shigella_afe0.7)
round((table(data_MALED$shigella_afe0.7)[2])/(dim(data_MALED)[1])*100,2)

table(data_MALED$agegroup)
table(data_MALED$shigella_afe0.3,data_MALED$agegroup)
round((table(data_MALED$shigella_afe0.3,data_MALED$agegroup)[2,])/(table(data_MALED$agegroup))*100,2)
table(data_MALED$shigella_afe0.5,data_MALED$agegroup)
round((table(data_MALED$shigella_afe0.5,data_MALED$agegroup)[2,])/(table(data_MALED$agegroup))*100,2)
table(data_MALED$shigella_afe0.7,data_MALED$agegroup)
round((table(data_MALED$shigella_afe0.7,data_MALED$agegroup)[2,])/(table(data_MALED$agegroup))*100,2)
table(data_MALED$cholera_afe0.3,data_MALED$agegroup)
round((table(data_MALED$cholera_afe0.3,data_MALED$agegroup)[2,])/(table(data_MALED$agegroup))*100,2)
table(data_MALED$cholera_afe0.5,data_MALED$agegroup)
round((table(data_MALED$cholera_afe0.5,data_MALED$agegroup)[2,])/(table(data_MALED$agegroup))*100,2)
table(data_MALED$cholera_afe0.7,data_MALED$agegroup)
round((table(data_MALED$cholera_afe0.7,data_MALED$agegroup)[2,])/(table(data_MALED$agegroup))*100,2)

################### external validation 2var: shigella0.5 MUAC: GEMS 0-59mo ####
GEMS_glm_shigella <- glm(shigella_afe0.5~base_age+f4a_drh_blood,
                      data=cases,family="binomial",control=glm.control(maxit=50))
summary(GEMS_glm_shigella)
# Estimate Std. Error z value Pr(>|z|)    
# (Intercept)    -2.516451   0.076900  -32.72   <2e-16 ***
#   base_age        0.038278   0.002761   13.87   <2e-16 ***
#   f4a_drh_blood1  2.254699   0.076578   29.44   <2e-16 ***
round(exp(coef(GEMS_glm_shigella)),4)
# (Intercept)       base_age f4a_drh_blood1 
# 0.0807         1.0390         9.5324 
#save(GEMS_glm_shigella, file = "")
round(exp(confint(GEMS_glm_shigella)),4)
# 2.5 %  97.5 %
#   (Intercept)    0.0693  0.0937
# base_age       1.0334  1.0447
# f4a_drh_blood1 8.2101 11.0850

GEMS_2var<-cases %>% select(shigella_afe0.5,base_age,f4a_drh_blood)
GEMS_2var$pred_glm <- as.numeric(predict(GEMS_glm_shigella,newdata=GEMS_2var,type="response"))
GEMS_AUC <- AUC(predictions=GEMS_2var$pred_glm,labels=GEMS_2var$shigella_afe0.5)
GEMS_AUC
# [1] 0.8031113
#using pROC package since need CI
GEMS_AUC <- roc(response=GEMS_2var$shigella_afe0.5,predictor=GEMS_2var$pred_glm)
paste(round(GEMS_AUC$auc,2)," (",
      round(ci.auc(GEMS_AUC)[1],2),", ",
      round(ci.auc(GEMS_AUC)[3],2),")",sep="")
# "0.8 (0.79, 0.82)"


GEMS_decilesCC <- GEMS_2var %>% mutate(decile_glm=ntile(pred_glm,10)) %>%
  group_by(decile_glm) %>% summarize(mean(shigella_afe0.5),mean(pred_glm))
#save(GEMS_decilesCC, file = "")
#load(file = "ShigellaMUAC_GEMS059_decilesCC.Rdata")


data_MALED$GEMS_pred_glm <- as.numeric(predict(GEMS_glm_shigella,newdata=data_MALED,type="response"))
MALED_AUC <- AUC(predictions=data_MALED$GEMS_pred_glm,labels=data_MALED$shigella_afe0.5)
MALED_AUC
# [1] 0.7817442
MALED_AUC <- roc(response=data_MALED$shigella_afe0.5,predictor=data_MALED$GEMS_pred_glm)
paste(round(MALED_AUC$auc,2)," (",
      round(ci.auc(MALED_AUC)[1],2),", ",
      round(ci.auc(MALED_AUC)[3],2),")",sep="")
# "0.78 (0.75, 0.81)"

MALED_decilesCC <- data_MALED %>% mutate(decile_glm=ntile(GEMS_pred_glm,10)) %>%
  group_by(decile_glm) %>% summarize(mean(shigella_afe0.5),mean(GEMS_pred_glm))
#save(MALED_decilesCC, file = "")
#load(file = "ShigellaMUAC_MALED_decilesCC_059.Rdata")


#jpeg("",width=600,height=480,quality=100)
#tiff("ShigellaMUAC_CC_ExternalVal_GEMS059.tif",units="px",width=2400,height=1920,res=300)
plot(x=seq(0,1,by=0.1),y=seq(0,1,by=0.1),type="l",
     xlab="Predicted Probability",ylab="Observed Proportion",
     main=expression(paste("Calibration Curve: Shigella AFe0.5")),
     xlim=c(0,1.0),ylim=c(0,1.0))
title("GEMS-derived CPR",line=0.7,font.main=1)
points(GEMS_decilesCC$`mean(pred_glm)`,GEMS_decilesCC$`mean(shigella_afe0.5)`,col="red",pch=1)
points(MALED_decilesCC$`mean(GEMS_pred_glm)`,MALED_decilesCC$`mean(shigella_afe0.5)`,col="blue",pch=2)
legend("topleft",col=c("red","blue"),c("GEMS data (0-59mo), AUC=0.80 (0.79, 0.82)","MALED data (0-23mo), AUC=0.78 (0.75, 0.81)"),pch=c(1,2))
dev.off()

#calib
intercept <- glm(shigella_afe0.5~1,offset=log(GEMS_pred_glm/(1-GEMS_pred_glm)),family="binomial",data=data_MALED)
summary(intercept)
confint(intercept)
slope <- glm(shigella_afe0.5~log(GEMS_pred_glm/(1-GEMS_pred_glm)),family="binomial",data=data_MALED)
summary(slope)
confint(slope)

# Estimate Std. Error z value Pr(>|z|)
# (Intercept) -0.02185    0.08667  -0.252    0.801
# 2.5 %     97.5 % 
#   -0.1950427  0.1448748 

# Estimate Std. Error z value Pr(>|z|)    
# (Intercept)                             0.02429    0.22532   0.108    0.914    
# log(GEMS_pred_glm/(1 - GEMS_pred_glm))  1.02597    0.11703   8.766   <2e-16 ***
#   2.5 %    97.5 %
#   (Intercept)                            -0.4149618 0.4721205
# log(GEMS_pred_glm/(1 - GEMS_pred_glm))  0.7985008 1.2588695



################### Se/Sp shigella0.5 MUAC look external performance: derive:GEMS 0-59mo assess:MALED ####
shigella<-SeSp.funct(data=cases,outcome="shigella_afe0.5",predictors=c("base_age","f4a_drh_blood"),new.data=data_MALED,
                     POC.Se=0.9,POC.Sp=0.9,min=0,max=1.0,increment=0.01)
str(shigella)
summary(shigella[["glm"]])
shigella[["results"]]
results<-shigella[["results"]]
results_dysentery<-shigella[["blood"]]
#save(shigella, file = "Shigella_GEMS059_SeSpFPRFNR_POC0.90.0_MALED.Rdata")
#load(file = "Shigella_GEMS059_SeSpFPRFNR_POC0.90.0_MALED.Rdata")


#how few need testing before start doing better than dysentery-only regimen
shigella[["blood"]]
# TN.bld  TP.bld  FN.bld  FP.bld  Se.bld  Sp.bld PPV.bld NPV.bld FPR.bld FNR.bld 
# 1067.00   26.00  120.00   46.00    0.18    0.96    0.36    0.90    0.04    0.82 
# TN.bld  TP.bld  FN.bld  FP.bld  Se.bld  Sp.bld PPV.bld NPV.bld FPR.bld FNR.bld 
# 1066.00   27.00  118.00   48.00    0.19    0.96    0.36    0.90    0.04    0.81 
temp <- results[which(results$Se>=0.18),]


# Estimate Std. Error z value Pr(>|z|)    
# (Intercept)    -2.516451   0.076900  -32.72   <2e-16 ***
#   base_age        0.038278   0.002761   13.87   <2e-16 ***
#   f4a_drh_blood1  2.254699   0.076578   29.44   <2e-16 ***

colors <- c("Se"="#7b3294", "Sp"="#c2a5cf", "FNR"="#008837", "FPR"="#a6dba0")
#tiff("Shigella_GEMS059_SeSpFPRFNR_POC0.90.9_MALED.tif",units="px",width=2400,height=1920,res=300)
ggplot(shigella[["results"]], aes(x=prop.tested)) + 
  theme_bw() + ylim(0,1) +
  ggtitle("Accuracy of Clinical Prediction Rule (CPR)-Guided Testing for Shigella \n point-of-care Se=0.90; point-of-care Sp=0.90 \n Derived in GEMS 0-59mo, assessed in MALED 0-23mo") +
  geom_smooth(aes(y=Se,color="Se",alpha = "CPR-derived"),method="loess",se=FALSE,size=1.2) + 
  geom_smooth(aes(y=Sp,color="Sp",alpha = "BD-derived"),method="loess",se=FALSE,size=1.2) +
  geom_smooth(aes(y=FNR,color="FNR"),method="loess",se=FALSE,size=1.2) +
  geom_smooth(aes(y=FPR,color="FPR"),method="loess",se=FALSE,size=1.2) +
  geom_hline(yintercept=shigella[["blood"]][["Se.bld"]], linetype="dashed",color="#7b3294",size=1.1) +
  geom_hline(yintercept=shigella[["blood"]][["Sp.bld"]], linetype="dashed",color="#c2a5cf",size=1.1) +
  geom_hline(yintercept=shigella[["blood"]][["FNR.bld"]], linetype="dashed",color="#008837",size=1.1) +
  geom_hline(yintercept=shigella[["blood"]][["FPR.bld"]], linetype="dashed",color="#a6dba0",size=1.1) +
  labs(x="Proportion Receiving POC Diagnostic Test",y="Diagnostic Regimen Accuracy")+
  scale_color_manual(name=NULL,
                     breaks=c("Se", "Sp", "FNR", "FPR"),
                     values=c("Se"="#7b3294", "Sp"="#c2a5cf", "FNR"="#008837", "FPR"="#a6dba0"),
                     guide = guide_legend(override.aes = list(size = 5))) +
  scale_alpha_manual(name = NULL,
                     labels = c("CPR regimen", "dysentery regimen"),
                     values = c(1, 1),
                     breaks = c("CPR-derived", "BD-derived"),
                     guide = guide_legend(override.aes = list(linetype = c(1, 3),
                                                              shape = c(NA, NA),
                                                              color = "black")))
dev.off()

glm<-glm(shigella_afe0.5~base_age+f4a_drh_blood, data=cases, family="binomial")
summary(glm)
# Estimate Std. Error z value Pr(>|z|)    
# (Intercept)    -2.516451   0.076900  -32.72   <2e-16 ***
#   base_age        0.038278   0.002761   13.87   <2e-16 ***
#   f4a_drh_blood1  2.254699   0.076578   29.44   <2e-16 ***

#example of each possible kid
representative=data.frame(base_age=c(seq(from=0,to=24,1),seq(from=0,to=24,1)),
                          f4a_drh_blood=as.factor(c(rep(0,25),rep(1,25))))
representative
str(representative)
representative$predicted <- predict(glm, newdata=representative, type="response")
representative$test <- ifelse(representative$predicted<0.13 | representative$predicted>0.55, 0,
                              1)
representative
representative.ordered<-representative[order(representative$base_age,representative$f4a_drh_blood),]
write.csv(representative.ordered, file = "representative.ordered.csv")




################### shigella overlapping histograms of predicted probabilities ####
val_data<-shigella[["val_data"]]
p1 <- hist(val_data[which(val_data$shigella_afe0.5==0),]$predicted, freq=F)
p2 <- hist(val_data[which(val_data$shigella_afe0.5==1),]$predicted, freq=F)
#by actual bloody diarrhea status
p3 <- hist(val_data[which(val_data$f4a_drh_blood==0),]$predicted, freq=F, density=15, angle=45, col="darkgray")
p4 <- hist(val_data[which(val_data$f4a_drh_blood==1),]$predicted, freq=F, density=15, angle=135, add=T)

#jpeg("",width=600,height=480)
plot( p1, col=rgb(1,0,0,1/4), freq=F, xlim=c(0,1.0),ylim=c(0,12),xlab="predicted probability of infection",main="Predicted Probability of \nShigella Infection, by True Status")
plot( p2, col=rgb(0,0,1,1/4), freq=F, xlim=c(0,1.0),ylim=c(0,12), add=T)
legend('topright',
       legend=c('Truly NOT Shigella','Truly Shigella'),
       fill = c(rgb(1,0,0,1/4),rgb(0,0,1,1/4)), 
       bty = 'n')
dev.off()


################### external validation 2var: shigella0.5 MUAC: GEMS 0-23mo ####
GEMS_glm_shigella <- glm(shigella_afe0.5~base_age+f4a_drh_blood,
                         data=cases_age4,family="binomial",control=glm.control(maxit=50))
summary(GEMS_glm_shigella)
# Estimate Std. Error z value Pr(>|z|)    
# (Intercept)    -4.120700   0.149313  -27.60   <2e-16 ***
#   base_age        0.165558   0.009065   18.26   <2e-16 ***
#   f4a_drh_blood1  2.198486   0.102454   21.46   <2e-16 ***
round(coef(GEMS_glm_shigella),4)
#save(GEMS_glm_shigella, file = "")
round(exp(coef(glm(shigella_afe0.5~base_age+f4a_drh_blood,
                   data=cases_age4,family="binomial",control=glm.control(maxit=50)))),4)
# (Intercept)       base_age f4a_drh_blood1 
# 0.0162         1.1801         9.0114

GEMS_2var<-cases %>% select(shigella_afe0.5,base_age,f4a_drh_blood)
GEMS_2var$pred_glm <- as.numeric(predict(GEMS_glm_shigella,newdata=GEMS_2var,type="response"))
GEMS_AUC <- AUC(predictions=GEMS_2var$pred_glm,labels=GEMS_2var$shigella_afe0.5)
GEMS_AUC
# [1] 0.7625685
#using pROC package since need CI
GEMS_AUC <- roc(response=GEMS_2var$shigella_afe0.5,predictor=GEMS_2var$pred_glm)
paste(round(GEMS_AUC$auc,2)," (",
      round(ci.auc(GEMS_AUC)[1],2),", ",
      round(ci.auc(GEMS_AUC)[3],2),")",sep="")
# "0.76 (0.75, 0.78)"


GEMS_decilesCC <- GEMS_2var %>% mutate(decile_glm=ntile(pred_glm,10)) %>%
  group_by(decile_glm) %>% summarize(mean(shigella_afe0.5),mean(pred_glm))
#save(GEMS_decilesCC, file = "")

data_MALED$GEMS_pred_glm <- as.numeric(predict(GEMS_glm_shigella,newdata=data_MALED,type="response"))
MALED_AUC <- AUC(predictions=data_MALED$GEMS_pred_glm,labels=data_MALED$shigella_afe0.5)
MALED_AUC
# [1] 0.7990255
MALED_AUC <- roc(response=data_MALED$shigella_afe0.5,predictor=data_MALED$GEMS_pred_glm)
paste(round(MALED_AUC$auc,2)," (",
      round(ci.auc(MALED_AUC)[1],2),", ",
      round(ci.auc(MALED_AUC)[3],2),")",sep="")
# "0.77 (0.74, 0.81)"
# "0.8 (0.77, 0.83)"

MALED_decilesCC <- data_MALED %>% mutate(decile_glm=ntile(GEMS_pred_glm,10)) %>%
  group_by(decile_glm) %>% summarize(mean(shigella_afe0.5),mean(GEMS_pred_glm))
#save(MALED_decilesCC, file = "")

#jpeg("",width=600,height=480,quality=100)
plot(x=seq(0,1,by=0.1),y=seq(0,1,by=0.1),type="l",
     xlab="Predicted Probability",ylab="Observed Proportion",
     main=expression(paste("Calibration Curve: Shigella AFe0.5")),
     xlim=c(0,1.0),ylim=c(0,1.0))
title("GEMS-derived CPR",line=0.7,font.main=1)
points(GEMS_decilesCC$`mean(pred_glm)`,GEMS_decilesCC$`mean(shigella_afe0.5)`,col="red",pch=1)
points(MALED_decilesCC$`mean(GEMS_pred_glm)`,MALED_decilesCC$`mean(shigella_afe0.5)`,col="blue",pch=2)
legend("topleft",col=c("red","blue"),c("GEMS data (0-23mo), AUC=0.76 (0.75, 0.78)","MALED data (0-23mo), AUC=0.77 (0.74, 0.81)"),pch=c(1,2))
dev.off()

#calib
intercept <- glm(shigella_afe0.5~1,offset=log(GEMS_pred_glm/(1-GEMS_pred_glm)),family="binomial",data=kilifi)
summary(intercept)
confint(intercept)
slope <- glm(shigella_afe0.5~log(GEMS_pred_glm/(1-GEMS_pred_glm)),family="binomial",data=kilifi)
summary(slope)
confint(slope)

# Estimate Std. Error z value Pr(>|z|)  
# (Intercept) -0.18071    0.09117  -1.982   0.0475 *
#   2.5 %      97.5 % 
#   -0.36313587 -0.00550114 
# Estimate Std. Error z value Pr(>|z|)    
# (Intercept)                             -0.3017     0.2044  -1.476     0.14    
# log(GEMS_pred_glm/(1 - GEMS_pred_glm))   0.9289     0.1074   8.647   <2e-16 ***
#   2.5 %     97.5 %
#   (Intercept)                            -0.7046998 0.09958622
# log(GEMS_pred_glm/(1 - GEMS_pred_glm))  0.7182880 1.14057865




################### Se/Sp shigella0.5 MUAC look external performance: derive:GEMS 0-23mo assess:MALED ####
shigella<-SeSp.funct(data=cases_age4,outcome="shigella_afe0.5",predictors=c("base_age","f4a_drh_blood"),new.data=data_MALED,
                     POC.Se=0.9,POC.Sp=0.9,min=0,max=1.0,increment=0.05)
str(shigella)
summary(shigella[["glm"]])
shigella[["results"]]
results<-shigella[["results"]]
results_dysentery<-shigella[["blood"]]
shigella[["blood"]]
# Estimate Std. Error z value Pr(>|z|)    
# (Intercept)    -4.120700   0.149313  -27.60   <2e-16 ***
#   base_age        0.165558   0.009065   18.26   <2e-16 ***
#   f4a_drh_blood1  2.198486   0.102454   21.46   <2e-16 ***
  
colors <- c("Se"="#7b3294", "Sp"="#c2a5cf", "FNR"="#008837", "FPR"="#a6dba0")
#jpeg("",width=600,height=480)
ggplot(shigella[["results"]], aes(x=prop.tested)) + 
  theme_bw() + ylim(0,1) +
  ggtitle("Accuracy of Clinical Prediction Rule (CPR)-Guided Testing for Shigella \n point-of-care Se=0.90; point-of-care Sp=0.90 \n Derived in GEMS 0-23mo, assessed in MALED 0-23mo") +
  geom_smooth(aes(y=Se,color="Se",alpha = "CPR-derived"),method="loess",se=FALSE,size=1.2) + #se removes CI
  geom_smooth(aes(y=Sp,color="Sp",alpha = "BD-derived"),method="loess",se=FALSE,size=1.2) +
  geom_smooth(aes(y=FNR,color="FNR"),method="loess",se=FALSE,size=1.2) +
  geom_smooth(aes(y=FPR,color="FPR"),method="loess",se=FALSE,size=1.2) +
  geom_hline(yintercept=shigella[["blood"]][["Se.bld"]], linetype="dashed",color="#7b3294",size=1.1) +
  geom_hline(yintercept=shigella[["blood"]][["Sp.bld"]], linetype="dashed",color="#c2a5cf",size=1.1) +
  geom_hline(yintercept=shigella[["blood"]][["FNR.bld"]], linetype="dashed",color="#008837",size=1.1) +
  geom_hline(yintercept=shigella[["blood"]][["FPR.bld"]], linetype="dashed",color="#a6dba0",size=1.1) +
  labs(x="Proportion Tested",y="Test Accuracy")+
  scale_color_manual(name="Legend",
                     breaks=c("Se", "Sp", "FNR", "FPR"),
                     values=c("Se"="#7b3294", "Sp"="#c2a5cf", "FNR"="#008837", "FPR"="#a6dba0"),
                     guide = guide_legend(override.aes = list(size = 5))) +
  scale_alpha_manual(name = NULL,
                     values = c(1, 1),
                     breaks = c("CPR-derived", "BD-derived"),
                     guide = guide_legend(override.aes = list(linetype = c(1, 3),
                                                              shape = c(NA, NA),
                                                              color = "black")))
dev.off()

####################
#additional GEMS sensitivity analyses 
####################
#shigella AFe>=0.5 only HAZ
################### var screening, AUC ####
shigella0.5HAZ <- CPR.funct(data=cases,outcome="shigella_afe0.5",iter=100,nvars_opts=c(2,5,10))
shigella0.5HAZ[["df_imps"]]
shigella0.5HAZ[["AUC_df"]]
shigella0.5HAZ[["calib"]]

# names    var_red
# 1               base_age 89.9527846
# 2          f4a_drh_blood 78.6048242
# 3                f4b_haz 34.8842720
# 4               f4b_resp 31.5690078
# 5               f4b_temp 29.5257957
# 6               f4b_eyes 27.5583261
# 7          f4a_ppl_house 23.6390769
# 8                   site 21.4553894
# 9           f4a_drh_days 15.9615761
# 10         f4a_slp_rooms 14.4498337
# 11        f4a_drh_strain 13.6439955
# 12         f4a_share_fac 13.4069449
# 13      f4a_yng_children 12.4553014
# 14         f4a_prim_schl 11.6621515
# 15         f4a_drh_vomit 11.5400569
# 16         f4a_breastfed 11.2696776
# 17        f4a_fuel_grass 10.7955309
# 18             f4b_mouth 10.0273774
# 19        f4a_max_stools  9.1047632
# 20        f4a_offr_drink  8.8086259
# 21         f4b_recommend  8.7856863
# 22         f4a_fac_waste  8.0461620
# 23         f4a_fuel_crop  7.9036261
# 24       f4a_store_water  7.6651874
# 25     f4a_drh_bellypain  7.5349246
# 26       f4a_cur_thirsty  7.5305264
# 27         f4a_drh_cough  7.1212277
# 28       f4a_water_avail  7.0486891
# 29        f4a_disp_feces  7.0462691
# 30        f4a_drh_thirst  6.7006662
# 31          f4a_dad_live  6.6839561
# 32        f4a_trt_method  6.6198215
# 33          f4a_ms_water  6.4935744
# 34            f4b_mental  5.5979812
# 35      f4a_cur_drymouth  5.4563062
# 36      f4a_hometrt_none  5.3069931
# 37          f4a_wash_use  5.2039806
# 38             f3_gender  5.1472849
# 39       f4a_hometrt_ors  5.0945606
# 40        f4a_wash_child  5.0142888
# 41  f4a_drh_lethrgy_miss  4.8796307
# 42         f4a_wash_cook  4.8485712
# 43        f4a_wash_nurse  4.8297888
# 44        f4a_house_tele  4.6753737
# 45              f4b_skin  4.6384028
# 46        f4a_house_bike  4.6382894
# 47      f4a_seek_outside  4.6057022
# 48          f4a_ani_fowl  4.6007058
# 49          f4a_wash_def  4.5402556
# 50      f4a_drh_restless  4.4637691
# 51         f4a_fuel_dung  4.4246426
# 52             f4a_floor  4.4092494
# 53      f4a_relationship  4.3952091
# 54     f4a_fuel_charcoal  4.3904581
# 55           f4a_ani_cat  4.3837830
# 56       f4a_house_phone  4.3610305
# 57          f4a_wash_eat  4.3139771
# 58    f4a_water_deepwell  4.2787287
# 59      f4a_cur_restless  4.2384285
# 60       f4a_house_radio  4.2297203
# 61             f4b_admit  4.2101804
# 62         f4a_fuel_wood  4.1634855
# 63          f4a_cur_skin  4.1513633
# 64           f3_drh_hosp  4.0801897
# 65           f4a_ani_dog  3.9473678
# 66       f4a_ani_rodents  3.9472327
# 67        f4a_house_elec  3.9438399
# 68        f4a_seek_pharm  3.9409486
# 69      f4a_water_pubtap  3.9303695
# 70       f4a_house_scoot  3.8620238
# 71      f4a_house_fridge  3.8285426
# 72      f4a_house_agland  3.7883378
# 73          f4a_ani_goat  3.6423968
# 74         f4a_trt_water  3.6053924
# 75     f4a_hometrt_othr1  3.4611225
# 76           f4a_ani_cow  3.4489496
# 77             f3_drh_iv  3.4309675
# 78     f4a_hometrt_maize  3.4106875
# 79        f4b_under_nutr  3.4074608
# 80         f3_drh_turgor  3.3179158
# 81    f4a_cur_fastbreath  3.2711334
# 82     f4a_drh_lessdrink  3.2114751
# 83        f4a_hometrt_ab  3.0127638
# 84       f4a_fuel_natgas  2.8788027
# 85   f4a_water_shallwell  2.8719863
# 86            f4a_ani_no  2.8457302
# 87         f4a_ani_sheep  2.5964012
# 88         f4a_wash_othr  2.5433945
# 89       f4a_wash_animal  2.5225045
# 90        f4a_water_yard  2.4469389
# 91      f4a_hometrt_herb  2.3811821
# 92         f4a_house_car  2.3630118
# 93          f4a_drh_conv  2.3603243
# 94         f4a_ani_other  2.2235486
# 95        f4a_house_cart  2.2156302
# 96       f4a_water_house  2.1563281
# 97      f4a_water_bought  2.1197096
# 98       f4a_seek_healer  1.9553701
# 99         f4a_fuel_kero  1.7754216
# 100    f4a_water_pubwell  1.7673294
# 101     f4a_hometrt_zinc  1.7017337
# 102     f4a_fuel_propane  1.5578058
# 103       f4a_water_bore  1.5248904
# 104       f4a_house_none  1.5165556
# 105     f4a_seek_privdoc  1.5110347
# 106       f4a_seek_remdy  1.4796681
# 107         f4b_abn_hair  1.4530767
# 108       f4a_water_well  1.4040951
# 109         f4a_seek_doc  1.3583551
# 110        f4a_fuel_coal  1.2827461
# 111       f4a_water_rain  1.2358480
# 112       f4a_fuel_other  1.1833578
# 113        f4a_drh_consc  1.1551589
# 114    f4a_hometrt_othr2  1.1534799
# 115        f4a_fuel_elec  1.1355864
# 116     f4a_drh_prolapse  1.1278585
# 117     f4a_hometrt_milk  1.1228978
# 118      f4b_chest_indrw  1.1171168
# 119       f4b_skin_flaky  1.1133770
# 120    f4a_water_covwell  1.0203741
# 121   f4a_water_covpwell  0.8496377
# 122      f4a_fuel_biogas  0.7949238
# 123  f4a_hometrt_othrliq  0.7578056
# 124      f4a_water_river  0.7408249
# 125       f4a_water_othr  0.7376825
# 126       f4a_house_boat  0.6878855
# 127       f4a_seek_other  0.6863021
# 128       f4a_water_pond  0.6061334
# 129          f4b_bipedal  0.5545301
# 130           f4b_rectal  0.4993921
# 131   f4a_water_unspring  0.3528463
# 132      f4a_seek_friend  0.2876824
# 133  f4a_water_prospring  0.1622578

# AUC          SE     lower     upper level Model nvar
# 1 0.8026314 0.001616921 0.7994623 0.8058005  0.95    LR    2
# 2 0.7923219 0.001679061 0.7890310 0.7956128  0.95    LR    5
# 3 0.7970040 0.001679637 0.7937119 0.8002960  0.95    LR   10
# 4 0.8063299 0.001623963 0.8031470 0.8095128  0.95    RF    2
# 5 0.8146255 0.001576203 0.8115362 0.8177148  0.95    RF    5
# 6 0.8337590 0.001494597 0.8308296 0.8366884  0.95    RF   10

# nvar     intc intc_LCI intc_UCI slope slope_LCI slope_UCI
# <dbl>    <dbl>    <dbl>    <dbl> <dbl>     <dbl>     <dbl>
#   1     2 -0.00364   -0.167    0.157 1.01      0.867      1.15
# 2     5 -0.00244   -0.167    0.159 0.996     0.859      1.14
# 3    10 -0.00336   -0.170    0.160 0.980     0.847      1.12

####################
#shigella AFe>=0.5 HAZ+MUAC
################### var screening, AUC ####
names<-append(x=names, values=c("f4b_muac"))

shigella0.5HAZMUAC <- CPR.funct(data=cases,outcome="shigella_afe0.5",iter=100,nvars_opts=c(2,5,10))
shigella0.5HAZMUAC[["df_imps"]]
shigella0.5HAZMUAC[["AUC_df"]]
shigella0.5HAZMUAC[["calib"]]

# names    var_red
# 1               base_age 86.2531197
# 2          f4a_drh_blood 80.0020013
# 3                f4b_haz 33.3409382
# 4               f4b_muac 31.6025091
# 5               f4b_resp 30.6090870
# 6               f4b_temp 28.5759929
# 7               f4b_eyes 25.6730699
# 8          f4a_ppl_house 21.9301753
# 9                   site 21.3055396
# 10          f4a_drh_days 15.5164642
# 11        f4a_drh_strain 14.3207340
# 12         f4a_slp_rooms 13.7252622
# 13         f4a_share_fac 12.8765973
# 14      f4a_yng_children 11.9472586
# 15         f4a_drh_vomit 11.1270171
# 16         f4a_prim_schl 10.9632425
# 17         f4a_breastfed 10.8371064
# 18        f4a_fuel_grass 10.5865508
# 19             f4b_mouth  9.9430325
# 20        f4a_max_stools  9.1513127
# 21         f4b_recommend  8.7750973
# 22        f4a_offr_drink  8.5033990
# 23       f4a_store_water  8.3505170
# 24         f4a_fuel_crop  7.8704382
# 25         f4a_fac_waste  7.4922003
# 26     f4a_drh_bellypain  7.4289011
# 27       f4a_water_avail  6.8137451
# 28         f4a_drh_cough  6.7830777
# 29        f4a_disp_feces  6.7415287
# 30       f4a_cur_thirsty  6.7026684
# 31        f4a_trt_method  6.5044626
# 32        f4a_drh_thirst  6.3658490
# 33          f4a_ms_water  6.3177998
# 34          f4a_dad_live  6.2888291
# 35            f4b_mental  5.3056588
# 36      f4a_cur_drymouth  5.1921130
# 37         f4a_fuel_dung  4.9852833
# 38      f4a_hometrt_none  4.9492145
# 39       f4a_hometrt_ors  4.8776267
# 40             f3_gender  4.8365022
# 41          f4a_wash_use  4.7811823
# 42         f4a_wash_cook  4.6969843
# 43        f4a_wash_child  4.6767505
# 44        f4a_wash_nurse  4.6463865
# 45              f4b_skin  4.4964213
# 46      f4a_seek_outside  4.4804456
# 47  f4a_drh_lethrgy_miss  4.4204358
# 48          f4a_wash_def  4.4071609
# 49        f4a_house_tele  4.4068749
# 50        f4a_house_bike  4.3762851
# 51          f4a_ani_fowl  4.2224055
# 52         f4a_fuel_wood  4.1721743
# 53       f4a_house_radio  4.1626565
# 54          f4a_wash_eat  4.1493203
# 55           f4a_ani_cat  4.1341867
# 56     f4a_fuel_charcoal  4.1237955
# 57    f4a_water_deepwell  4.1172660
# 58       f4a_house_phone  4.1098052
# 59      f4a_drh_restless  4.1097047
# 60      f4a_cur_restless  4.0673006
# 61          f4a_cur_skin  4.0266994
# 62             f4b_admit  4.0147363
# 63           f3_drh_hosp  4.0140016
# 64      f4a_relationship  4.0004044
# 65             f4a_floor  3.9715646
# 66           f4a_ani_dog  3.8576045
# 67        f4a_seek_pharm  3.8191624
# 68        f4a_house_elec  3.7640291
# 69       f4a_ani_rodents  3.7449318
# 70      f4a_water_pubtap  3.7372806
# 71       f4a_house_scoot  3.6176749
# 72      f4a_house_agland  3.5121521
# 73      f4a_house_fridge  3.4798417
# 74          f4a_ani_goat  3.4496366
# 75         f4a_trt_water  3.4424920
# 76           f4a_ani_cow  3.3952526
# 77     f4a_hometrt_othr1  3.3286709
# 78             f3_drh_iv  3.2956143
# 79     f4a_drh_lessdrink  3.2929433
# 80         f3_drh_turgor  3.2751254
# 81     f4a_hometrt_maize  3.2577734
# 82    f4a_cur_fastbreath  3.2483491
# 83        f4b_under_nutr  2.9998599
# 84        f4a_hometrt_ab  2.9826831
# 85       f4a_fuel_natgas  2.7380709
# 86            f4a_ani_no  2.7133353
# 87   f4a_water_shallwell  2.6755615
# 88         f4a_wash_othr  2.6317929
# 89         f4a_ani_sheep  2.5755524
# 90      f4a_hometrt_herb  2.3880747
# 91       f4a_wash_animal  2.3836709
# 92        f4a_water_yard  2.3408081
# 93          f4a_drh_conv  2.3271923
# 94         f4a_house_car  2.3138960
# 95         f4a_ani_other  2.2052495
# 96      f4a_water_bought  2.1583895
# 97        f4a_house_cart  2.1199812
# 98       f4a_water_house  2.1060394
# 99       f4a_seek_healer  1.8719775
# 100        f4a_fuel_kero  1.7875650
# 101    f4a_water_pubwell  1.6984183
# 102     f4a_seek_privdoc  1.5805483
# 103     f4a_fuel_propane  1.5173932
# 104     f4a_hometrt_zinc  1.5159222
# 105       f4a_house_none  1.4567644
# 106       f4a_water_bore  1.4481770
# 107       f4a_seek_remdy  1.4322805
# 108       f4a_water_well  1.3468202
# 109         f4b_abn_hair  1.3034272
# 110         f4a_seek_doc  1.2935629
# 111     f4a_hometrt_milk  1.2254630
# 112       f4a_fuel_other  1.2161605
# 113       f4a_water_rain  1.1880162
# 114        f4a_drh_consc  1.1497728
# 115        f4a_fuel_coal  1.1130611
# 116        f4a_fuel_elec  1.0927107
# 117    f4a_hometrt_othr2  1.0804274
# 118     f4a_drh_prolapse  1.0406255
# 119      f4b_chest_indrw  1.0197249
# 120       f4b_skin_flaky  0.9879977
# 121    f4a_water_covwell  0.9755378
# 122   f4a_water_covpwell  0.8187217
# 123       f4a_seek_other  0.7466230
# 124      f4a_water_river  0.7393747
# 125  f4a_hometrt_othrliq  0.7310826
# 126      f4a_fuel_biogas  0.7297437
# 127       f4a_water_othr  0.6785253
# 128       f4a_house_boat  0.5844313
# 129       f4a_water_pond  0.5600813
# 130           f4b_rectal  0.4422468
# 131          f4b_bipedal  0.4369654
# 132      f4a_seek_friend  0.3005243
# 133   f4a_water_unspring  0.2924222
# 134  f4a_water_prospring  0.1335331

# AUC          SE     lower     upper level Model nvar
# 1 0.8022894 0.001610594 0.7991327 0.8054461  0.95    LR    2
# 2 0.7988069 0.001606188 0.7956589 0.8019550  0.95    LR    5
# 3 0.7988490 0.001658525 0.7955983 0.8020996  0.95    LR   10
# 4 0.8059127 0.001617789 0.8027419 0.8090835  0.95    RF    2
# 5 0.8080782 0.001584222 0.8049732 0.8111833  0.95    RF    5
# 6 0.8340135 0.001489852 0.8310934 0.8369336  0.95    RF   10

# nvar   intc intc_LCI intc_UCI slope slope_LCI slope_UCI
# <dbl>  <dbl>    <dbl>    <dbl> <dbl>     <dbl>     <dbl>
#   1     2 0.0211   -0.142    0.181 1.00      0.865      1.15
# 2     5 0.0218   -0.142    0.182 1.00      0.863      1.14
# 3    10 0.0216   -0.144    0.184 0.986     0.853      1.13

####################
#shigella AFe>=0.5 0-11mo
################### var screening, AUC ####
names<-append(x=names, values=c("f4b_muac"))
names <- names[!names %in% c("f4b_haz")]

shigella0.5MUAC.011 <- CPR.funct(data=cases_age1,outcome="shigella_afe0.5",iter=100,nvars_opts=c(2,5,10))
shigella0.5MUAC.011[["df_imps"]]
shigella0.5MUAC.011[["AUC_df"]]
shigella0.5MUAC.011[["calib"]]

# names     var_red
# 1               base_age 6.875053963
# 2          f4a_drh_blood 5.930870938
# 3               f4b_muac 4.974277244
# 4               f4b_resp 4.827065538
# 5               f4b_temp 4.414632255
# 6          f4a_ppl_house 4.077667949
# 7           f4a_drh_days 3.292963455
# 8          f4a_slp_rooms 2.330792448
# 9          f4a_share_fac 2.087144861
# 10        f4a_offr_drink 1.931412712
# 11      f4a_yng_children 1.900575049
# 12         f4a_prim_schl 1.882336592
# 13       f4a_water_avail 1.680083777
# 14        f4a_drh_strain 1.603453863
# 15         f4a_breastfed 1.561997649
# 16                  site 1.490704701
# 17        f4a_max_stools 1.472628277
# 18          f4a_ms_water 1.425606197
# 19         f4a_fac_waste 1.415493079
# 20         f4b_recommend 1.354207737
# 21              f4b_eyes 1.291458032
# 22        f4a_trt_method 1.189874242
# 23       f4a_cur_thirsty 1.130419516
# 24     f4a_hometrt_maize 1.109758698
# 25          f4a_wash_use 1.068797414
# 26        f4a_drh_thirst 1.068300832
# 27       f4a_fuel_natgas 1.067476659
# 28        f4a_disp_feces 1.049399956
# 29          f4a_dad_live 1.006004990
# 30     f4a_drh_bellypain 0.996619009
# 31          f4a_wash_eat 0.991938601
# 32       f4a_hometrt_ors 0.946680498
# 33        f4a_seek_pharm 0.917117278
# 34            f4b_mental 0.915828710
# 35       f4a_house_phone 0.907985648
# 36         f4a_drh_cough 0.905836934
# 37      f4a_hometrt_milk 0.897480279
# 38      f4a_cur_drymouth 0.869140855
# 39         f4a_wash_cook 0.866546757
# 40         f4a_drh_vomit 0.866356280
# 41     f4a_hometrt_othr1 0.832934376
# 42       f4a_store_water 0.807782072
# 43        f4a_wash_child 0.806958369
# 44          f4a_wash_def 0.797526761
# 45      f4a_seek_outside 0.792463756
# 46             f3_gender 0.785385921
# 47        f4a_wash_nurse 0.778273650
# 48         f4a_trt_water 0.777584997
# 49             f4b_mouth 0.766020736
# 50      f4a_house_fridge 0.765913490
# 51      f4a_cur_restless 0.749222073
# 52        f4a_house_tele 0.749064691
# 53        f4a_house_bike 0.740285485
# 54      f4a_hometrt_none 0.737370652
# 55         f4a_fuel_crop 0.734404136
# 56              f4b_skin 0.718727973
# 57      f4a_drh_restless 0.706367969
# 58  f4a_drh_lethrgy_miss 0.705555176
# 59           f4a_ani_cow 0.703401417
# 60           f4a_ani_cat 0.696667717
# 61        f4a_fuel_grass 0.683567956
# 62         f4a_fuel_wood 0.658838239
# 63          f4a_ani_fowl 0.650124680
# 64       f4a_ani_rodents 0.647294168
# 65          f4a_cur_skin 0.647249393
# 66          f4a_ani_goat 0.647072580
# 67        f4a_house_elec 0.635914965
# 68             f4a_floor 0.634885150
# 69            f4a_ani_no 0.629420225
# 70             f4b_admit 0.628165216
# 71       f4a_house_radio 0.618091898
# 72       f4a_water_house 0.617908579
# 73    f4a_cur_fastbreath 0.608105269
# 74           f3_drh_hosp 0.584624736
# 75      f4a_seek_privdoc 0.569946170
# 76           f4a_ani_dog 0.555143674
# 77     f4a_drh_lessdrink 0.553309849
# 78         f4a_wash_othr 0.549932398
# 79      f4a_house_agland 0.539349428
# 80      f4a_water_pubtap 0.528251380
# 81        f4a_fuel_other 0.526459712
# 82        f4a_hometrt_ab 0.518840573
# 83         f3_drh_turgor 0.518164114
# 84      f4a_water_bought 0.516265122
# 85        f4a_house_boat 0.514378772
# 86         f4a_house_car 0.509281940
# 87      f4a_hometrt_zinc 0.506543628
# 88      f4a_drh_prolapse 0.503669117
# 89        f4a_water_yard 0.472285480
# 90   f4a_water_shallwell 0.471823671
# 91       f4a_house_scoot 0.466099799
# 92       f4a_wash_animal 0.463561832
# 93         f4a_fuel_dung 0.447262510
# 94        f4a_seek_remdy 0.417999532
# 95        f4a_house_none 0.412484901
# 96      f4a_hometrt_herb 0.406472231
# 97             f3_drh_iv 0.404488272
# 98    f4a_water_deepwell 0.398126572
# 99        f4b_under_nutr 0.397491136
# 100    f4a_fuel_charcoal 0.394149847
# 101         f4b_abn_hair 0.392905251
# 102        f4a_ani_sheep 0.365129533
# 103    f4a_hometrt_othr2 0.352539165
# 104        f4a_ani_other 0.347299898
# 105      f4b_chest_indrw 0.337227046
# 106         f4a_seek_doc 0.319062879
# 107       f4b_skin_flaky 0.313580488
# 108       f4a_water_well 0.304134516
# 109       f4a_house_cart 0.274626000
# 110           f4b_rectal 0.260267883
# 111      f4a_fuel_biogas 0.257293619
# 112    f4a_water_pubwell 0.252660388
# 113        f4a_fuel_kero 0.240024913
# 114        f4a_drh_consc 0.237884091
# 115     f4a_relationship 0.228153873
# 116      f4a_seek_healer 0.227858534
# 117  f4a_hometrt_othrliq 0.193969051
# 118       f4a_water_bore 0.180745971
# 119        f4a_fuel_coal 0.172020896
# 120       f4a_water_othr 0.169012387
# 121          f4b_bipedal 0.168310075
# 122     f4a_fuel_propane 0.162472075
# 123   f4a_water_covpwell 0.160150316
# 124       f4a_water_pond 0.148038704
# 125       f4a_water_rain 0.132076261
# 126      f4a_water_river 0.128481624
# 127   f4a_water_unspring 0.087967353
# 128       f4a_seek_other 0.071260774
# 129    f4a_water_covwell 0.070613221
# 130      f4a_seek_friend 0.011410715
# 131         f4a_drh_conv 0.008879731
# 132  f4a_water_prospring 0.004420541
# 133        f4a_fuel_elec 0.002048995

# AUC          SE     lower     upper level Model nvar
# 1 0.7433524 0.005089053 0.7333781 0.7533268  0.95    LR    2
# 2 0.7780541 0.004483613 0.7692664 0.7868419  0.95    LR    5
# 3 0.7714799 0.004554937 0.7625524 0.7804074  0.95    LR   10
# 4 0.7225813 0.005335343 0.7121242 0.7330384  0.95    RF    2
# 5 0.7183281 0.005123347 0.7082865 0.7283697  0.95    RF    5
# 6 0.7292662 0.005032925 0.7194018 0.7391305  0.95    RF   10

# nvar     intc intc_LCI intc_UCI slope slope_LCI slope_UCI
# <dbl>    <dbl>    <dbl>    <dbl> <dbl>     <dbl>     <dbl>
#   1     2 -0.0106    -0.421    0.365 0.948     0.568      1.35
# 2     5 -0.00784   -0.425    0.375 0.962     0.629      1.32
# 3    10 -0.00496   -0.425    0.381 0.886     0.569      1.23

####################
#shigella AFe>=0.5 12-23mo
################### var screening, AUC ####
names<-append(x=names, values=c("f4b_muac"))
names <- names[!names %in% c("f4b_haz")]

shigella0.5MUAC.1223 <- CPR.funct(data=cases_age2,outcome="shigella_afe0.5",iter=100,nvars_opts=c(2,5,10))
shigella0.5MUAC.1223[["df_imps"]]
shigella0.5MUAC.1223[["AUC_df"]]
shigella0.5MUAC.1223[["calib"]]

# names     var_red
# 1          f4a_drh_blood 32.21659673
# 2               f4b_muac 13.95602274
# 3               f4b_resp 13.91927634
# 4               f4b_temp 12.42099092
# 5               base_age 11.46012625
# 6               f4b_eyes 11.35739955
# 7          f4a_ppl_house 11.05951218
# 8                   site  8.26282432
# 9           f4a_drh_days  7.14848951
# 10        f4a_drh_strain  6.95810014
# 11         f4a_share_fac  6.95806891
# 12         f4a_slp_rooms  6.82630318
# 13         f4a_drh_vomit  6.21888679
# 14      f4a_yng_children  5.60489674
# 15             f4b_mouth  5.37564974
# 16         f4a_prim_schl  5.13145451
# 17       f4a_store_water  4.36563643
# 18         f4b_recommend  4.20406024
# 19        f4a_offr_drink  4.02418950
# 20        f4a_max_stools  4.00764862
# 21        f4a_fuel_grass  3.79308164
# 22         f4a_breastfed  3.70930973
# 23       f4a_water_avail  3.35607293
# 24     f4a_drh_bellypain  3.14070790
# 25        f4a_trt_method  3.06384239
# 26         f4a_fac_waste  3.00219949
# 27            f4b_mental  2.99926304
# 28        f4a_drh_thirst  2.96962563
# 29          f4a_dad_live  2.92501888
# 30      f4a_cur_drymouth  2.90151936
# 31          f4a_ms_water  2.81712190
# 32        f4a_disp_feces  2.78981850
# 33         f4a_drh_cough  2.52197843
# 34  f4a_drh_lethrgy_miss  2.44361594
# 35       f4a_cur_thirsty  2.41128053
# 36       f4a_hometrt_ors  2.37764232
# 37        f4a_wash_child  2.30461473
# 38         f4a_wash_cook  2.27075673
# 39             f4b_admit  2.26996431
# 40          f4a_ani_fowl  2.24354036
# 41              f4b_skin  2.22332853
# 42         f4a_fuel_crop  2.17777504
# 43             f3_gender  2.14709785
# 44        f4a_wash_nurse  2.13734804
# 45          f4a_wash_use  2.13387612
# 46        f4a_house_bike  2.09925135
# 47           f4a_ani_cat  2.09716171
# 48      f4a_hometrt_none  2.09398315
# 49         f4a_fuel_wood  2.08189009
# 50          f4a_wash_def  2.05827005
# 51      f4a_drh_restless  2.05012147
# 52        f4a_house_tele  2.04706974
# 53       f4a_house_scoot  2.02456203
# 54      f4a_cur_restless  2.01693578
# 55          f4a_cur_skin  1.97740114
# 56      f4a_seek_outside  1.94134705
# 57     f4a_fuel_charcoal  1.89181256
# 58           f4a_ani_dog  1.88229540
# 59       f4a_house_phone  1.85432855
# 60      f4a_water_pubtap  1.82199300
# 61       f4a_ani_rodents  1.80665811
# 62       f4a_house_radio  1.80386014
# 63        f4b_under_nutr  1.78713672
# 64    f4a_water_deepwell  1.78111956
# 65          f4a_ani_goat  1.74316410
# 66         f4a_fuel_dung  1.74085263
# 67    f4a_cur_fastbreath  1.72948839
# 68             f3_drh_iv  1.72806132
# 69           f3_drh_hosp  1.69683725
# 70           f4a_ani_cow  1.67928715
# 71        f4a_house_elec  1.67702260
# 72          f4a_wash_eat  1.67547795
# 73         f4a_trt_water  1.65750182
# 74         f3_drh_turgor  1.64031519
# 75             f4a_floor  1.57076604
# 76      f4a_house_fridge  1.54503915
# 77      f4a_relationship  1.49835033
# 78     f4a_hometrt_maize  1.46395140
# 79      f4a_house_agland  1.43795127
# 80     f4a_drh_lessdrink  1.41365248
# 81     f4a_hometrt_othr1  1.37282913
# 82            f4a_ani_no  1.37046174
# 83       f4a_fuel_natgas  1.35407898
# 84       f4a_seek_healer  1.30790355
# 85         f4a_ani_sheep  1.29335608
# 86       f4a_wash_animal  1.25054048
# 87         f4a_ani_other  1.24434909
# 88      f4a_hometrt_herb  1.18901436
# 89         f4a_fuel_kero  1.16854211
# 90        f4a_seek_pharm  1.15956557
# 91        f4a_hometrt_ab  1.15291152
# 92        f4a_water_yard  1.14480500
# 93        f4a_house_cart  1.13060536
# 94       f4a_water_house  1.12322677
# 95         f4a_wash_othr  1.09056967
# 96         f4a_house_car  1.08742573
# 97     f4a_water_pubwell  1.07021832
# 98        f4a_water_bore  0.97296415
# 99      f4a_water_bought  0.96544908
# 100    f4a_water_covwell  0.81430338
# 101  f4a_water_shallwell  0.80827196
# 102       f4a_water_rain  0.73558267
# 103     f4a_hometrt_zinc  0.72246948
# 104       f4a_house_none  0.70819324
# 105         f4b_abn_hair  0.69455358
# 106     f4a_seek_privdoc  0.67854413
# 107      f4b_chest_indrw  0.64954793
# 108       f4a_fuel_other  0.64103465
# 109     f4a_fuel_propane  0.61982258
# 110         f4a_drh_conv  0.58276273
# 111       f4b_skin_flaky  0.56393550
# 112       f4a_water_well  0.49861792
# 113       f4a_seek_remdy  0.47934489
# 114        f4a_fuel_coal  0.47464014
# 115    f4a_hometrt_othr2  0.46262239
# 116         f4a_seek_doc  0.45643766
# 117      f4a_water_river  0.41315584
# 118        f4a_drh_consc  0.39684530
# 119       f4a_seek_other  0.37608795
# 120        f4a_fuel_elec  0.31749713
# 121  f4a_hometrt_othrliq  0.30323353
# 122      f4a_fuel_biogas  0.26724846
# 123   f4a_water_covpwell  0.26316857
# 124       f4a_water_othr  0.25555258
# 125     f4a_drh_prolapse  0.23950756
# 126          f4b_bipedal  0.21239664
# 127       f4a_house_boat  0.18352019
# 128     f4a_hometrt_milk  0.16588997
# 129       f4a_water_pond  0.16398975
# 130      f4a_seek_friend  0.12739020
# 131   f4a_water_unspring  0.12635516
# 132           f4b_rectal  0.08736565
# 133  f4a_water_prospring  0.01108000

# AUC          SE     lower     upper level Model nvar
# 1 0.7216550 0.003102217 0.7155747 0.7277352  0.95    LR    2
# 2 0.7346508 0.002997997 0.7287749 0.7405268  0.95    LR    5
# 3 0.7642473 0.002782977 0.7587927 0.7697018  0.95    LR   10
# 4 0.7178200 0.003107149 0.7117301 0.7239099  0.95    RF    2
# 5 0.7367445 0.002928728 0.7310043 0.7424847  0.95    RF    5
# 6 0.7597467 0.002798965 0.7542609 0.7652326  0.95    RF   10

# nvar     intc intc_LCI intc_UCI slope slope_LCI slope_UCI
# <dbl>    <dbl>    <dbl>    <dbl> <dbl>     <dbl>     <dbl>
#   1     2 -0.00590   -0.261    0.243 0.985     0.753      1.23
# 2     5 -0.00630   -0.263    0.244 0.971     0.743      1.21
# 3    10 -0.00619   -0.267    0.249 0.937     0.721      1.17

####################
#shigella AFe>=0.5 0-23mo
################### var screening, AUC ####
names<-append(x=names, values=c("f4b_muac"))
names <- names[!names %in% c("f4b_haz")]

shigella0.5MUAC.023 <- CPR.funct(data=cases_age4,outcome="shigella_afe0.5",iter=100,nvars_opts=c(2,5,10))
shigella0.5MUAC.023[["df_imps"]]
shigella0.5MUAC.023[["AUC_df"]]
shigella0.5MUAC.023[["calib"]]

# names     var_red
# 1               base_age 58.62760983
# 2          f4a_drh_blood 42.34691644
# 3               f4b_resp 21.20085135
# 4               f4b_muac 20.39777531
# 5               f4b_temp 17.94649443
# 6          f4a_ppl_house 15.64660340
# 7               f4b_eyes 11.17708525
# 8          f4a_breastfed 10.74931269
# 9           f4a_drh_days 10.68773754
# 10         f4a_slp_rooms  9.77928162
# 11                  site  9.01339170
# 12         f4a_share_fac  8.92239198
# 13      f4a_yng_children  7.72530541
# 14         f4a_prim_schl  7.56203840
# 15        f4a_drh_strain  7.49885884
# 16         f4a_drh_vomit  7.07472039
# 17        f4a_max_stools  6.16676330
# 18        f4a_offr_drink  5.94462482
# 19         f4b_recommend  5.92667987
# 20         f4a_fac_waste  4.92263940
# 21             f4b_mouth  4.85121637
# 22       f4a_water_avail  4.60820200
# 23          f4a_dad_live  4.47447670
# 24     f4a_drh_bellypain  4.26483857
# 25          f4a_ms_water  4.24084452
# 26        f4a_disp_feces  4.23868974
# 27        f4a_fuel_grass  4.17728072
# 28        f4a_trt_method  4.14958236
# 29         f4a_drh_cough  3.91182605
# 30      f4a_cur_drymouth  3.78278312
# 31        f4a_drh_thirst  3.72201526
# 32       f4a_hometrt_ors  3.41018650
# 33        f4a_wash_child  3.40037903
# 34             f3_gender  3.39451872
# 35       f4a_cur_thirsty  3.36067215
# 36       f4a_store_water  3.32818336
# 37          f4a_wash_use  3.29729137
# 38         f4a_wash_cook  3.23144490
# 39      f4a_hometrt_none  3.23072543
# 40              f4b_skin  3.16610418
# 41            f4b_mental  3.16073579
# 42         f4a_fuel_crop  3.13218536
# 43        f4a_house_bike  3.08408355
# 44        f4a_wash_nurse  3.06867227
# 45  f4a_drh_lethrgy_miss  3.04727415
# 46      f4a_seek_outside  3.00483999
# 47          f4a_ani_fowl  2.99630325
# 48      f4a_cur_restless  2.97914474
# 49      f4a_drh_restless  2.96724721
# 50          f4a_wash_def  2.96125415
# 51        f4a_house_tele  2.94902411
# 52         f4a_fuel_wood  2.90472400
# 53       f4a_house_phone  2.83175207
# 54             f4b_admit  2.80227130
# 55             f4a_floor  2.79719630
# 56          f4a_wash_eat  2.77642625
# 57           f4a_ani_cat  2.74019705
# 58       f4a_house_scoot  2.73620118
# 59          f4a_cur_skin  2.71086896
# 60           f4a_ani_dog  2.65434652
# 61    f4a_cur_fastbreath  2.62210695
# 62       f4a_house_radio  2.61673673
# 63      f4a_house_agland  2.61146181
# 64       f4a_ani_rodents  2.59443473
# 65           f3_drh_hosp  2.54434498
# 66           f4a_ani_cow  2.53132981
# 67      f4a_house_fridge  2.45947211
# 68       f4a_fuel_natgas  2.45835383
# 69        f4a_house_elec  2.45077269
# 70          f4a_ani_goat  2.43201401
# 71     f4a_fuel_charcoal  2.42252340
# 72             f3_drh_iv  2.38890295
# 73         f4a_trt_water  2.37511916
# 74     f4a_hometrt_othr1  2.34387734
# 75        f4b_under_nutr  2.32501445
# 76     f4a_hometrt_maize  2.30395035
# 77        f4a_seek_pharm  2.29699719
# 78         f3_drh_turgor  2.22822408
# 79      f4a_water_pubtap  2.22414769
# 80      f4a_relationship  2.10543528
# 81     f4a_drh_lessdrink  2.02066572
# 82    f4a_water_deepwell  2.01653821
# 83         f4a_wash_othr  1.99624645
# 84       f4a_wash_animal  1.95903140
# 85            f4a_ani_no  1.94811232
# 86        f4a_hometrt_ab  1.89259641
# 87         f4a_fuel_dung  1.81064363
# 88      f4a_hometrt_herb  1.73991771
# 89       f4a_water_house  1.72554886
# 90         f4a_ani_other  1.71922205
# 91         f4a_ani_sheep  1.67194257
# 92         f4a_house_car  1.53652649
# 93        f4a_water_yard  1.51186518
# 94        f4a_house_cart  1.48682837
# 95      f4a_water_bought  1.46752918
# 96      f4a_hometrt_zinc  1.32421288
# 97       f4a_seek_healer  1.31738232
# 98         f4a_fuel_kero  1.30914612
# 99     f4a_water_pubwell  1.28329081
# 100  f4a_water_shallwell  1.26051730
# 101         f4b_abn_hair  1.16362664
# 102     f4a_seek_privdoc  1.15082940
# 103       f4a_house_none  1.07671591
# 104       f4a_water_bore  1.06554289
# 105       f4a_seek_remdy  1.00719558
# 106     f4a_hometrt_milk  0.98636027
# 107       f4a_fuel_other  0.94595453
# 108         f4a_seek_doc  0.93541513
# 109         f4a_drh_conv  0.91908083
# 110      f4b_chest_indrw  0.88413111
# 111     f4a_fuel_propane  0.86343584
# 112    f4a_water_covwell  0.84660340
# 113       f4a_water_well  0.83528171
# 114    f4a_hometrt_othr2  0.83449694
# 115       f4a_water_rain  0.80305463
# 116       f4b_skin_flaky  0.77818889
# 117        f4a_drh_consc  0.71886096
# 118     f4a_drh_prolapse  0.69916965
# 119        f4a_fuel_coal  0.63013736
# 120       f4a_house_boat  0.57599572
# 121  f4a_hometrt_othrliq  0.51981292
# 122       f4a_water_othr  0.51393396
# 123      f4a_water_river  0.49412405
# 124      f4a_fuel_biogas  0.47575032
# 125   f4a_water_covpwell  0.45696884
# 126       f4a_seek_other  0.45144645
# 127        f4a_fuel_elec  0.42278766
# 128           f4b_rectal  0.35331843
# 129          f4b_bipedal  0.34757069
# 130       f4a_water_pond  0.29536788
# 131      f4a_seek_friend  0.23181194
# 132   f4a_water_unspring  0.14555465
# 133  f4a_water_prospring  0.03086628

# AUC          SE     lower     upper level Model nvar
# 1 0.8199190 0.001986181 0.8160261 0.8238118  0.95    LR    2
# 2 0.8216321 0.001901198 0.8179058 0.8253584  0.95    LR    5
# 3 0.8266144 0.001870984 0.8229473 0.8302815  0.95    LR   10
# 4 0.8127537 0.002038897 0.8087575 0.8167499  0.95    RF    2
# 5 0.8092479 0.002004396 0.8053193 0.8131764  0.95    RF    5
# 6 0.8301809 0.001869511 0.8265168 0.8338451  0.95    RF   10

# nvar      intc intc_LCI intc_UCI slope slope_LCI slope_UCI
# <dbl>     <dbl>    <dbl>    <dbl> <dbl>     <dbl>     <dbl>
#   1     2 -0.00187    -0.214    0.204 1.02      0.847      1.20
# 2     5 -0.00151    -0.214    0.206 1.01      0.843      1.19
# 3    10 -0.000706   -0.215    0.208 0.999     0.834      1.18

####################
#shigella AFe>=0.5 24-59mo
################### var screening, AUC ####
names<-append(x=names, values=c("f4b_muac"))
names <- names[!names %in% c("f4b_haz")]

shigella0.5MUAC.2459 <- CPR.funct(data=cases_age3,outcome="shigella_afe0.5",iter=100,nvars_opts=c(2,5,10))
shigella0.5MUAC.2459[["df_imps"]]
shigella0.5MUAC.2459[["AUC_df"]]
shigella0.5MUAC.2459[["calib"]]

# names      var_red
# 1          f4a_drh_blood 28.944513787
# 2               f4b_eyes 16.247438586
# 3                   site 14.730750958
# 4               f4b_muac 10.066631508
# 5               f4b_temp 10.027035657
# 6               f4b_resp  9.502113757
# 7         f4a_fuel_grass  8.507841107
# 8               base_age  8.434136053
# 9          f4a_ppl_house  7.845541813
# 10       f4a_store_water  6.999218978
# 11         f4a_fuel_crop  6.263540819
# 12        f4a_drh_strain  5.896632971
# 13             f4b_mouth  5.834510227
# 14         f4a_share_fac  4.899782360
# 15         f4a_prim_schl  4.899467862
# 16          f4a_drh_days  4.805611511
# 17         f4a_slp_rooms  4.465885833
# 18      f4a_yng_children  4.326674653
# 19       f4a_cur_thirsty  4.017523048
# 20         f4a_fuel_dung  3.860297654
# 21         f4a_drh_vomit  3.589483641
# 22        f4a_drh_thirst  3.504496978
# 23         f4b_recommend  3.091962083
# 24        f4a_offr_drink  3.064893204
# 25        f4a_max_stools  2.946762996
# 26         f4a_fac_waste  2.860041181
# 27        f4a_trt_method  2.802054875
# 28        f4a_disp_feces  2.605561726
# 29          f4a_ms_water  2.553855825
# 30            f4b_mental  2.463717175
# 31         f4a_breastfed  2.430045977
# 32       f4a_water_avail  2.340265981
# 33         f4a_drh_cough  2.248495091
# 34    f4a_water_deepwell  2.242947648
# 35        f4a_seek_pharm  2.123133486
# 36     f4a_drh_bellypain  2.119355994
# 37      f4a_cur_drymouth  2.089679181
# 38      f4a_seek_outside  2.081884043
# 39          f4a_wash_eat  2.044353173
# 40        f4a_wash_nurse  1.982519237
# 41      f4a_relationship  1.947097097
# 42     f4a_fuel_charcoal  1.912632480
# 43      f4a_hometrt_none  1.908044251
# 44          f4a_dad_live  1.907192995
# 45       f4a_hometrt_ors  1.861689267
# 46             f3_gender  1.824859559
# 47   f4a_water_shallwell  1.806960538
# 48      f4a_water_pubtap  1.749812016
# 49           f4a_ani_cat  1.733681994
# 50  f4a_drh_lethrgy_miss  1.719072513
# 51        f4a_house_tele  1.694685812
# 52       f4a_house_phone  1.602291617
# 53         f4a_wash_cook  1.592574471
# 54              f4b_skin  1.584039700
# 55          f4a_cur_skin  1.580127776
# 56          f4a_wash_def  1.575707789
# 57        f4a_house_bike  1.545825311
# 58             f4a_floor  1.540411386
# 59          f4a_wash_use  1.535812639
# 60        f4a_house_elec  1.528546890
# 61       f4a_house_radio  1.496421967
# 62       f4a_ani_rodents  1.446309698
# 63        f4a_wash_child  1.422217109
# 64           f4a_ani_dog  1.404759447
# 65          f4a_ani_fowl  1.388599845
# 66      f4a_drh_restless  1.353251189
# 67             f4b_admit  1.322861497
# 68         f4a_trt_water  1.315612001
# 69         f4a_fuel_wood  1.299161667
# 70      f4a_cur_restless  1.295595440
# 71     f4a_hometrt_maize  1.253736941
# 72        f4a_hometrt_ab  1.251323265
# 73         f3_drh_turgor  1.241740014
# 74     f4a_drh_lessdrink  1.226857465
# 75             f3_drh_iv  1.218798184
# 76          f4a_ani_goat  1.198839302
# 77      f4a_house_agland  1.184402760
# 78      f4a_house_fridge  1.160079858
# 79           f3_drh_hosp  1.159430909
# 80          f4a_drh_conv  1.065148165
# 81           f4a_ani_cow  1.057564393
# 82       f4a_house_scoot  1.056529419
# 83         f4a_ani_sheep  1.044809922
# 84         f4a_fuel_kero  1.021878819
# 85            f4a_ani_no  1.002141551
# 86     f4a_hometrt_othr1  0.967892731
# 87        f4a_water_yard  0.964365484
# 88        f4b_under_nutr  0.959556237
# 89         f4a_fuel_elec  0.935863022
# 90      f4a_water_bought  0.746094099
# 91        f4a_house_cart  0.745391390
# 92       f4a_water_house  0.743799827
# 93         f4a_house_car  0.700328506
# 94      f4a_hometrt_herb  0.692876129
# 95         f4a_ani_other  0.671758350
# 96      f4a_fuel_propane  0.664123726
# 97    f4a_cur_fastbreath  0.654648791
# 98         f4a_wash_othr  0.627273928
# 99       f4a_fuel_natgas  0.624121315
# 100        f4a_fuel_coal  0.575489561
# 101      f4a_wash_animal  0.533081910
# 102       f4a_water_well  0.523577903
# 103       f4a_water_bore  0.495033023
# 104    f4a_water_pubwell  0.465959307
# 105       f4a_house_none  0.445938804
# 106       f4a_water_rain  0.442960429
# 107         f4a_seek_doc  0.418104550
# 108   f4a_water_covpwell  0.406014104
# 109      f4a_seek_healer  0.402046196
# 110         f4b_abn_hair  0.383616267
# 111     f4a_hometrt_zinc  0.380110903
# 112       f4a_seek_other  0.334497606
# 113       f4a_seek_remdy  0.331760201
# 114     f4a_seek_privdoc  0.312501107
# 115     f4a_drh_prolapse  0.304204017
# 116        f4a_drh_consc  0.299061409
# 117       f4a_water_othr  0.294035910
# 118       f4a_water_pond  0.287173184
# 119    f4a_hometrt_othr2  0.283584118
# 120       f4a_fuel_other  0.271138216
# 121      f4a_water_river  0.243593671
# 122    f4a_water_covwell  0.233759513
# 123       f4b_skin_flaky  0.202772216
# 124      f4a_fuel_biogas  0.195338288
# 125     f4a_hometrt_milk  0.194626603
# 126       f4a_house_boat  0.172256479
# 127  f4a_hometrt_othrliq  0.130681606
# 128   f4a_water_unspring  0.129542114
# 129           f4b_rectal  0.122178314
# 130          f4b_bipedal  0.118897413
# 131  f4a_water_prospring  0.100593860
# 132      f4b_chest_indrw  0.060968919
# 133      f4a_seek_friend  0.008535635

# AUC          SE     lower     upper level Model nvar
# 1 0.7999075 0.004315140 0.7914500 0.8083650  0.95    LR    2
# 2 0.8211246 0.002790859 0.8156547 0.8265946  0.95    LR    5
# 3 0.8193452 0.002812279 0.8138332 0.8248571  0.95    LR   10
# 4 0.8002734 0.004311913 0.7918222 0.8087246  0.95    RF    2
# 5 0.8173376 0.002757242 0.8119335 0.8227417  0.95    RF    5
# 6 0.8164928 0.002818351 0.8109689 0.8220166  0.95    RF   10

# nvar     intc intc_LCI intc_UCI slope slope_LCI slope_UCI
# <dbl>    <dbl>    <dbl>    <dbl> <dbl>     <dbl>     <dbl>
#   1     2 -0.00314   -0.300    0.287 1.00      0.786      1.24
# 2     5 -0.00115   -0.309    0.299 0.980     0.768      1.22
# 3    10 -0.00267   -0.311    0.298 0.963     0.753      1.20


####################
#shigella AFe>=0.3 only MUAC
################### var screening, AUC ####
names<-append(x=names, values=c("f4b_muac"))
names <- names[!names %in% c("f4b_haz")]

shigella0.3 <- CPR.funct(data=cases,outcome="shigella_afe0.3",iter=100,nvars_opts=c(2,5,10))
shigella0.3[["df_imps"]]
shigella0.3[["AUC_df"]]
shigella0.3[["calib"]]

# names     var_red
# 1               base_age 101.7168703
# 2          f4a_drh_blood  75.4269335
# 3               f4b_muac  38.7763062
# 4               f4b_resp  35.3351048
# 5               f4b_temp  32.9673697
# 6               f4b_eyes  26.5883921
# 7          f4a_ppl_house  25.9720623
# 8           f4a_drh_days  18.5523075
# 9                   site  17.7637739
# 10         f4a_slp_rooms  16.6789352
# 11         f4a_share_fac  15.2390545
# 12         f4a_breastfed  14.5262194
# 13      f4a_yng_children  14.2396443
# 14         f4a_prim_schl  13.2983031
# 15        f4a_drh_strain  12.5081911
# 16         f4a_drh_vomit  11.9042131
# 17             f4b_mouth  11.1161514
# 18        f4a_max_stools  10.0595242
# 19        f4a_offr_drink   9.8473339
# 20         f4b_recommend   9.7071019
# 21        f4a_fuel_grass   8.9943136
# 22         f4a_fac_waste   8.5005942
# 23     f4a_drh_bellypain   8.3291072
# 24         f4a_drh_cough   7.7390943
# 25          f4a_ms_water   7.5847943
# 26        f4a_disp_feces   7.5169662
# 27        f4a_trt_method   7.4868682
# 28       f4a_water_avail   7.4171774
# 29          f4a_dad_live   7.3375862
# 30        f4a_drh_thirst   6.9286375
# 31       f4a_store_water   6.9217170
# 32         f4a_fuel_crop   6.7890778
# 33       f4a_cur_thirsty   6.7548472
# 34            f4b_mental   6.5348961
# 35             f3_gender   6.0132738
# 36      f4a_hometrt_none   5.8770189
# 37         f4a_wash_cook   5.7748680
# 38      f4a_cur_drymouth   5.7499342
# 39          f4a_wash_use   5.6947043
# 40  f4a_drh_lethrgy_miss   5.4906472
# 41        f4a_wash_child   5.4648023
# 42        f4a_wash_nurse   5.4511910
# 43        f4a_house_tele   5.4444060
# 44       f4a_hometrt_ors   5.3044028
# 45        f4a_house_bike   5.2599137
# 46          f4a_cur_skin   5.2220867
# 47      f4a_relationship   5.2110530
# 48          f4a_wash_def   5.1115965
# 49          f4a_ani_fowl   5.0514228
# 50              f4b_skin   5.0381578
# 51      f4a_seek_outside   5.0019552
# 52           f4a_ani_cat   4.9140805
# 53       f4a_house_phone   4.9049396
# 54       f4a_house_radio   4.8620916
# 55     f4a_fuel_charcoal   4.8085588
# 56      f4a_drh_restless   4.7999567
# 57      f4a_cur_restless   4.7512191
# 58         f4a_fuel_wood   4.6395641
# 59          f4a_wash_eat   4.6019624
# 60           f4a_ani_dog   4.5610354
# 61             f4a_floor   4.4521565
# 62             f4b_admit   4.3640854
# 63       f4a_ani_rodents   4.3525186
# 64       f4a_house_scoot   4.3098676
# 65         f4a_fuel_dung   4.2906776
# 66        f4a_house_elec   4.2869309
# 67         f4a_trt_water   4.2755648
# 68           f3_drh_hosp   4.2529971
# 69        f4b_under_nutr   4.2367334
# 70          f4a_ani_goat   4.2268882
# 71      f4a_water_pubtap   4.2066612
# 72    f4a_water_deepwell   4.0164433
# 73    f4a_cur_fastbreath   3.9809467
# 74     f4a_drh_lessdrink   3.9672211
# 75           f4a_ani_cow   3.9670564
# 76      f4a_house_agland   3.8850036
# 77      f4a_house_fridge   3.8708834
# 78         f3_drh_turgor   3.8697131
# 79        f4a_seek_pharm   3.8487228
# 80     f4a_hometrt_maize   3.8395457
# 81     f4a_hometrt_othr1   3.8362756
# 82             f3_drh_iv   3.7291781
# 83        f4a_hometrt_ab   3.4643201
# 84            f4a_ani_no   3.3440218
# 85         f4a_ani_sheep   3.1523738
# 86       f4a_wash_animal   3.0380525
# 87       f4a_fuel_natgas   3.0127811
# 88        f4a_water_yard   2.8179141
# 89      f4a_hometrt_herb   2.8171894
# 90   f4a_water_shallwell   2.8066023
# 91         f4a_house_car   2.7875930
# 92         f4a_wash_othr   2.7501549
# 93        f4a_house_cart   2.6606709
# 94       f4a_water_house   2.6364212
# 95         f4a_ani_other   2.4815257
# 96      f4a_water_bought   2.2765283
# 97          f4a_drh_conv   2.2502724
# 98         f4a_fuel_kero   2.1460443
# 99     f4a_water_pubwell   2.0982902
# 100      f4a_seek_healer   1.8553971
# 101     f4a_fuel_propane   1.7927007
# 102     f4a_hometrt_zinc   1.7224112
# 103       f4a_water_bore   1.6774496
# 104     f4a_seek_privdoc   1.6427591
# 105       f4a_house_none   1.5905114
# 106       f4a_seek_remdy   1.5039739
# 107         f4b_abn_hair   1.5013035
# 108       f4a_water_well   1.4937163
# 109       f4a_water_rain   1.4712073
# 110        f4a_fuel_coal   1.4410233
# 111      f4b_chest_indrw   1.3108850
# 112    f4a_hometrt_othr2   1.3051618
# 113         f4a_seek_doc   1.2971851
# 114       f4a_fuel_other   1.2567057
# 115       f4b_skin_flaky   1.2380002
# 116    f4a_water_covwell   1.2070187
# 117        f4a_drh_consc   1.1952860
# 118     f4a_drh_prolapse   1.1379648
# 119      f4a_water_river   1.0673503
# 120     f4a_hometrt_milk   0.9890136
# 121        f4a_fuel_elec   0.8949267
# 122   f4a_water_covpwell   0.8543794
# 123       f4a_water_othr   0.8435466
# 124      f4a_fuel_biogas   0.7939326
# 125       f4a_water_pond   0.7671402
# 126  f4a_hometrt_othrliq   0.7252243
# 127           f4b_rectal   0.6994717
# 128       f4a_house_boat   0.6895300
# 129       f4a_seek_other   0.6700535
# 130      f4a_seek_friend   0.5025518
# 131          f4b_bipedal   0.4927638
# 132   f4a_water_unspring   0.4301730
# 133  f4a_water_prospring   0.1810077

# AUC          SE     lower     upper level Model nvar
# 1 0.7890271 0.001573267 0.7859436 0.7921107  0.95    LR    2
# 2 0.7815780 0.001598544 0.7784449 0.7847111  0.95    LR    5
# 3 0.7816230 0.001624505 0.7784390 0.7848070  0.95    LR   10
# 4 0.7915991 0.001578589 0.7885051 0.7946931  0.95    RF    2
# 5 0.7911628 0.001567637 0.7880903 0.7942353  0.95    RF    5
# 6 0.8175435 0.001473719 0.8146550 0.8204319  0.95    RF   10

# nvar    intc intc_LCI intc_UCI slope slope_LCI slope_UCI
# <dbl>   <dbl>    <dbl>    <dbl> <dbl>     <dbl>     <dbl>
#   1     2 -0.0138   -0.169    0.139 1.00      0.859      1.15
# 2     5 -0.0143   -0.170    0.139 0.995     0.853      1.14
# 3    10 -0.0154   -0.172    0.139 0.982     0.843      1.13


####################
#shigella AFe>=0.7 only MUAC
################### var screening, AUC ####
names<-append(x=names, values=c("f4b_muac"))
names <- names[!names %in% c("f4b_haz")]

shigella0.7 <- CPR.funct(data=cases,outcome="shigella_afe0.7",iter=100,nvars_opts=c(2,5,10))
shigella0.7[["df_imps"]]
shigella0.7[["AUC_df"]]
shigella0.7[["calib"]]

# names     var_red
# 1               base_age 74.79350765
# 2          f4a_drh_blood 65.37540899
# 3               f4b_eyes 28.97197800
# 4               f4b_muac 28.87074158
# 5               f4b_resp 25.48608734
# 6               f4b_temp 25.11426360
# 7                   site 22.76536284
# 8          f4a_ppl_house 19.95543775
# 9         f4a_drh_strain 14.49290906
# 10        f4a_fuel_grass 14.16975119
# 11          f4a_drh_days 13.32618626
# 12         f4a_slp_rooms 12.42111271
# 13         f4a_share_fac 10.96545850
# 14         f4a_fuel_crop 10.52052891
# 15      f4a_yng_children 10.26434610
# 16         f4a_prim_schl  9.91439235
# 17             f4b_mouth  9.87612462
# 18       f4a_store_water  9.39928577
# 19         f4a_drh_vomit  8.60511387
# 20         f4a_breastfed  8.58095042
# 21         f4a_fac_waste  8.40688954
# 22        f4a_max_stools  8.34577297
# 23         f4b_recommend  8.19517813
# 24        f4a_offr_drink  7.65233594
# 25       f4a_cur_thirsty  6.32889561
# 26       f4a_water_avail  6.32032242
# 27         f4a_drh_cough  6.30718747
# 28     f4a_drh_bellypain  6.08957151
# 29        f4a_trt_method  6.07322453
# 30        f4a_disp_feces  6.07049729
# 31          f4a_dad_live  5.84804644
# 32        f4a_drh_thirst  5.39802328
# 33          f4a_ms_water  5.24508540
# 34    f4a_water_deepwell  4.96955859
# 35         f4a_fuel_dung  4.84770681
# 36            f4b_mental  4.78674221
# 37      f4a_cur_drymouth  4.74759445
# 38       f4a_hometrt_ors  4.60282653
# 39      f4a_hometrt_none  4.50419946
# 40             f3_gender  4.50194732
# 41       f4a_house_radio  4.32457156
# 42          f4a_wash_use  4.26807783
# 43         f4a_wash_cook  4.15392913
# 44      f4a_seek_outside  4.14127381
# 45  f4a_drh_lethrgy_miss  4.12272253
# 46        f4a_wash_child  4.08634173
# 47        f4a_house_tele  4.05964478
# 48              f4b_skin  4.03542378
# 49          f4a_wash_def  3.96862151
# 50        f4a_wash_nurse  3.96775206
# 51     f4a_fuel_charcoal  3.95410795
# 52          f4a_ani_fowl  3.91908005
# 53         f4a_fuel_wood  3.90143907
# 54        f4a_house_bike  3.90026579
# 55             f4b_admit  3.84272148
# 56      f4a_drh_restless  3.78424182
# 57   f4a_water_shallwell  3.76289074
# 58          f4a_wash_eat  3.73214936
# 59      f4a_cur_restless  3.71961580
# 60             f4a_floor  3.67387394
# 61           f3_drh_hosp  3.64643859
# 62      f4a_relationship  3.55902575
# 63       f4a_house_phone  3.50997462
# 64           f4a_ani_cat  3.48621656
# 65        f4a_seek_pharm  3.44173800
# 66       f4a_house_scoot  3.43501827
# 67           f4a_ani_dog  3.40273524
# 68      f4a_house_fridge  3.33485335
# 69        f4a_house_elec  3.32533243
# 70      f4a_water_pubtap  3.27887912
# 71          f4a_cur_skin  3.27780056
# 72      f4a_house_agland  3.26570361
# 73       f4a_ani_rodents  3.17047849
# 74     f4a_drh_lessdrink  3.13473516
# 75             f3_drh_iv  3.02409222
# 76           f4a_ani_cow  2.98131584
# 77         f3_drh_turgor  2.97832985
# 78         f4a_trt_water  2.97508087
# 79    f4a_cur_fastbreath  2.94269020
# 80          f4a_ani_goat  2.90176054
# 81     f4a_hometrt_othr1  2.82961433
# 82        f4a_hometrt_ab  2.75026729
# 83        f4b_under_nutr  2.70126251
# 84            f4a_ani_no  2.60264959
# 85          f4a_drh_conv  2.42525335
# 86       f4a_fuel_natgas  2.35249591
# 87     f4a_hometrt_maize  2.30128885
# 88       f4a_wash_animal  2.23125441
# 89         f4a_ani_sheep  2.17530176
# 90        f4a_water_yard  2.10077490
# 91         f4a_wash_othr  2.08126254
# 92         f4a_house_car  1.97545806
# 93      f4a_hometrt_herb  1.97249862
# 94        f4a_house_cart  1.88993798
# 95         f4a_ani_other  1.81744880
# 96       f4a_water_house  1.77907797
# 97      f4a_hometrt_zinc  1.74293936
# 98      f4a_water_bought  1.64842970
# 99      f4a_fuel_propane  1.62077129
# 100    f4a_water_pubwell  1.61169110
# 101        f4a_fuel_kero  1.38831321
# 102      f4a_seek_healer  1.35257500
# 103       f4a_seek_remdy  1.29989557
# 104         f4a_seek_doc  1.24060978
# 105     f4a_seek_privdoc  1.21324651
# 106       f4a_water_bore  1.20238519
# 107       f4a_water_well  1.13951629
# 108         f4b_abn_hair  1.13585638
# 109     f4a_drh_prolapse  1.05419858
# 110        f4a_drh_consc  1.02556775
# 111    f4a_hometrt_othr2  0.97831895
# 112       f4a_water_rain  0.85682160
# 113      f4b_chest_indrw  0.85625153
# 114   f4a_water_covpwell  0.84504856
# 115       f4b_skin_flaky  0.84215845
# 116        f4a_fuel_coal  0.81693420
# 117       f4a_fuel_other  0.80626775
# 118       f4a_house_none  0.79870438
# 119    f4a_water_covwell  0.70778672
# 120       f4a_house_boat  0.66818445
# 121       f4a_water_othr  0.59368417
# 122          f4b_bipedal  0.57992015
# 123  f4a_hometrt_othrliq  0.56203451
# 124      f4a_water_river  0.54128256
# 125      f4a_fuel_biogas  0.52166636
# 126        f4a_fuel_elec  0.51903747
# 127       f4a_seek_other  0.51688391
# 128       f4a_water_pond  0.48245345
# 129           f4b_rectal  0.45152973
# 130   f4a_water_unspring  0.33543575
# 131     f4a_hometrt_milk  0.27369604
# 132      f4a_seek_friend  0.06273676
# 133  f4a_water_prospring  0.05915352

# AUC          SE     lower     upper level Model nvar
# 1 0.8218085 0.001631988 0.8186099 0.8250071  0.95    LR    2
# 2 0.8235018 0.001651319 0.8202652 0.8267383  0.95    LR    5
# 3 0.8244628 0.001698164 0.8211344 0.8277911  0.95    LR   10
# 4 0.8251692 0.001636312 0.8219620 0.8283763  0.95    RF    2
# 5 0.8411276 0.001583031 0.8380249 0.8442303  0.95    RF    5
# 6 0.8549177 0.001519202 0.8519401 0.8578953  0.95    RF   10

# nvar     intc intc_LCI intc_UCI slope slope_LCI slope_UCI
# <dbl>    <dbl>    <dbl>    <dbl> <dbl>     <dbl>     <dbl>
#   1     2 -0.00638   -0.183    0.166 1.00      0.864      1.15
# 2     5 -0.00571   -0.185    0.170 0.997     0.863      1.14
# 3    10 -0.00503   -0.187    0.173 0.984     0.852      1.12


####################
#shigella AFe>=0.5 only MUAC among bloody diarrhea only
################### var screening, AUC ####
names<-append(x=names, values=c("f4b_muac"))
names <- names[!names %in% c("f4b_haz")]

shigella0.5MUAC.bloody <- CPR.funct(data=cases.bloody,outcome="shigella_afe0.5",iter=100,nvars_opts=c(2,5,10))
shigella0.5MUAC.bloody[["df_imps"]]
shigella0.5MUAC.bloody[["AUC_df"]]
shigella0.5MUAC.bloody[["calib"]]

# names      var_red
# 1               base_age 49.487381728
# 2               f4b_temp 11.815091696
# 3               f4b_resp 11.391826013
# 4               f4b_muac 10.783980958
# 5          f4a_ppl_house  6.698996576
# 6           f4a_drh_days  5.378919468
# 7          f4a_slp_rooms  4.797682384
# 8         f4a_max_stools  4.662456184
# 9          f4a_share_fac  4.462623712
# 10         f4a_prim_schl  4.288641231
# 11      f4a_yng_children  3.720360942
# 12                  site  3.576883191
# 13     f4a_fuel_charcoal  3.446664433
# 14         f4a_fac_waste  2.983785575
# 15         f4a_breastfed  2.924300078
# 16      f4a_hometrt_none  2.862879286
# 17        f4a_offr_drink  2.738041041
# 18         f4b_recommend  2.735695651
# 19       f4a_hometrt_ors  2.538606948
# 20          f4a_dad_live  2.515533341
# 21        f4a_disp_feces  2.316160549
# 22         f4a_drh_cough  2.242855955
# 23        f4a_drh_strain  2.157689874
# 24      f4a_seek_outside  2.048896255
# 25              f4b_eyes  2.043008239
# 26          f4a_ms_water  2.008948699
# 27        f4a_trt_method  1.911228455
# 28         f4a_drh_vomit  1.850511066
# 29       f4a_water_avail  1.822725894
# 30          f4a_cur_skin  1.791526218
# 31             f3_gender  1.791011815
# 32     f4a_drh_bellypain  1.779143786
# 33          f4a_wash_use  1.771409666
# 34      f4a_house_agland  1.709198766
# 35       f4a_house_phone  1.698751741
# 36       f4a_cur_thirsty  1.691442439
# 37         f4a_wash_cook  1.673173129
# 38      f4a_cur_restless  1.671039733
# 39          f4a_ani_fowl  1.655469612
# 40       f4a_house_radio  1.632271809
# 41        f4a_wash_child  1.617130506
# 42      f4a_cur_drymouth  1.602303900
# 43          f4a_wash_eat  1.593218127
# 44        f4a_house_bike  1.567398573
# 45             f4b_mouth  1.550966360
# 46        f4a_wash_nurse  1.542761399
# 47        f4a_house_elec  1.531068362
# 48        f4a_house_tele  1.526918135
# 49        f4a_fuel_grass  1.513954635
# 50     f4a_hometrt_othr1  1.482459435
# 51           f4a_ani_cat  1.479264239
# 52       f4a_house_scoot  1.476024429
# 53           f4a_ani_dog  1.465615844
# 54          f4a_wash_def  1.438022396
# 55      f4a_drh_restless  1.428393874
# 56             f4a_floor  1.422328575
# 57         f4a_wash_othr  1.417793138
# 58          f4a_ani_goat  1.410786117
# 59  f4a_drh_lethrgy_miss  1.409493402
# 60        f4a_drh_thirst  1.409057540
# 61       f4a_ani_rodents  1.392446161
# 62            f4b_mental  1.352833432
# 63        f4a_seek_pharm  1.306012750
# 64         f4a_fuel_wood  1.295941900
# 65       f4a_store_water  1.284614785
# 66           f4a_ani_cow  1.239985818
# 67             f4b_admit  1.233322936
# 68       f4a_wash_animal  1.193356791
# 69         f4a_fuel_dung  1.158174893
# 70      f4a_house_fridge  1.157387158
# 71     f4a_drh_lessdrink  1.152247175
# 72         f4a_fuel_crop  1.146867955
# 73    f4a_water_deepwell  1.137543227
# 74      f4a_water_pubtap  1.084320889
# 75        f4a_hometrt_ab  1.067489630
# 76         f4a_trt_water  1.055577979
# 77            f4a_ani_no  0.981788677
# 78           f3_drh_hosp  0.957647077
# 79   f4a_water_shallwell  0.933240171
# 80         f4a_ani_other  0.898335241
# 81         f4a_ani_sheep  0.865474918
# 82      f4a_seek_privdoc  0.861617332
# 83         f3_drh_turgor  0.826566752
# 84              f4b_skin  0.817841949
# 85    f4a_cur_fastbreath  0.813933451
# 86         f4a_house_car  0.805773487
# 87       f4a_fuel_natgas  0.801004336
# 88     f4a_hometrt_maize  0.773925165
# 89      f4a_relationship  0.747106156
# 90        f4b_under_nutr  0.695133902
# 91        f4a_water_yard  0.688653566
# 92       f4a_water_house  0.657451773
# 93      f4a_water_bought  0.647254616
# 94      f4a_hometrt_herb  0.635806610
# 95        f4a_house_cart  0.635329076
# 96      f4a_hometrt_zinc  0.607433718
# 97        f4a_water_rain  0.585357023
# 98             f3_drh_iv  0.545870369
# 99         f4a_fuel_kero  0.489953382
# 100       f4a_seek_remdy  0.475543937
# 101     f4a_fuel_propane  0.474119268
# 102    f4a_water_pubwell  0.472893918
# 103      f4a_seek_healer  0.452571230
# 104       f4a_house_none  0.436124655
# 105     f4a_drh_prolapse  0.406407161
# 106         f4a_seek_doc  0.401294955
# 107       f4a_water_bore  0.382009714
# 108        f4a_drh_consc  0.351966228
# 109    f4a_hometrt_othr2  0.345179633
# 110  f4a_water_prospring  0.343125020
# 111           f4b_rectal  0.307479239
# 112  f4a_hometrt_othrliq  0.298868968
# 113       f4a_seek_other  0.294479548
# 114       f4a_fuel_other  0.290021599
# 115        f4a_fuel_coal  0.289583402
# 116      f4a_water_river  0.265629406
# 117         f4a_drh_conv  0.256378798
# 118      f4b_chest_indrw  0.246384630
# 119      f4a_seek_friend  0.212251315
# 120       f4a_water_well  0.209262472
# 121       f4a_water_othr  0.193368112
# 122       f4a_house_boat  0.173322590
# 123    f4a_water_covwell  0.171029642
# 124        f4a_fuel_elec  0.167140230
# 125   f4a_water_covpwell  0.160906207
# 126         f4b_abn_hair  0.155143227
# 127   f4a_water_unspring  0.112348322
# 128       f4a_water_pond  0.111901612
# 129      f4a_fuel_biogas  0.107754902
# 130       f4b_skin_flaky  0.081860772
# 131     f4a_hometrt_milk  0.036568350
# 132          f4b_bipedal  0.007521521
# 133        f4a_drh_blood  0.000000000

# AUC          SE     lower     upper level Model nvar
# 1 0.7284832 0.003535680 0.7215533 0.7354130  0.95    LR    2
# 2 0.7329859 0.003411624 0.7262992 0.7396726  0.95    LR    5
# 3 0.7423445 0.003249645 0.7359753 0.7487136  0.95    LR   10
# 4 0.7489354 0.003335767 0.7423974 0.7554733  0.95    RF    2
# 5 0.7792569 0.003099256 0.7731825 0.7853313  0.95    RF    5
# 6 0.8101189 0.002844068 0.8045447 0.8156932  0.95    RF   10

# nvar    intc intc_LCI intc_UCI slope slope_LCI slope_UCI
# <dbl>   <dbl>    <dbl>    <dbl> <dbl>     <dbl>     <dbl>
#   1     2 0.00665   -0.266    0.284 1.02      0.647      1.44
# 2     5 0.00245   -0.273    0.283 1.02      0.668      1.41
# 3    10 0.00141   -0.282    0.289 0.908     0.608      1.24

####################
#shigella AFe>=0.5 only MUAC among non-bloody diarrhea only
################### var screening, AUC ####
names<-append(x=names, values=c("f4b_muac"))
names <- names[!names %in% c("f4b_haz")]

shigella0.5MUAC.nonbloody <- CPR.funct(data=cases.nonbloody,outcome="shigella_afe0.5",iter=100,nvars_opts=c(2,5,10))
shigella0.5MUAC.nonbloody[["df_imps"]]
shigella0.5MUAC.nonbloody[["AUC_df"]]
shigella0.5MUAC.nonbloody[["calib"]]

# names     var_red
# 1               base_age 30.17277442
# 2               f4b_muac 20.83554668
# 3               f4b_resp 18.73251870
# 4               f4b_temp 17.85804082
# 5          f4a_ppl_house 15.15095029
# 6           f4a_drh_days  9.74143766
# 7          f4a_slp_rooms  8.81493571
# 8          f4a_share_fac  8.66178235
# 9       f4a_yng_children  8.03894598
# 10         f4a_breastfed  7.39967688
# 11         f4a_prim_schl  6.99034570
# 12                  site  5.51108513
# 13        f4a_offr_drink  5.42350525
# 14         f4a_drh_vomit  5.13851133
# 15       f4a_water_avail  4.85952288
# 16        f4a_max_stools  4.58562514
# 17        f4a_trt_method  4.53543894
# 18          f4a_ms_water  4.33117567
# 19        f4a_disp_feces  4.10676091
# 20         f4b_recommend  3.88856211
# 21          f4a_dad_live  3.86732956
# 22         f4a_fac_waste  3.81072661
# 23             f4b_mouth  3.64184244
# 24         f4a_drh_cough  3.58751567
# 25      f4a_relationship  3.56607235
# 26            f4b_mental  3.50456879
# 27              f4b_skin  3.32599229
# 28     f4a_drh_bellypain  3.22065183
# 29             f3_gender  3.09814066
# 30          f4a_wash_use  3.06054762
# 31      f4a_cur_drymouth  3.00001540
# 32          f4a_drh_conv  2.95661366
# 33         f4a_wash_cook  2.94857980
# 34        f4a_house_tele  2.94738343
# 35        f4a_drh_thirst  2.93869351
# 36        f4a_wash_child  2.92727096
# 37        f4a_wash_nurse  2.91323775
# 38  f4a_drh_lethrgy_miss  2.89358872
# 39          f4a_wash_def  2.88594822
# 40       f4a_cur_thirsty  2.86102249
# 41        f4a_house_bike  2.83557731
# 42       f4a_hometrt_ors  2.82496920
# 43      f4a_hometrt_none  2.69184557
# 44          f4a_ani_fowl  2.68332089
# 45              f4b_eyes  2.64026105
# 46      f4a_drh_restless  2.63079283
# 47         f4a_fuel_wood  2.61624660
# 48      f4a_cur_restless  2.53247565
# 49        f4b_under_nutr  2.51894939
# 50    f4a_cur_fastbreath  2.48636203
# 51      f4a_seek_outside  2.48394180
# 52           f4a_ani_cat  2.46185636
# 53          f4a_cur_skin  2.45636841
# 54          f4a_wash_eat  2.43699368
# 55       f4a_store_water  2.43327393
# 56      f4a_water_pubtap  2.41464311
# 57       f4a_house_scoot  2.40378160
# 58       f4a_house_phone  2.39181826
# 59       f4a_ani_rodents  2.35221285
# 60         f4a_trt_water  2.31331493
# 61           f3_drh_hosp  2.31203550
# 62       f4a_house_radio  2.29746932
# 63          f4a_ani_goat  2.26740703
# 64        f4a_house_elec  2.22768464
# 65     f4a_hometrt_maize  2.21038920
# 66           f4a_ani_dog  2.19048194
# 67      f4a_house_fridge  2.12706452
# 68     f4a_drh_lessdrink  2.03076441
# 69           f4a_ani_cow  2.02186010
# 70             f3_drh_iv  1.99292967
# 71         f3_drh_turgor  1.97386212
# 72     f4a_fuel_charcoal  1.96486891
# 73             f4a_floor  1.94452498
# 74      f4a_house_agland  1.80644284
# 75        f4a_drh_strain  1.76100631
# 76        f4a_seek_pharm  1.73625338
# 77             f4b_admit  1.72925437
# 78         f4a_ani_sheep  1.69494022
# 79     f4a_hometrt_othr1  1.69466263
# 80        f4a_water_yard  1.68192790
# 81       f4a_fuel_natgas  1.64770559
# 82        f4a_hometrt_ab  1.58101635
# 83         f4a_house_car  1.57200417
# 84       f4a_seek_healer  1.56221776
# 85       f4a_water_house  1.55127485
# 86    f4a_water_deepwell  1.53476423
# 87      f4a_hometrt_herb  1.52832462
# 88            f4a_ani_no  1.46392084
# 89        f4a_fuel_grass  1.41824521
# 90         f4a_ani_other  1.40509401
# 91      f4a_water_bought  1.40401741
# 92        f4a_house_cart  1.36553559
# 93         f4a_wash_othr  1.35362323
# 94      f4a_hometrt_milk  1.29194613
# 95         f4a_fuel_crop  1.25754077
# 96     f4a_water_pubwell  1.24403570
# 97         f4a_fuel_kero  1.22439389
# 98       f4a_wash_animal  1.13765630
# 99        f4a_water_bore  1.12141232
# 100        f4a_fuel_elec  1.10654821
# 101       f4a_water_well  1.08917805
# 102         f4b_abn_hair  1.07494154
# 103       f4a_house_none  1.01071857
# 104     f4a_fuel_propane  0.99505337
# 105       f4b_skin_flaky  0.97740737
# 106        f4a_drh_consc  0.91249512
# 107       f4a_seek_remdy  0.90615202
# 108       f4a_fuel_other  0.90497547
# 109    f4a_water_covwell  0.86223323
# 110      f4b_chest_indrw  0.84349524
# 111        f4a_fuel_coal  0.79834885
# 112     f4a_hometrt_zinc  0.75257074
# 113   f4a_water_covpwell  0.72875479
# 114     f4a_seek_privdoc  0.72802097
# 115         f4a_seek_doc  0.71258341
# 116       f4a_water_rain  0.70588323
# 117    f4a_hometrt_othr2  0.69975474
# 118  f4a_water_shallwell  0.66899303
# 119      f4a_fuel_biogas  0.64428563
# 120        f4a_fuel_dung  0.58842562
# 121       f4a_water_pond  0.56858007
# 122       f4a_seek_other  0.53218073
# 123       f4a_water_othr  0.50651894
# 124      f4a_water_river  0.49027640
# 125       f4a_house_boat  0.45728936
# 126          f4b_bipedal  0.42565874
# 127  f4a_hometrt_othrliq  0.38209705
# 128     f4a_drh_prolapse  0.35834345
# 129      f4a_seek_friend  0.24301874
# 130           f4b_rectal  0.13860520
# 131   f4a_water_unspring  0.12243901
# 132  f4a_water_prospring  0.02085684
# 133        f4a_drh_blood  0.00000000

# AUC          SE     lower     upper level Model nvar
# 1 0.6582715 0.002368217 0.6536299 0.6629131  0.95    LR    2
# 2 0.6413242 0.002526771 0.6363718 0.6462765  0.95    LR    5
# 3 0.6355184 0.002645600 0.6303332 0.6407037  0.95    LR   10
# 4 0.6450741 0.002569766 0.6400374 0.6501107  0.95    RF    2
# 5 0.6716931 0.002498501 0.6667962 0.6765901  0.95    RF    5
# 6 0.6959667 0.002461913 0.6911415 0.7007920  0.95    RF   10

# nvar    intc intc_LCI intc_UCI slope slope_LCI slope_UCI
# <dbl>   <dbl>    <dbl>    <dbl> <dbl>     <dbl>     <dbl>
#   1     2 0.00611   -0.202    0.205 1.01      0.510      1.49
# 2     5 0.00607   -0.202    0.205 0.965     0.485      1.44
# 3    10 0.00416   -0.205    0.204 0.839     0.444      1.24

####################
#shigella AFe>=0.5 w/o Bangladesh
################### var screening, AUC ####
names<-append(x=names, values=c("f4b_muac"))
names <- names[!names %in% c("f4b_haz")]

shigella0.5MUAC.noBdesh <- CPR.funct(data=cases.noBdesh,outcome="shigella_afe0.5",iter=100,nvars_opts=c(2,5,10))
shigella0.5MUAC.noBdesh[["df_imps"]]
shigella0.5MUAC.noBdesh[["AUC_df"]]
shigella0.5MUAC.noBdesh[["calib"]]

# names    var_red
# 1          f4a_drh_blood 52.2065849
# 2               base_age 43.6650581
# 3               f4b_muac 25.9598038
# 4               f4b_resp 25.3412828
# 5               f4b_temp 22.3077370
# 6          f4a_ppl_house 19.8126421
# 7           f4a_drh_days 12.4915731
# 8          f4a_slp_rooms 12.0306318
# 9          f4a_share_fac 11.2593332
# 10      f4a_yng_children 11.1388318
# 11         f4a_breastfed  9.6270325
# 12         f4a_prim_schl  9.2480730
# 13              f4b_eyes  8.5194449
# 14         f4a_drh_vomit  8.4465978
# 15         f4b_recommend  8.2262351
# 16        f4a_offr_drink  7.4949236
# 17                  site  7.0384947
# 18       f4a_water_avail  6.4830556
# 19        f4a_max_stools  6.2432749
# 20          f4a_ms_water  5.9491976
# 21        f4a_trt_method  5.8274137
# 22        f4a_disp_feces  5.4929251
# 23             f4b_mouth  5.0540020
# 24          f4a_dad_live  5.0398182
# 25         f4a_drh_cough  4.8860129
# 26            f4b_mental  4.6860704
# 27     f4a_drh_bellypain  4.3459305
# 28         f4a_fac_waste  4.2724675
# 29      f4a_relationship  4.1741407
# 30         f4a_wash_cook  4.1258181
# 31        f4a_wash_nurse  4.1125712
# 32             f3_gender  4.0769465
# 33       f4a_cur_thirsty  3.9971378
# 34        f4a_wash_child  3.9079218
# 35              f4b_skin  3.8956037
# 36      f4a_cur_drymouth  3.8676687
# 37          f4a_wash_use  3.8428914
# 38      f4a_hometrt_none  3.8212841
# 39  f4a_drh_lethrgy_miss  3.8181242
# 40       f4a_hometrt_ors  3.7268017
# 41        f4a_drh_thirst  3.6699235
# 42          f4a_wash_def  3.6435918
# 43      f4a_drh_restless  3.5848727
# 44        f4a_house_bike  3.5639236
# 45        f4a_house_tele  3.5590188
# 46          f4a_cur_skin  3.5456021
# 47          f4a_ani_fowl  3.4569767
# 48      f4a_water_pubtap  3.4017561
# 49      f4a_cur_restless  3.3563213
# 50         f4a_fuel_wood  3.3528501
# 51           f4a_ani_cat  3.3260002
# 52          f4a_wash_eat  3.3225554
# 53       f4a_house_phone  3.3157733
# 54         f4a_trt_water  3.2597931
# 55     f4a_fuel_charcoal  3.2517361
# 56     f4a_hometrt_maize  3.2402650
# 57       f4a_house_radio  3.2067234
# 58     f4a_drh_lessdrink  3.1603415
# 59      f4a_seek_outside  3.1025779
# 60       f4a_ani_rodents  3.0646878
# 61      f4a_house_fridge  3.0278800
# 62        f4b_under_nutr  3.0101159
# 63        f4a_drh_strain  2.9660151
# 64             f4b_admit  2.9658089
# 65       f4a_store_water  2.9631840
# 66       f4a_house_scoot  2.9428409
# 67           f3_drh_hosp  2.9091631
# 68          f4a_ani_goat  2.8918828
# 69    f4a_cur_fastbreath  2.8635917
# 70           f4a_ani_dog  2.8611451
# 71             f4a_floor  2.7240862
# 72         f3_drh_turgor  2.6795889
# 73        f4a_house_elec  2.6084529
# 74       f4a_fuel_natgas  2.5698005
# 75             f3_drh_iv  2.4653186
# 76           f4a_ani_cow  2.4157896
# 77     f4a_hometrt_othr1  2.3770606
# 78      f4a_house_agland  2.3684504
# 79         f4a_ani_sheep  2.3078578
# 80        f4a_water_yard  2.2554141
# 81            f4a_ani_no  2.2119755
# 82        f4a_house_cart  2.0973108
# 83         f4a_house_car  2.0922153
# 84      f4a_water_bought  2.0916307
# 85         f4a_ani_other  2.0789678
# 86      f4a_hometrt_herb  2.0707278
# 87       f4a_water_house  1.9514288
# 88        f4a_hometrt_ab  1.8439642
# 89       f4a_seek_healer  1.8125246
# 90        f4a_seek_pharm  1.8011877
# 91         f4a_wash_othr  1.7978902
# 92         f4a_fuel_kero  1.7482165
# 93     f4a_water_pubwell  1.6492410
# 94       f4a_wash_animal  1.5592966
# 95        f4a_water_bore  1.5111221
# 96    f4a_water_deepwell  1.4553908
# 97      f4a_fuel_propane  1.4345655
# 98          f4b_abn_hair  1.4101345
# 99      f4a_seek_privdoc  1.3533073
# 100       f4a_water_well  1.3366042
# 101     f4a_hometrt_milk  1.2533459
# 102        f4a_fuel_coal  1.1420934
# 103        f4a_fuel_elec  1.1010953
# 104       f4b_skin_flaky  1.0952646
# 105       f4a_water_rain  1.0700378
# 106        f4a_drh_consc  1.0425869
# 107      f4b_chest_indrw  1.0389090
# 108       f4a_fuel_other  1.0133479
# 109         f4a_drh_conv  0.9977284
# 110    f4a_water_covwell  0.9684282
# 111    f4a_hometrt_othr2  0.9555697
# 112     f4a_drh_prolapse  0.9540295
# 113         f4a_seek_doc  0.9120035
# 114   f4a_water_covpwell  0.8297272
# 115       f4a_house_none  0.7976912
# 116      f4a_fuel_biogas  0.7518731
# 117       f4a_water_othr  0.7256207
# 118       f4a_seek_other  0.7028342
# 119       f4a_seek_remdy  0.6957630
# 120      f4a_water_river  0.6840705
# 121     f4a_hometrt_zinc  0.6416151
# 122  f4a_hometrt_othrliq  0.6287140
# 123       f4a_water_pond  0.5506557
# 124           f4b_rectal  0.5083807
# 125          f4b_bipedal  0.4884048
# 126       f4a_fuel_grass  0.4844065
# 127       f4a_house_boat  0.4489716
# 128        f4a_fuel_crop  0.3509313
# 129   f4a_water_unspring  0.3225936
# 130      f4a_seek_friend  0.2764423
# 131  f4a_water_shallwell  0.2158289
# 132        f4a_fuel_dung  0.2045654
# 133  f4a_water_prospring  0.1197469

# AUC          SE     lower     upper level Model nvar
# 1 0.7484874 0.002087551 0.7443958 0.7525789  0.95    LR    2
# 2 0.7356321 0.002173127 0.7313728 0.7398913  0.95    LR    5
# 3 0.7349915 0.002198713 0.7306821 0.7393009  0.95    LR   10
# 4 0.7661968 0.002074901 0.7621301 0.7702635  0.95    RF    2
# 5 0.7584421 0.002082867 0.7543598 0.7625245  0.95    RF    5
# 6 0.7707181 0.002021571 0.7667559 0.7746803  0.95    RF   10

# nvar     intc intc_LCI intc_UCI slope slope_LCI slope_UCI
# <dbl>    <dbl>    <dbl>    <dbl> <dbl>     <dbl>     <dbl>
#   1     2 -0.00701   -0.194    0.175 0.991     0.804      1.18
# 2     5 -0.00805   -0.196    0.174 0.985     0.799      1.18
# 3    10 -0.00843   -0.196    0.174 0.973     0.790      1.16

####################
#shigella AFe>=0.5 only MUAC + clinically observed stool variables
################### var screening, AUC ####
names<-append(x=names, values=c("f4b_muac","f11_consistency","f11_blood","f11_pus","f11_mucus"))
names <- names[!names %in% c("f4b_haz")]

shigella0.5MUAC.stool <- CPR.funct(data=cases,outcome="shigella_afe0.5",iter=100,nvars_opts=c(2,5,10))
shigella0.5MUAC.stool[["df_imps"]]
shigella0.5MUAC.stool[["AUC_df"]]
shigella0.5MUAC.stool[["calib"]]

# names    var_red
# 1               base_age 84.8940886
# 2          f4a_drh_blood 68.0304613
# 3              f11_blood 33.5455531
# 4               f4b_muac 31.4190176
# 5               f4b_resp 29.4989160
# 6               f4b_temp 27.6304745
# 7               f4b_eyes 22.6578714
# 8          f4a_ppl_house 21.6348432
# 9                   site 20.1199321
# 10       f11_consistency 19.0553981
# 11          f4a_drh_days 15.2254599
# 12         f4a_slp_rooms 13.6391656
# 13         f4a_share_fac 12.9247254
# 14      f4a_yng_children 11.9206418
# 15         f4a_prim_schl 11.3639218
# 16         f4a_breastfed 10.9663492
# 17        f4a_drh_strain 10.9341634
# 18             f11_mucus 10.2588061
# 19             f4b_mouth  9.7640158
# 20         f4a_drh_vomit  9.7386210
# 21        f4a_max_stools  8.9402345
# 22        f4a_fuel_grass  8.7607881
# 23        f4a_offr_drink  8.4996279
# 24         f4b_recommend  7.8021227
# 25         f4a_fac_waste  7.7519986
# 26       f4a_store_water  7.7456300
# 27         f4a_fuel_crop  7.2998706
# 28        f4a_disp_feces  6.8408566
# 29     f4a_drh_bellypain  6.8339377
# 30         f4a_drh_cough  6.7842757
# 31       f4a_water_avail  6.6077142
# 32       f4a_cur_thirsty  6.3540809
# 33        f4a_trt_method  6.2588152
# 34          f4a_dad_live  6.2155704
# 35          f4a_ms_water  6.1320915
# 36        f4a_drh_thirst  5.8812858
# 37      f4a_hometrt_none  5.2880669
# 38      f4a_cur_drymouth  5.2078846
# 39            f4b_mental  5.0316368
# 40       f4a_hometrt_ors  4.9100420
# 41             f3_gender  4.8884740
# 42          f4a_wash_use  4.7673742
# 43        f4a_wash_child  4.7272823
# 44         f4a_wash_cook  4.6859558
# 45     f4a_fuel_charcoal  4.6122177
# 46        f4a_wash_nurse  4.5946775
# 47        f4a_house_tele  4.4600747
# 48          f4a_ani_fowl  4.4400416
# 49         f4a_fuel_dung  4.4317356
# 50          f4a_wash_def  4.4132702
# 51        f4a_house_bike  4.3814029
# 52  f4a_drh_lethrgy_miss  4.3573125
# 53      f4a_seek_outside  4.3377964
# 54      f4a_cur_restless  4.2778236
# 55      f4a_drh_restless  4.2351720
# 56              f4b_skin  4.1787027
# 57      f4a_relationship  4.1336765
# 58       f4a_house_phone  4.1193343
# 59       f4a_house_radio  4.1053081
# 60           f4a_ani_cat  4.0801870
# 61          f4a_cur_skin  4.0793158
# 62         f4a_fuel_wood  4.0612032
# 63             f4a_floor  3.9877476
# 64          f4a_wash_eat  3.9681063
# 65        f4a_house_elec  3.7773451
# 66             f4b_admit  3.7700322
# 67           f4a_ani_dog  3.7153623
# 68        f4a_seek_pharm  3.7065811
# 69       f4a_house_scoot  3.6784115
# 70       f4a_ani_rodents  3.6693579
# 71      f4a_house_agland  3.6690116
# 72    f4a_water_deepwell  3.6167810
# 73      f4a_house_fridge  3.6087098
# 74           f3_drh_hosp  3.5731243
# 75          f4a_ani_goat  3.5410815
# 76      f4a_water_pubtap  3.5405641
# 77         f4a_trt_water  3.3949147
# 78           f4a_ani_cow  3.3157871
# 79        f4b_under_nutr  3.2721904
# 80     f4a_hometrt_othr1  3.1597172
# 81     f4a_drh_lessdrink  3.1443312
# 82    f4a_cur_fastbreath  3.1407418
# 83     f4a_hometrt_maize  3.1216194
# 84               f11_pus  3.0964598
# 85             f3_drh_iv  3.0670281
# 86         f3_drh_turgor  2.9916026
# 87            f4a_ani_no  2.8585503
# 88       f4a_fuel_natgas  2.7954899
# 89         f4a_wash_othr  2.7428239
# 90        f4a_hometrt_ab  2.7390971
# 91   f4a_water_shallwell  2.5815699
# 92         f4a_ani_sheep  2.5794589
# 93       f4a_wash_animal  2.3988908
# 94      f4a_hometrt_herb  2.3819973
# 95        f4a_water_yard  2.3414487
# 96       f4a_water_house  2.2724614
# 97         f4a_house_car  2.2282069
# 98          f4a_drh_conv  2.1670609
# 99         f4a_ani_other  2.1468895
# 100       f4a_house_cart  2.1362091
# 101      f4a_seek_healer  1.9834395
# 102     f4a_water_bought  1.9744014
# 103    f4a_water_pubwell  1.7054432
# 104        f4a_fuel_kero  1.6398149
# 105       f4a_house_none  1.5581442
# 106     f4a_hometrt_zinc  1.5387312
# 107       f4a_seek_remdy  1.5230154
# 108     f4a_fuel_propane  1.4797276
# 109       f4a_water_bore  1.4611356
# 110     f4a_seek_privdoc  1.4345997
# 111         f4b_abn_hair  1.3285919
# 112         f4a_seek_doc  1.2805315
# 113       f4a_water_well  1.2421427
# 114       f4a_fuel_other  1.2257344
# 115        f4a_drh_consc  1.1812549
# 116     f4a_hometrt_milk  1.1752984
# 117        f4a_fuel_coal  1.1516957
# 118       f4a_water_rain  1.0818431
# 119        f4a_fuel_elec  1.0673963
# 120      f4b_chest_indrw  1.0381236
# 121    f4a_hometrt_othr2  1.0071367
# 122       f4b_skin_flaky  0.9685606
# 123    f4a_water_covwell  0.8878524
# 124     f4a_drh_prolapse  0.8674993
# 125   f4a_water_covpwell  0.7948301
# 126      f4a_water_river  0.7586936
# 127      f4a_fuel_biogas  0.7475219
# 128       f4a_water_othr  0.7067706
# 129       f4a_seek_other  0.6903315
# 130  f4a_hometrt_othrliq  0.6455968
# 131       f4a_water_pond  0.6220991
# 132       f4a_house_boat  0.6117178
# 133          f4b_bipedal  0.4620638
# 134           f4b_rectal  0.3936521
# 135      f4a_seek_friend  0.3034177
# 136   f4a_water_unspring  0.2839810
# 137  f4a_water_prospring  0.1614644

# AUC          SE     lower     upper level Model nvar
# 1 0.8053309 0.001603341 0.8021884 0.8084734  0.95    LR    2
# 2 0.8026321 0.001586853 0.7995220 0.8057423  0.95    LR    5
# 3 0.8050455 0.001636291 0.8018384 0.8082525  0.95    LR   10
# 4 0.8091665 0.001608834 0.8060133 0.8123198  0.95    RF    2
# 5 0.8144310 0.001566885 0.8113600 0.8175021  0.95    RF    5
# 6 0.8346665 0.001502414 0.8317218 0.8376112  0.95    RF   10

# nvar     intc intc_LCI intc_UCI slope slope_LCI slope_UCI
# <dbl>    <dbl>    <dbl>    <dbl> <dbl>     <dbl>     <dbl>
#   1     2 -0.00195   -0.165    0.159  1.02     0.881      1.16
# 2     5 -0.00399   -0.168    0.157  1.01     0.876      1.16
# 3    10 -0.00403   -0.171    0.160  1.00     0.866      1.14

####################
#shigella AFe>=0.5 only MUAC + season
################### var screening, AUC ####
names<-append(x=names, values=c("f4b_muac","season"))
names <- names[!names %in% c("f4b_haz")]

shigella0.5MUAC.season <- CPR.funct(data=cases,outcome="shigella_afe0.5",iter=100,nvars_opts=c(2,5,10))
shigella0.5MUAC.season[["df_imps"]]
shigella0.5MUAC.season[["AUC_df"]]
shigella0.5MUAC.season[["calib"]]

# names    var_red
# 1               base_age 89.4565712
# 2          f4a_drh_blood 79.7900902
# 3               f4b_muac 32.8619983
# 4               f4b_resp 31.0262245
# 5               f4b_temp 29.0704308
# 6               f4b_eyes 26.9344592
# 7          f4a_ppl_house 22.9093575
# 8                   site 20.9258502
# 9           f4a_drh_days 15.8976535
# 10         f4a_slp_rooms 14.1816129
# 11        f4a_drh_strain 13.6344893
# 12         f4a_share_fac 13.3438681
# 13      f4a_yng_children 12.4285789
# 14         f4a_prim_schl 11.6174517
# 15        f4a_fuel_grass 11.0932880
# 16         f4a_breastfed 10.9644904
# 17         f4a_drh_vomit 10.8875194
# 18             f4b_mouth 10.1691545
# 19        f4a_max_stools  9.2604571
# 20        f4a_offr_drink  8.9764289
# 21         f4a_fuel_crop  8.7695395
# 22         f4b_recommend  8.6434577
# 23       f4a_store_water  8.0080214
# 24         f4a_fac_waste  7.9534170
# 25     f4a_drh_bellypain  7.9001414
# 26         f4a_drh_cough  7.1370604
# 27       f4a_water_avail  7.1025667
# 28        f4a_disp_feces  6.8104177
# 29          f4a_dad_live  6.6706645
# 30                season  6.6378030
# 31        f4a_trt_method  6.6179745
# 32        f4a_drh_thirst  6.5693402
# 33       f4a_cur_thirsty  6.5336820
# 34          f4a_ms_water  6.4407988
# 35            f4b_mental  5.6616772
# 36      f4a_cur_drymouth  5.4069090
# 37      f4a_hometrt_none  5.3131044
# 38       f4a_hometrt_ors  5.1249094
# 39             f3_gender  5.0539985
# 40         f4a_wash_cook  4.8719461
# 41        f4a_wash_nurse  4.8340216
# 42  f4a_drh_lethrgy_miss  4.8244624
# 43        f4a_wash_child  4.7601142
# 44          f4a_wash_use  4.7551690
# 45          f4a_ani_fowl  4.6173799
# 46        f4a_house_tele  4.5997375
# 47        f4a_house_bike  4.5970785
# 48      f4a_seek_outside  4.5910291
# 49          f4a_wash_def  4.5740378
# 50      f4a_drh_restless  4.5443830
# 51              f4b_skin  4.5059879
# 52         f4a_fuel_dung  4.4772278
# 53         f4a_fuel_wood  4.3275091
# 54     f4a_fuel_charcoal  4.3130009
# 55           f4a_ani_cat  4.2956927
# 56          f4a_cur_skin  4.2932979
# 57          f4a_wash_eat  4.2712383
# 58      f4a_cur_restless  4.2649484
# 59      f4a_relationship  4.2464349
# 60       f4a_house_phone  4.2379168
# 61             f4b_admit  4.1503681
# 62       f4a_house_radio  4.1413561
# 63           f3_drh_hosp  4.0296011
# 64             f4a_floor  4.0293165
# 65    f4a_water_deepwell  3.9965145
# 66       f4a_ani_rodents  3.9328588
# 67        f4a_house_elec  3.8966068
# 68      f4a_water_pubtap  3.8331577
# 69           f4a_ani_dog  3.8162612
# 70       f4a_house_scoot  3.7388303
# 71        f4a_seek_pharm  3.7276415
# 72          f4a_ani_goat  3.7144952
# 73      f4a_house_fridge  3.6507694
# 74      f4a_house_agland  3.5673063
# 75         f4a_trt_water  3.5632302
# 76           f4a_ani_cow  3.5335184
# 77         f3_drh_turgor  3.4515547
# 78             f3_drh_iv  3.4116157
# 79     f4a_hometrt_maize  3.3911572
# 80     f4a_hometrt_othr1  3.3475625
# 81     f4a_drh_lessdrink  3.3005182
# 82    f4a_cur_fastbreath  3.2782899
# 83        f4b_under_nutr  3.2355792
# 84        f4a_hometrt_ab  3.0043446
# 85   f4a_water_shallwell  2.9520358
# 86       f4a_fuel_natgas  2.8968924
# 87            f4a_ani_no  2.8942771
# 88         f4a_wash_othr  2.7144409
# 89         f4a_ani_sheep  2.6801475
# 90       f4a_wash_animal  2.4405806
# 91          f4a_drh_conv  2.4390954
# 92        f4a_water_yard  2.4314343
# 93      f4a_hometrt_herb  2.3749419
# 94         f4a_house_car  2.3210905
# 95        f4a_house_cart  2.2891471
# 96       f4a_water_house  2.2663657
# 97         f4a_ani_other  2.1713774
# 98      f4a_water_bought  2.1601224
# 99       f4a_seek_healer  1.9852496
# 100        f4a_fuel_kero  1.9355259
# 101    f4a_water_pubwell  1.8273409
# 102     f4a_hometrt_zinc  1.6248922
# 103     f4a_seek_privdoc  1.6129964
# 104       f4a_water_bore  1.5626634
# 105     f4a_fuel_propane  1.5467135
# 106       f4a_house_none  1.5298889
# 107         f4b_abn_hair  1.4361094
# 108       f4a_seek_remdy  1.3813345
# 109       f4a_water_well  1.3524969
# 110         f4a_seek_doc  1.3258216
# 111     f4a_hometrt_milk  1.2665022
# 112        f4a_fuel_coal  1.2506745
# 113       f4a_fuel_other  1.2478298
# 114        f4a_drh_consc  1.2122703
# 115    f4a_hometrt_othr2  1.2039403
# 116       f4a_water_rain  1.1649155
# 117        f4a_fuel_elec  1.1318355
# 118       f4b_skin_flaky  1.1211832
# 119      f4b_chest_indrw  1.0671366
# 120     f4a_drh_prolapse  1.0301407
# 121    f4a_water_covwell  0.9777440
# 122   f4a_water_covpwell  0.8318790
# 123       f4a_seek_other  0.7543726
# 124      f4a_water_river  0.7454699
# 125       f4a_water_othr  0.7171799
# 126      f4a_fuel_biogas  0.6966764
# 127  f4a_hometrt_othrliq  0.6941170
# 128       f4a_house_boat  0.6742493
# 129       f4a_water_pond  0.5987644
# 130           f4b_rectal  0.4840726
# 131          f4b_bipedal  0.4839404
# 132      f4a_seek_friend  0.3039388
# 133   f4a_water_unspring  0.3011117
# 134  f4a_water_prospring  0.1610185

# AUC          SE     lower     upper level Model nvar
# 1 0.8026307 0.001610827 0.7994735 0.8057879  0.95    LR    2
# 2 0.7928096 0.001672486 0.7895316 0.7960877  0.95    LR    5
# 3 0.7985253 0.001669045 0.7952540 0.8017965  0.95    LR   10
# 4 0.8062690 0.001617086 0.8030996 0.8094385  0.95    RF    2
# 5 0.8093772 0.001585937 0.8062688 0.8124855  0.95    RF    5
# 6 0.8341626 0.001494842 0.8312327 0.8370924  0.95    RF   10

# nvar    intc intc_LCI intc_UCI slope slope_LCI slope_UCI
# <dbl>   <dbl>    <dbl>    <dbl> <dbl>     <dbl>     <dbl>
#   1     2 0.00528   -0.158    0.166 1.00      0.863      1.14
# 2     5 0.00450   -0.160    0.166 0.994     0.858      1.14
# 3    10 0.00598   -0.160    0.169 0.986     0.852      1.13

####################
#shigella AFe>=0.5 only MUAC Gambia only
################### var screening, AUC ####
names<-append(x=names, values=c("f4b_muac")) #add MUAC to list
names <- names[!names %in% c("f4b_haz")] #take HAZ off list

summary(Gambia$f4b_recommend)
Gambia_complete <- Gambia %>% filter(f4b_recommend!=2 & f4b_recommend!=3)
summary(Gambia_complete$f4b_recommend)

shigella0.5MUAC.Gambia <- CPR.funct(data=Gambia_complete,outcome="shigella_afe0.5",iter=100,nvars_opts=c(2,5,10))
shigella0.5MUAC.Gambia[["df_imps"]]
shigella0.5MUAC.Gambia[["AUC_df"]]
shigella0.5MUAC.Gambia[["calib"]]

# names      var_red
# 1          f4a_drh_blood 1.131027e+01
# 2               base_age 7.561273e+00
# 3               f4b_muac 4.423404e+00
# 4          f4a_ppl_house 4.275836e+00
# 5               f4b_resp 4.021998e+00
# 6               f4b_temp 3.688919e+00
# 7          f4a_slp_rooms 2.990772e+00
# 8       f4a_yng_children 2.483914e+00
# 9               f4b_eyes 2.433228e+00
# 10         f4b_recommend 2.388908e+00
# 11          f4a_drh_days 2.237372e+00
# 12        f4a_offr_drink 1.994945e+00
# 13        f4a_drh_strain 1.667715e+00
# 14         f4a_breastfed 1.627909e+00
# 15         f4a_drh_vomit 1.256424e+00
# 16          f4a_dad_live 1.172836e+00
# 17          f4a_ms_water 1.150703e+00
# 18             f4b_mouth 1.112269e+00
# 19         f4a_prim_schl 1.094145e+00
# 20       f4a_water_avail 1.033758e+00
# 21              f4b_skin 9.220463e-01
# 22      f4a_hometrt_none 8.149884e-01
# 23         f4a_drh_cough 8.145235e-01
# 24     f4a_drh_bellypain 8.095507e-01
# 25         f4a_ani_other 7.946037e-01
# 26        f4a_wash_child 7.783255e-01
# 27        f4a_max_stools 7.780915e-01
# 28          f4a_ani_goat 7.706920e-01
# 29  f4a_drh_lethrgy_miss 7.472377e-01
# 30        f4a_wash_nurse 7.191540e-01
# 31      f4a_drh_restless 7.148776e-01
# 32         f4a_house_car 7.145207e-01
# 33          f4a_wash_def 7.092625e-01
# 34     f4a_water_pubwell 6.992615e-01
# 35         f4a_wash_cook 6.930506e-01
# 36             f3_gender 6.885352e-01
# 37      f4a_cur_drymouth 6.836483e-01
# 38       f4a_cur_thirsty 6.829579e-01
# 39           f4a_ani_dog 6.827154e-01
# 40      f4a_water_pubtap 6.745814e-01
# 41         f4a_trt_water 6.642383e-01
# 42        f4a_drh_thirst 6.547712e-01
# 43       f4a_ani_rodents 6.481016e-01
# 44          f4a_wash_use 6.309852e-01
# 45           f4a_ani_cow 6.264302e-01
# 46        f4a_trt_method 6.256436e-01
# 47         f4a_share_fac 6.233195e-01
# 48    f4a_cur_fastbreath 6.222193e-01
# 49           f4a_ani_cat 5.992938e-01
# 50     f4a_hometrt_maize 5.900121e-01
# 51       f4a_house_scoot 5.864444e-01
# 52            f4b_mental 5.823919e-01
# 53       f4a_store_water 5.705545e-01
# 54        f4a_house_tele 5.659718e-01
# 55        f4a_house_elec 5.643208e-01
# 56          f4a_cur_skin 5.544361e-01
# 57         f4a_fac_waste 5.337936e-01
# 58         f4a_ani_sheep 5.277502e-01
# 59      f4a_cur_restless 5.268483e-01
# 60             f3_drh_iv 5.207144e-01
# 61     f4a_drh_lessdrink 5.170350e-01
# 62      f4a_house_agland 5.095124e-01
# 63         f3_drh_turgor 5.065955e-01
# 64      f4a_house_fridge 5.057024e-01
# 65           f3_drh_hosp 5.025558e-01
# 66        f4a_water_yard 4.993523e-01
# 67      f4a_hometrt_herb 4.974813e-01
# 68        f4a_hometrt_ab 4.935845e-01
# 69             f4b_admit 4.927744e-01
# 70        f4b_under_nutr 4.752115e-01
# 71     f4a_hometrt_othr1 4.722977e-01
# 72    f4a_water_deepwell 4.653252e-01
# 73        f4a_house_bike 4.531512e-01
# 74             f4a_floor 4.386000e-01
# 75        f4a_house_cart 4.322310e-01
# 76          f4a_ani_fowl 4.129849e-01
# 77        f4a_seek_pharm 4.093354e-01
# 78        f4a_water_well 3.979474e-01
# 79      f4a_seek_outside 3.800646e-01
# 80      f4a_relationship 3.186331e-01
# 81       f4a_hometrt_ors 2.941194e-01
# 82       f4a_wash_animal 2.803570e-01
# 83        f4b_skin_flaky 2.776176e-01
# 84          f4b_abn_hair 2.575490e-01
# 85          f4a_drh_conv 2.528080e-01
# 86       f4a_house_phone 2.463488e-01
# 87        f4a_water_othr 2.341670e-01
# 88     f4a_fuel_charcoal 2.291105e-01
# 89       f4a_house_radio 2.246577e-01
# 90         f4a_fuel_wood 2.095956e-01
# 91    f4a_water_covpwell 2.072847e-01
# 92         f4a_wash_othr 2.070791e-01
# 93       f4a_seek_healer 2.044898e-01
# 94     f4a_hometrt_othr2 2.024120e-01
# 95        f4a_disp_feces 1.929627e-01
# 96          f4a_wash_eat 1.907069e-01
# 97       f4b_chest_indrw 1.818974e-01
# 98     f4a_water_covwell 1.818126e-01
# 99   f4a_hometrt_othrliq 1.654484e-01
# 100       f4a_seek_remdy 1.416206e-01
# 101      f4a_seek_friend 1.242310e-01
# 102       f4a_water_bore 1.148825e-01
# 103     f4a_drh_prolapse 1.077459e-01
# 104        f4a_fuel_coal 8.431875e-02
# 105           f4a_ani_no 8.380848e-02
# 106           f4b_rectal 6.875105e-02
# 107         f4a_seek_doc 6.189284e-02
# 108     f4a_hometrt_milk 4.787203e-02
# 109        f4a_drh_consc 3.700233e-02
# 110      f4a_water_house 1.933168e-02
# 111       f4a_seek_other 1.621570e-02
# 112     f4a_seek_privdoc 1.261830e-02
# 113          f4b_bipedal 1.235621e-02
# 114  f4a_water_shallwell 3.585400e-03
# 115       f4a_water_pond 3.018954e-03
# 116       f4a_fuel_other 2.038095e-03
# 117       f4a_house_none 1.405830e-03
# 118       f4a_water_rain 7.500000e-04
# 119     f4a_hometrt_zinc 4.450549e-04
# 120                 site 0.000000e+00
# 121       f4a_house_boat 0.000000e+00
# 122        f4a_fuel_elec 0.000000e+00
# 123      f4a_fuel_biogas 0.000000e+00
# 124       f4a_fuel_grass 0.000000e+00
# 125     f4a_fuel_propane 0.000000e+00
# 126        f4a_fuel_dung 0.000000e+00
# 127      f4a_fuel_natgas 0.000000e+00
# 128        f4a_fuel_crop 0.000000e+00
# 129        f4a_fuel_kero 0.000000e+00
# 130  f4a_water_prospring 0.000000e+00
# 131   f4a_water_unspring 0.000000e+00
# 132      f4a_water_river 0.000000e+00
# 133     f4a_water_bought 0.000000e+00

# AUC          SE     lower     upper level Model nvar
# 1 0.7781132 0.004697605 0.7689060 0.7873203  0.95    LR    2
# 2 0.7744907 0.004948161 0.7647925 0.7841889  0.95    LR    5
# 3 0.7779051 0.004958946 0.7681858 0.7876245  0.95    LR   10
# 4 0.8185335 0.004330090 0.8100467 0.8270203  0.95    RF    2
# 5 0.8076314 0.004246390 0.7993086 0.8159541  0.95    RF    5
# 6 0.8103131 0.004319644 0.8018467 0.8187794  0.95    RF   10

# nvar    intc intc_LCI intc_UCI slope slope_LCI slope_UCI
# <dbl>   <dbl>    <dbl>    <dbl> <dbl>     <dbl>     <dbl>
#   1     2 -0.0229   -0.498    0.423 1.02      0.609      1.47
# 2     5 -0.0248   -0.507    0.428 0.993     0.598      1.43
# 3    10 -0.0227   -0.515    0.440 0.908     0.543      1.32


####################
#shigella AFe>=0.5 only MUAC Mali only
################### var screening, AUC ####
names<-append(x=names, values=c("f4b_muac")) #add MUAC to list
names <- names[!names %in% c("f4b_haz")] #take HAZ off list

summary(Mali$f4a_prim_schl)
Mali_complete <- Mali
Mali_complete$f4a_prim_schl <- as.factor(ifelse(Mali_complete$f4a_prim_schl==5,4,Mali_complete$f4a_prim_schl))
summary(Mali_complete$f4a_prim_schl)

shigella0.5MUAC.Mali <- CPR.funct(data=Mali_complete,outcome="shigella_afe0.5",iter=100,nvars_opts=c(2,5,10))
shigella0.5MUAC.Mali[["df_imps"]]
shigella0.5MUAC.Mali[["AUC_df"]]
shigella0.5MUAC.Mali[["calib"]]

# names      var_red
# 1               base_age 7.9963285826
# 2               f4b_muac 5.6872058575
# 3               f4b_resp 5.4093202487
# 4               f4b_temp 5.2787194071
# 5          f4a_ppl_house 4.4563067868
# 6          f4a_slp_rooms 3.2877691626
# 7       f4a_yng_children 2.7648506286
# 8           f4a_drh_days 2.6712308952
# 9          f4a_share_fac 2.3482474022
# 10         f4a_prim_schl 1.9438191402
# 11         f4a_breastfed 1.7686784104
# 12         f4a_drh_blood 1.6782066927
# 13        f4a_offr_drink 1.4007542517
# 14          f4a_dad_live 1.3419060923
# 15      f4a_relationship 1.1949869658
# 16         f4a_drh_vomit 1.1292056713
# 17      f4a_cur_drymouth 1.1132143827
# 18        f4a_max_stools 1.0403446884
# 19        f4a_house_bike 1.0232523637
# 20              f4b_eyes 1.0086475560
# 21             f4b_mouth 0.9848196157
# 22          f4a_ani_fowl 0.9431115206
# 23             f3_gender 0.9349261508
# 24          f4a_ms_water 0.9207318760
# 25         f4a_drh_cough 0.9054621319
# 26         f4a_wash_cook 0.9033641525
# 27        f4a_wash_nurse 0.8941305872
# 28        f4a_wash_child 0.8901197288
# 29         f4a_fuel_wood 0.8812309618
# 30      f4a_hometrt_none 0.8622581814
# 31     f4a_fuel_charcoal 0.8450907857
# 32          f4a_wash_def 0.8438281774
# 33     f4a_drh_bellypain 0.8435073104
# 34       f4a_hometrt_ors 0.8269013158
# 35              f4b_skin 0.8213129950
# 36       f4a_ani_rodents 0.8118751799
# 37       f4a_seek_healer 0.8001766406
# 38        f4a_disp_feces 0.7905943084
# 39        f4a_house_tele 0.7838406223
# 40          f4a_wash_use 0.7836803326
# 41  f4a_drh_lethrgy_miss 0.7774943768
# 42      f4a_house_fridge 0.7733152594
# 43      f4a_drh_restless 0.7670832799
# 44      f4a_water_pubtap 0.7668972828
# 45        f4a_house_cart 0.7585250929
# 46      f4a_cur_restless 0.7443390457
# 47        f4a_house_elec 0.7413816263
# 48       f4a_store_water 0.7413420613
# 49       f4a_house_scoot 0.7233620145
# 50         f4a_fuel_elec 0.7170952761
# 51          f4a_cur_skin 0.6998146216
# 52      f4a_seek_outside 0.6916119714
# 53        f4b_under_nutr 0.6887233800
# 54      f4a_water_bought 0.6853213408
# 55        f4a_drh_thirst 0.6740129830
# 56      f4a_hometrt_herb 0.6582906603
# 57        f4a_fuel_other 0.6201866615
# 58         f4a_house_car 0.6194354947
# 59         f3_drh_turgor 0.6103506698
# 60            f4b_mental 0.6071671632
# 61          f4a_ani_goat 0.5934164474
# 62        f4a_drh_strain 0.5928160704
# 63         f4a_ani_sheep 0.5853476615
# 64       f4a_water_avail 0.5722207703
# 65    f4a_cur_fastbreath 0.5678504092
# 66       f4a_cur_thirsty 0.5450300477
# 67          f4a_wash_eat 0.5436167320
# 68       f4a_house_radio 0.5269080068
# 69     f4a_drh_lessdrink 0.4927584118
# 70           f4a_ani_cat 0.4790536124
# 71     f4a_hometrt_othr1 0.4719691271
# 72            f4a_ani_no 0.4687964985
# 73        f4a_hometrt_ab 0.4426168775
# 74     f4a_hometrt_maize 0.4171617115
# 75     f4a_water_covwell 0.4116188455
# 76      f4a_hometrt_milk 0.3979342333
# 77       f4a_fuel_natgas 0.3873221299
# 78         f4a_fac_waste 0.3850577870
# 79        f4a_water_well 0.3843620478
# 80        f4a_seek_remdy 0.3730336904
# 81         f4a_ani_other 0.3718578951
# 82           f4a_ani_dog 0.3671253100
# 83       f4a_house_phone 0.3202664204
# 84          f4b_abn_hair 0.3094300969
# 85        f4a_house_none 0.3038627564
# 86      f4a_house_agland 0.2820429714
# 87        f4a_water_yard 0.2711498393
# 88        f4a_seek_pharm 0.2631595765
# 89         f4b_recommend 0.2520731450
# 90       f4b_chest_indrw 0.2297910458
# 91        f4b_skin_flaky 0.2292420795
# 92             f3_drh_iv 0.2187994399
# 93      f4a_drh_prolapse 0.1975286717
# 94      f4a_seek_privdoc 0.1900577710
# 95    f4a_water_covpwell 0.1874994874
# 96   f4a_hometrt_othrliq 0.1830255033
# 97     f4a_water_pubwell 0.1776399616
# 98       f4a_fuel_biogas 0.1754489154
# 99        f4a_seek_other 0.1697277284
# 100       f4a_trt_method 0.1651421077
# 101        f4a_trt_water 0.1589157187
# 102          f4a_ani_cow 0.1430671409
# 103            f4a_floor 0.1415405488
# 104    f4a_hometrt_othr2 0.1316339242
# 105      f4a_seek_friend 0.1154626809
# 106        f4a_wash_othr 0.1033847358
# 107          f4b_bipedal 0.0744400757
# 108      f4a_water_house 0.0646016826
# 109            f4b_admit 0.0616805820
# 110      f4a_wash_animal 0.0609110628
# 111          f3_drh_hosp 0.0450064586
# 112           f4b_rectal 0.0403367006
# 113       f4a_water_bore 0.0208629397
# 114        f4a_drh_consc 0.0100986885
# 115     f4a_fuel_propane 0.0100136412
# 116         f4a_seek_doc 0.0052901069
# 117  f4a_water_prospring 0.0044713349
# 118   f4a_water_deepwell 0.0039103541
# 119       f4a_house_boat 0.0029781642
# 120        f4a_fuel_coal 0.0022035213
# 121     f4a_hometrt_zinc 0.0009828394
# 122         f4a_drh_conv 0.0009351398
# 123                 site 0.0000000000
# 124       f4a_fuel_grass 0.0000000000
# 125        f4a_fuel_dung 0.0000000000
# 126        f4a_fuel_crop 0.0000000000
# 127        f4a_fuel_kero 0.0000000000
# 128   f4a_water_unspring 0.0000000000
# 129      f4a_water_river 0.0000000000
# 130       f4a_water_pond 0.0000000000
# 131       f4a_water_rain 0.0000000000
# 132  f4a_water_shallwell 0.0000000000
# 133       f4a_water_othr 0.0000000000

# AUC          SE     lower     upper level Model nvar
# 1 0.6294747 0.004792995 0.6200806 0.6388688  0.95    LR    2
# 2 0.6030409 0.005028438 0.5931853 0.6128964  0.95    LR    5
# 3 0.5721534 0.005454497 0.5614628 0.5828441  0.95    LR   10
# 4 0.6161968 0.005269247 0.6058693 0.6265243  0.95    RF    2
# 5 0.6650714 0.004936347 0.6553964 0.6747465  0.95    RF    5
# 6 0.6536583 0.004993425 0.6438713 0.6634452  0.95    RF   10

# nvar    intc intc_LCI intc_UCI slope slope_LCI slope_UCI
# <dbl>   <dbl>    <dbl>    <dbl> <dbl>     <dbl>     <dbl>
#   1     2 0.00918   -0.410    0.395 1.08     -0.177      2.33
# 2     5 0.0102    -0.410    0.396 0.854    -0.339      2.04
# 3    10 0.0146    -0.409    0.405 0.513    -0.322      1.37

####################
#shigella AFe>=0.5 only MUAC Mozambique only
################### var screening, AUC ####
names<-append(x=names, values=c("f4b_muac")) #add MUAC to list
names <- names[!names %in% c("f4b_haz")] #take HAZ off list

shigella0.5MUAC.Moz <- CPR.funct(data=Moz,outcome="shigella_afe0.5",iter=100,nvars_opts=c(2,5,10))
shigella0.5MUAC.Moz[["df_imps"]]
shigella0.5MUAC.Moz[["AUC_df"]]
shigella0.5MUAC.Moz[["calib"]]

# names      var_red
# 1               base_age 6.709180e+00
# 2          f4a_drh_blood 4.414815e+00
# 3               f4b_muac 3.581613e+00
# 4               f4b_resp 3.390990e+00
# 5          f4a_breastfed 3.031378e+00
# 6               f4b_temp 2.209460e+00
# 7          f4a_ppl_house 1.611049e+00
# 8          f4a_drh_vomit 1.516598e+00
# 9              f4b_admit 1.232497e+00
# 10          f4a_drh_days 1.125592e+00
# 11        f4a_disp_feces 1.042129e+00
# 12           f3_drh_hosp 9.825698e-01
# 13          f4a_dad_live 8.995046e-01
# 14         f4a_slp_rooms 8.742338e-01
# 15      f4a_yng_children 8.455611e-01
# 16         f4b_recommend 8.097498e-01
# 17         f4a_prim_schl 7.190606e-01
# 18       f4a_water_avail 6.832524e-01
# 19             f3_drh_iv 6.631651e-01
# 20         f4a_drh_cough 6.429754e-01
# 21             f4b_mouth 6.184932e-01
# 22        f4a_offr_drink 6.114234e-01
# 23        f4a_water_bore 6.108237e-01
# 24            f4b_mental 6.054933e-01
# 25       f4a_house_radio 5.894011e-01
# 26        f4a_max_stools 5.696734e-01
# 27              f4b_eyes 5.469227e-01
# 28          f4a_ms_water 5.452641e-01
# 29  f4a_drh_lethrgy_miss 5.289237e-01
# 30              f4b_skin 5.118924e-01
# 31       f4a_house_phone 5.045858e-01
# 32      f4a_water_pubtap 4.932110e-01
# 33      f4a_relationship 4.689347e-01
# 34           f4a_ani_cat 4.610535e-01
# 35      f4a_house_agland 4.429353e-01
# 36      f4a_cur_drymouth 4.253771e-01
# 37      f4a_house_fridge 4.099721e-01
# 38     f4a_water_pubwell 4.050811e-01
# 39     f4a_drh_bellypain 3.898472e-01
# 40          f4a_wash_use 3.859456e-01
# 41      f4a_seek_outside 3.821266e-01
# 42        f4a_wash_nurse 3.818524e-01
# 43        f4a_wash_child 3.805513e-01
# 44             f3_gender 3.740100e-01
# 45    f4a_cur_fastbreath 3.669869e-01
# 46          f4a_ani_goat 3.651216e-01
# 47          f4a_wash_eat 3.639085e-01
# 48       f4a_hometrt_ors 3.605887e-01
# 49          f4a_drh_conv 3.574695e-01
# 50        f4b_under_nutr 3.559630e-01
# 51         f4a_wash_cook 3.497197e-01
# 52          f4a_ani_fowl 3.349105e-01
# 53            f4a_ani_no 3.342447e-01
# 54      f4a_hometrt_none 3.320845e-01
# 55          f4a_cur_skin 3.317282e-01
# 56             f4a_floor 3.287704e-01
# 57        f4a_house_elec 3.234377e-01
# 58         f4a_fac_waste 3.206457e-01
# 59     f4a_fuel_charcoal 3.167533e-01
# 60       f4a_ani_rodents 3.152612e-01
# 61        f4a_house_bike 3.081217e-01
# 62     f4a_hometrt_othr2 3.035354e-01
# 63      f4a_hometrt_herb 2.865196e-01
# 64         f4a_house_car 2.755815e-01
# 65         f3_drh_turgor 2.734089e-01
# 66      f4a_cur_restless 2.678362e-01
# 67        f4a_drh_thirst 2.660049e-01
# 68        f4a_water_yard 2.591227e-01
# 69        f4a_house_tele 2.590097e-01
# 70       f4a_cur_thirsty 2.559276e-01
# 71           f4a_ani_dog 2.547841e-01
# 72          f4b_abn_hair 2.490032e-01
# 73      f4a_drh_restless 2.484223e-01
# 74     f4a_hometrt_othr1 2.411519e-01
# 75       f4b_chest_indrw 2.380268e-01
# 76        f4a_drh_strain 2.308602e-01
# 77           f4b_bipedal 2.214174e-01
# 78        f4a_trt_method 2.099276e-01
# 79          f4a_wash_def 2.082045e-01
# 80    f4a_water_covpwell 2.049738e-01
# 81     f4a_water_covwell 1.767135e-01
# 82         f4a_fuel_wood 1.653463e-01
# 83     f4a_drh_lessdrink 1.587104e-01
# 84         f4a_fuel_elec 1.573527e-01
# 85         f4a_trt_water 1.553366e-01
# 86           f4a_ani_cow 1.480055e-01
# 87         f4a_share_fac 1.143732e-01
# 88        f4a_water_rain 1.050997e-01
# 89        f4b_skin_flaky 9.226831e-02
# 90      f4a_drh_prolapse 7.776564e-02
# 91         f4a_ani_other 7.592323e-02
# 92       f4a_wash_animal 7.032749e-02
# 93        f4a_house_boat 6.796956e-02
# 94       f4a_fuel_natgas 5.980988e-02
# 95        f4a_house_cart 4.906125e-02
# 96        f4a_water_othr 4.754430e-02
# 97   f4a_hometrt_othrliq 4.604952e-02
# 98        f4a_water_well 4.330392e-02
# 99     f4a_hometrt_maize 4.228351e-02
# 100      f4a_store_water 4.067331e-02
# 101       f4a_house_none 3.018712e-02
# 102      f4a_house_scoot 2.668359e-02
# 103     f4a_water_bought 2.559959e-02
# 104        f4a_drh_consc 2.505468e-02
# 105      f4a_water_house 2.246797e-02
# 106   f4a_water_deepwell 1.530018e-02
# 107      f4a_seek_friend 1.090590e-02
# 108        f4a_ani_sheep 3.086195e-03
# 109       f4a_hometrt_ab 2.583952e-03
# 110       f4a_fuel_grass 1.648974e-03
# 111       f4a_fuel_other 1.454261e-03
# 112       f4a_water_pond 1.051669e-03
# 113           f4b_rectal 6.339706e-04
# 114      f4a_water_river 4.000000e-04
# 115     f4a_hometrt_zinc 3.354037e-04
# 116       f4a_seek_other 1.000000e-04
# 117     f4a_fuel_propane 3.308824e-05
# 118                 site 0.000000e+00
# 119      f4a_fuel_biogas 0.000000e+00
# 120        f4a_fuel_coal 0.000000e+00
# 121        f4a_fuel_dung 0.000000e+00
# 122        f4a_fuel_crop 0.000000e+00
# 123        f4a_fuel_kero 0.000000e+00
# 124  f4a_water_prospring 0.000000e+00
# 125   f4a_water_unspring 0.000000e+00
# 126  f4a_water_shallwell 0.000000e+00
# 127        f4a_wash_othr 0.000000e+00
# 128     f4a_hometrt_milk 0.000000e+00
# 129       f4a_seek_pharm 0.000000e+00
# 130      f4a_seek_healer 0.000000e+00
# 131         f4a_seek_doc 0.000000e+00
# 132     f4a_seek_privdoc 0.000000e+00
# 133       f4a_seek_remdy 0.000000e+00

# AUC          SE     lower     upper level Model nvar
# 1 0.8189387 0.005198175 0.8087504 0.8291269  0.95    LR    2
# 2 0.8171563 0.005731073 0.8059236 0.8283890  0.95    LR    5
# 3 0.8135040 0.005869763 0.8019995 0.8250086  0.95    LR   10
# 4 0.8209365 0.005540957 0.8100764 0.8317966  0.95    RF    2
# 5 0.8308213 0.005342062 0.8203510 0.8412915  0.95    RF    5
# 6 0.8519900 0.004926546 0.8423341 0.8616458  0.95    RF   10

# nvar     intc intc_LCI intc_UCI slope slope_LCI slope_UCI
# <dbl>    <dbl>    <dbl>    <dbl> <dbl>     <dbl>     <dbl>
#   1     2  0.0148    -0.617    0.592 0.961     0.465      1.53
# 2     5 -0.00459   -0.660    0.597 0.958     0.514      1.50
# 3    10 -0.0144    -0.687    0.604 0.845     0.446      1.33

####################
#shigella AFe>=0.5 only MUAC Kenya only
################### var screening, AUC ####
names<-append(x=names, values=c("f4b_muac")) #add MUAC to list
names <- names[!names %in% c("f4b_haz")] #take HAZ off list

shigella0.5MUAC.Kenya <- CPR.funct(data=Kenya,outcome="shigella_afe0.5",iter=100,nvars_opts=c(2,5,10))
shigella0.5MUAC.Kenya[["df_imps"]]
shigella0.5MUAC.Kenya[["AUC_df"]]
shigella0.5MUAC.Kenya[["calib"]]

# names      var_red
# 1               base_age 3.3549978916
# 2               f4b_resp 2.8100923474
# 3               f4b_muac 2.7758939878
# 4          f4a_drh_blood 2.6696938710
# 5               f4b_temp 2.5850439953
# 6          f4a_ppl_house 1.7898061270
# 7          f4a_share_fac 1.4536119456
# 8           f4a_drh_days 1.4302811989
# 9          f4a_prim_schl 1.0934178721
# 10        f4a_max_stools 1.0844612505
# 11        f4a_offr_drink 1.0708984841
# 12         f4a_drh_cough 0.8017013642
# 13         f4a_wash_othr 0.7802127478
# 14        f4a_trt_method 0.7756905438
# 15           f4a_ani_cow 0.7595671124
# 16      f4a_yng_children 0.7549435695
# 17      f4a_water_pubtap 0.7044067107
# 18          f4a_dad_live 0.6898027442
# 19         f4a_drh_vomit 0.6839825747
# 20         f4a_slp_rooms 0.6227132549
# 21      f4a_relationship 0.6179016790
# 22         f4a_breastfed 0.6121337171
# 23       f4a_hometrt_ors 0.5857847919
# 24           f4a_ani_cat 0.5851262167
# 25             f4b_mouth 0.5763618363
# 26           f4a_ani_dog 0.5723921512
# 27        f4b_under_nutr 0.5669481322
# 28            f4b_mental 0.5486622824
# 29     f4a_hometrt_maize 0.5484194089
# 30       f4a_house_phone 0.5395964191
# 31        f4a_disp_feces 0.5199809965
# 32          f4a_ms_water 0.5196804093
# 33          f4a_wash_eat 0.5155833214
# 34          f4a_wash_def 0.5082864333
# 35        f4a_seek_pharm 0.4932658350
# 36      f4a_cur_drymouth 0.4897849016
# 37       f4a_cur_thirsty 0.4837099246
# 38     f4a_drh_bellypain 0.4679491889
# 39             f3_gender 0.4657452337
# 40        f4a_drh_strain 0.4556241788
# 41         f4b_recommend 0.4520243304
# 42        f4a_house_bike 0.4479319533
# 43       f4a_ani_rodents 0.4477066960
# 44          f4a_wash_use 0.4342483503
# 45          f4a_ani_goat 0.4281401928
# 46        f4a_water_rain 0.4257166523
# 47     f4a_fuel_charcoal 0.4219699198
# 48              f4b_skin 0.4166079804
# 49  f4a_drh_lethrgy_miss 0.4110516691
# 50        f4a_drh_thirst 0.4099069048
# 51      f4a_cur_restless 0.4087381355
# 52         f4a_ani_sheep 0.4072798582
# 53      f4a_drh_restless 0.3995635691
# 54        f4a_wash_child 0.3986437417
# 55      f4a_hometrt_none 0.3952428387
# 56         f4a_trt_water 0.3912849216
# 57      f4a_seek_outside 0.3878129872
# 58       f4a_fuel_biogas 0.3853186065
# 59         f4a_wash_cook 0.3669728944
# 60       f4a_water_river 0.3557441750
# 61          f4b_abn_hair 0.3529547124
# 62        f4a_wash_nurse 0.3445590302
# 63         f3_drh_turgor 0.3301715958
# 64             f4a_floor 0.3290006245
# 65          f4a_ani_fowl 0.3286309815
# 66     f4a_hometrt_othr1 0.3273015279
# 67       f4a_house_radio 0.3213122952
# 68    f4a_cur_fastbreath 0.3096890416
# 69       f4a_seek_healer 0.3075337665
# 70     f4a_drh_lessdrink 0.3017845661
# 71    f4a_water_deepwell 0.2985470004
# 72        f4a_hometrt_ab 0.2957198929
# 73        f4a_water_bore 0.2945116198
# 74       f4a_store_water 0.2742944397
# 75         f4a_drh_consc 0.2678631110
# 76          f4a_cur_skin 0.2599187878
# 77        f4a_house_tele 0.2561923667
# 78       f4a_water_avail 0.2558433740
# 79      f4a_hometrt_herb 0.2529445600
# 80       f4a_wash_animal 0.2527321912
# 81        f4a_water_pond 0.2488929906
# 82         f4a_fuel_wood 0.2386644407
# 83      f4a_house_agland 0.2367715285
# 84         f4a_fac_waste 0.2286369238
# 85        f4a_house_cart 0.2273836760
# 86        f4a_water_well 0.2177649994
# 87             f3_drh_iv 0.2130925775
# 88        f4a_house_none 0.2001126214
# 89           f4b_bipedal 0.1941335829
# 90        f4b_skin_flaky 0.1915164565
# 91           f3_drh_hosp 0.1866446155
# 92    f4a_water_unspring 0.1825856966
# 93        f4a_seek_other 0.1805010568
# 94         f4a_fuel_crop 0.1711096632
# 95      f4a_drh_prolapse 0.1671066325
# 96         f4a_ani_other 0.1632919003
# 97          f4a_drh_conv 0.1511366109
# 98        f4a_house_elec 0.1437896788
# 99             f4b_admit 0.1414519133
# 100           f4b_rectal 0.1411350227
# 101      f4a_house_scoot 0.1370268290
# 102        f4a_fuel_kero 0.1311995421
# 103    f4a_water_pubwell 0.1282207949
# 104  f4a_water_shallwell 0.1033334652
# 105    f4a_hometrt_othr2 0.0957790717
# 106           f4a_ani_no 0.0791366996
# 107         f4a_seek_doc 0.0683048846
# 108       f4a_water_yard 0.0656506615
# 109       f4a_fuel_grass 0.0653394275
# 110             f4b_eyes 0.0639433928
# 111     f4a_hometrt_zinc 0.0639067014
# 112  f4a_water_prospring 0.0448796149
# 113       f4a_seek_remdy 0.0287595413
# 114      f4b_chest_indrw 0.0109463376
# 115  f4a_hometrt_othrliq 0.0094754540
# 116   f4a_water_covpwell 0.0070058758
# 117     f4a_hometrt_milk 0.0057960014
# 118        f4a_house_car 0.0038544833
# 119      f4a_seek_friend 0.0027423077
# 120      f4a_water_house 0.0012181287
# 121     f4a_house_fridge 0.0011428571
# 122      f4a_fuel_natgas 0.0007111111
# 123    f4a_water_covwell 0.0001469139
# 124                 site 0.0000000000
# 125       f4a_house_boat 0.0000000000
# 126        f4a_fuel_elec 0.0000000000
# 127     f4a_fuel_propane 0.0000000000
# 128        f4a_fuel_coal 0.0000000000
# 129        f4a_fuel_dung 0.0000000000
# 130       f4a_fuel_other 0.0000000000
# 131     f4a_water_bought 0.0000000000
# 132       f4a_water_othr 0.0000000000
# 133     f4a_seek_privdoc 0.0000000000

# AUC          SE     lower     upper level Model nvar
# 1 0.6148808 0.006970174 0.6012195 0.6285421  0.95    LR    2
# 2 0.6790027 0.007452718 0.6643956 0.6936097  0.95    LR    5
# 3 0.6667509 0.007476920 0.6520964 0.6814054  0.95    LR   10
# 4 0.5660948 0.007025022 0.5523260 0.5798636  0.95    RF    2
# 5 0.5900477 0.007615777 0.5751211 0.6049744  0.95    RF    5
# 6 0.6149701 0.008129173 0.5990372 0.6309030  0.95    RF   10

# nvar    intc intc_LCI intc_UCI slope slope_LCI slope_UCI
# <dbl>   <dbl>    <dbl>    <dbl> <dbl>     <dbl>     <dbl>
#   1     2 -0.0309   -0.604    0.474 0.901    -0.194      2.04
# 2     5 -0.0220   -0.611    0.501 0.926     0.261      1.63
# 3    10 -0.0297   -0.626    0.501 0.657     0.116      1.24

####################
#shigella AFe>=0.5 only MUAC India only
################### var screening, AUC ####
names<-append(x=names, values=c("f4b_muac")) #add MUAC to list
names <- names[!names %in% c("f4b_haz")] #take HAZ off list

summary(India$f4a_disp_feces)
# 1   2   3   4   6 
# 1   0 197 609  24 
India_complete <- India
India_complete$f4a_disp_feces <- as.factor(as.numeric(ifelse(India_complete$f4a_disp_feces==1 | India_complete$f4a_disp_feces==6,6,India_complete$f4a_disp_feces)))
summary(India_complete$f4a_disp_feces)
table(India_complete$f4a_disp_feces)

shigella0.5MUAC.India <- CPR.funct(data=India_complete,outcome="shigella_afe0.5",iter=100,nvars_opts=c(2,5,10))
shigella0.5MUAC.India[["df_imps"]]
shigella0.5MUAC.India[["AUC_df"]]
shigella0.5MUAC.India[["calib"]]

# names      var_red
# 1          f4a_drh_blood 1.209865e+01
# 2               base_age 9.170974e+00
# 3               f4b_muac 5.457868e+00
# 4               f4b_resp 5.169854e+00
# 5               f4b_temp 4.611209e+00
# 6          f4a_share_fac 4.453055e+00
# 7          f4a_ppl_house 3.647720e+00
# 8           f4a_drh_days 2.329913e+00
# 9          f4a_prim_schl 2.233619e+00
# 10        f4a_disp_feces 1.823665e+00
# 11              f4b_eyes 1.628991e+00
# 12        f4a_max_stools 1.624034e+00
# 13         f4a_slp_rooms 1.553071e+00
# 14        f4a_offr_drink 1.466536e+00
# 15        f4a_trt_method 1.425752e+00
# 16         f4a_drh_vomit 1.383632e+00
# 17      f4a_yng_children 1.345078e+00
# 18     f4a_drh_bellypain 1.196301e+00
# 19          f4a_wash_eat 1.173465e+00
# 20            f4b_mental 1.166368e+00
# 21         f4a_drh_cough 1.147781e+00
# 22       f4a_house_radio 1.102121e+00
# 23       f4a_cur_thirsty 1.084468e+00
# 24         f4a_breastfed 1.069211e+00
# 25        f4a_wash_child 9.900556e-01
# 26             f3_gender 9.568297e-01
# 27     f4a_hometrt_maize 9.214548e-01
# 28          f4a_ani_fowl 9.166637e-01
# 29        f4a_wash_nurse 9.028373e-01
# 30         f4a_fuel_coal 9.001560e-01
# 31          f4a_wash_use 8.875829e-01
# 32    f4a_cur_fastbreath 8.713636e-01
# 33           f4a_ani_cat 8.706951e-01
# 34      f4a_fuel_propane 8.612767e-01
# 35          f4a_cur_skin 8.554520e-01
# 36         f4a_trt_water 8.285016e-01
# 37         f4b_recommend 8.158931e-01
# 38          f4a_wash_def 8.138564e-01
# 39       f4a_house_phone 8.070219e-01
# 40             f4b_mouth 7.999542e-01
# 41         f4a_wash_cook 7.936963e-01
# 42         f4a_fuel_wood 7.919703e-01
# 43      f4a_hometrt_none 7.811355e-01
# 44      f4a_cur_restless 7.778519e-01
# 45        f4a_drh_thirst 7.712362e-01
# 46      f4a_water_pubtap 7.601688e-01
# 47        f4a_house_tele 7.371620e-01
# 48      f4a_drh_restless 7.297104e-01
# 49     f4a_drh_lessdrink 7.180885e-01
# 50       f4a_hometrt_ors 7.105445e-01
# 51           f4a_ani_dog 7.075202e-01
# 52        f4a_house_bike 7.072020e-01
# 53          f4a_ms_water 6.758052e-01
# 54  f4a_drh_lethrgy_miss 6.637874e-01
# 55        f4a_water_yard 6.580353e-01
# 56      f4a_relationship 6.275462e-01
# 57         f4a_fac_waste 6.003261e-01
# 58         f4a_fuel_kero 5.938841e-01
# 59             f4b_admit 5.672675e-01
# 60      f4a_seek_outside 5.568263e-01
# 61           f3_drh_hosp 5.328202e-01
# 62        f4b_under_nutr 5.283653e-01
# 63       f4a_house_scoot 5.056396e-01
# 64     f4a_hometrt_othr1 4.724351e-01
# 65        f4a_drh_strain 4.701945e-01
# 66           f4a_ani_cow 4.588045e-01
# 67      f4a_house_fridge 4.240082e-01
# 68          f4a_ani_goat 4.053966e-01
# 69         f3_drh_turgor 4.039969e-01
# 70    f4a_water_deepwell 3.970510e-01
# 71      f4a_cur_drymouth 3.945689e-01
# 72      f4a_seek_privdoc 3.929830e-01
# 73         f4a_ani_other 3.914016e-01
# 74      f4a_hometrt_milk 3.802061e-01
# 75              f4b_skin 3.607467e-01
# 76        f4a_hometrt_ab 3.584808e-01
# 77             f3_drh_iv 3.262506e-01
# 78       f4a_wash_animal 3.194617e-01
# 79          f4a_seek_doc 3.012258e-01
# 80             f4a_floor 2.931777e-01
# 81       f4a_ani_rodents 2.725774e-01
# 82        f4a_house_none 2.674242e-01
# 83      f4a_drh_prolapse 2.671043e-01
# 84       f4a_store_water 2.406874e-01
# 85        f4a_house_elec 2.321911e-01
# 86        f4a_seek_other 2.178878e-01
# 87       f4a_water_house 1.906823e-01
# 88        f4a_water_well 1.763408e-01
# 89      f4a_hometrt_herb 1.704205e-01
# 90       f4a_water_avail 1.681678e-01
# 91        f4a_seek_pharm 1.647639e-01
# 92            f4b_rectal 1.474742e-01
# 93     f4a_hometrt_othr2 1.337709e-01
# 94       f4b_chest_indrw 1.207314e-01
# 95         f4a_house_car 1.168267e-01
# 96         f4a_wash_othr 1.153541e-01
# 97          f4a_dad_live 1.138087e-01
# 98            f4a_ani_no 7.082293e-02
# 99          f4b_abn_hair 5.659165e-02
# 100       f4a_fuel_grass 5.370013e-02
# 101        f4a_fuel_elec 5.260574e-02
# 102       f4a_water_othr 5.258288e-02
# 103     f4a_hometrt_zinc 4.638096e-02
# 104  f4a_hometrt_othrliq 4.583967e-02
# 105         f4a_drh_conv 3.063045e-02
# 106    f4a_water_covwell 3.003726e-02
# 107      f4a_fuel_natgas 2.785264e-02
# 108    f4a_fuel_charcoal 2.006478e-02
# 109  f4a_water_shallwell 1.372528e-02
# 110       f4a_seek_remdy 3.811378e-03
# 111     f4a_house_agland 3.415440e-03
# 112     f4a_water_bought 3.038186e-03
# 113      f4a_seek_friend 2.580478e-03
# 114        f4a_ani_sheep 7.111111e-04
# 115        f4a_fuel_dung 6.240494e-04
# 116                 site 0.000000e+00
# 117       f4a_house_cart 0.000000e+00
# 118       f4a_house_boat 0.000000e+00
# 119      f4a_fuel_biogas 0.000000e+00
# 120        f4a_fuel_crop 0.000000e+00
# 121       f4a_fuel_other 0.000000e+00
# 122   f4a_water_covpwell 0.000000e+00
# 123  f4a_water_prospring 0.000000e+00
# 124   f4a_water_unspring 0.000000e+00
# 125    f4a_water_pubwell 0.000000e+00
# 126      f4a_water_river 0.000000e+00
# 127       f4a_water_pond 0.000000e+00
# 128       f4a_water_rain 0.000000e+00
# 129       f4a_water_bore 0.000000e+00
# 130        f4a_drh_consc 0.000000e+00
# 131      f4a_seek_healer 0.000000e+00
# 132          f4b_bipedal 0.000000e+00
# 133       f4b_skin_flaky 0.000000e+00

# AUC          SE     lower     upper level Model nvar
# 1 0.8042034 0.004433769 0.7955134 0.8128934  0.95    LR    2
# 2 0.7932302 0.004599968 0.7842145 0.8022460  0.95    LR    5
# 3 0.7747260 0.005098984 0.7647322 0.7847198  0.95    LR   10
# 4 0.7964279 0.004649952 0.7873142 0.8055417  0.95    RF    2
# 5 0.7869301 0.004554500 0.7780035 0.7958568  0.95    RF    5
# 6 0.8070417 0.004387937 0.7984415 0.8156419  0.95    RF   10

# nvar    intc intc_LCI intc_UCI slope slope_LCI slope_UCI
# <dbl>   <dbl>    <dbl>    <dbl> <dbl>     <dbl>     <dbl>
#   1     2 -0.0899   -0.560    0.347 1.04      0.673      1.46
# 2     5 -0.0883   -0.560    0.351 1.01      0.651      1.42
# 3    10 -0.0945   -0.571    0.350 0.864     0.528      1.24


####################
#shigella AFe>=0.5 only MUAC Bangladesh only
################### var screening, AUC ####
names<-append(x=names, values=c("f4b_muac")) #add MUAC to list
names <- names[!names %in% c("f4b_haz")] #take HAZ off list

shigella0.5MUAC.Bdesh <- CPR.funct(data=Bdesh,outcome="shigella_afe0.5",iter=100,nvars_opts=c(2,5,10))
shigella0.5MUAC.Bdesh[["df_imps"]]
shigella0.5MUAC.Bdesh[["AUC_df"]]
shigella0.5MUAC.Bdesh[["calib"]]

# Pak <- cases %>% filter(site == "Pakistan")

# names      var_red
# 1               base_age 45.577269268
# 2               f4b_temp  8.330509740
# 3               f4b_muac  8.201859297
# 4          f4a_drh_blood  7.970851729
# 5               f4b_resp  7.427337751
# 6               f4b_eyes  5.695252871
# 7         f4a_drh_strain  3.851761911
# 8           f4a_drh_days  3.730951486
# 9          f4a_ppl_house  3.538972495
# 10     f4a_drh_bellypain  3.487280723
# 11         f4a_drh_vomit  3.345739370
# 12         f4a_slp_rooms  2.836990566
# 13         f4a_fac_waste  2.641428828
# 14         f4a_breastfed  2.563365835
# 15        f4a_max_stools  2.549762022
# 16         f4a_share_fac  2.173024504
# 17         f4a_prim_schl  2.124198859
# 18             f4b_admit  1.881340666
# 19         f4b_recommend  1.846553289
# 20      f4a_hometrt_none  1.820093006
# 21          f4a_dad_live  1.718253038
# 22             f3_drh_iv  1.704935150
# 23        f4a_drh_thirst  1.674617382
# 24         f4a_drh_cough  1.666612104
# 25       f4a_cur_thirsty  1.603573573
# 26        f4a_offr_drink  1.550710883
# 27       f4a_hometrt_ors  1.379996450
# 28        f4a_disp_feces  1.370535431
# 29      f4a_seek_outside  1.340755235
# 30      f4a_yng_children  1.325538383
# 31      f4a_house_agland  1.281802125
# 32        f4a_wash_child  1.202065641
# 33        f4a_house_elec  1.183116082
# 34           f3_drh_hosp  1.182454352
# 35             f3_gender  1.153481125
# 36        f4a_hometrt_ab  1.151080414
# 37        f4a_house_tele  1.146668843
# 38          f4a_ani_fowl  1.128229572
# 39          f4a_wash_use  1.127533393
# 40         f4a_fuel_dung  1.111964820
# 41          f4a_drh_conv  1.096863893
# 42      f4a_hometrt_zinc  1.096233893
# 43        f4a_house_bike  1.094633965
# 44         f4a_wash_cook  1.092847557
# 45          f4a_wash_eat  1.090153177
# 46       f4a_house_phone  1.088320398
# 47        f4a_seek_pharm  1.085755318
# 48      f4a_drh_restless  1.042782157
# 49       f4a_house_radio  1.025869871
# 50    f4a_water_deepwell  1.020899196
# 51         f4a_fuel_crop  1.014235011
# 52     f4a_hometrt_othr1  1.002294300
# 53         f4a_fuel_wood  0.998942511
# 54   f4a_water_shallwell  0.987804045
# 55      f4a_house_fridge  0.987508865
# 56           f4a_ani_cow  0.973493273
# 57      f4a_cur_drymouth  0.935425258
# 58       f4a_wash_animal  0.934670756
# 59      f4a_cur_restless  0.922328346
# 60             f4a_floor  0.920308810
# 61        f4a_fuel_grass  0.911082504
# 62        f4a_wash_nurse  0.909839193
# 63         f4a_wash_othr  0.905465383
# 64          f4a_wash_def  0.872727116
# 65           f4a_ani_dog  0.847647125
# 66  f4a_drh_lethrgy_miss  0.822762977
# 67       f4a_ani_rodents  0.791407557
# 68           f4a_ani_cat  0.791277045
# 69       f4a_house_scoot  0.781702880
# 70          f4a_ani_goat  0.716586737
# 71         f3_drh_turgor  0.707514190
# 72          f4a_cur_skin  0.680571992
# 73        f4a_seek_remdy  0.679433799
# 74            f4a_ani_no  0.570357195
# 75        f4a_house_none  0.562198682
# 76       f4a_fuel_natgas  0.480274375
# 77          f4a_seek_doc  0.453485218
# 78            f4b_mental  0.430566627
# 79      f4a_hometrt_herb  0.424068571
# 80              f4b_skin  0.397961698
# 81     f4a_drh_lessdrink  0.369975208
# 82        f4a_trt_method  0.361434034
# 83         f4a_trt_water  0.340138986
# 84         f4a_ani_sheep  0.324173013
# 85    f4a_cur_fastbreath  0.302606560
# 86             f4b_mouth  0.300047747
# 87         f4a_house_car  0.267115946
# 88         f4a_drh_consc  0.250741119
# 89        f4a_house_boat  0.236963319
# 90     f4a_hometrt_othr2  0.211286238
# 91      f4a_seek_privdoc  0.194088788
# 92     f4a_hometrt_maize  0.187574056
# 93      f4a_drh_prolapse  0.161926905
# 94       f4a_seek_healer  0.148267436
# 95         f4a_fuel_elec  0.142385124
# 96        f4b_under_nutr  0.139866396
# 97      f4a_fuel_propane  0.120133050
# 98   f4a_hometrt_othrliq  0.098359572
# 99        f4a_house_cart  0.081592142
# 100      f4a_store_water  0.081341077
# 101       f4a_water_well  0.079154708
# 102         f4a_ms_water  0.077125377
# 103       f4a_fuel_other  0.055468890
# 104     f4a_relationship  0.054102954
# 105      f4a_seek_friend  0.033029163
# 106      f4b_chest_indrw  0.018629164
# 107    f4a_water_covwell  0.017350934
# 108           f4b_rectal  0.016419901
# 109      f4a_water_avail  0.008504195
# 110    f4a_water_pubwell  0.002353338
# 111        f4a_fuel_kero  0.002172997
# 112      f4a_water_house  0.001466351
# 113        f4a_ani_other  0.001310224
# 114                 site  0.000000000
# 115      f4a_fuel_biogas  0.000000000
# 116        f4a_fuel_coal  0.000000000
# 117    f4a_fuel_charcoal  0.000000000
# 118       f4a_water_yard  0.000000000
# 119   f4a_water_covpwell  0.000000000
# 120     f4a_water_pubtap  0.000000000
# 121  f4a_water_prospring  0.000000000
# 122   f4a_water_unspring  0.000000000
# 123      f4a_water_river  0.000000000
# 124       f4a_water_pond  0.000000000
# 125       f4a_water_rain  0.000000000
# 126     f4a_water_bought  0.000000000
# 127       f4a_water_othr  0.000000000
# 128       f4a_water_bore  0.000000000
# 129     f4a_hometrt_milk  0.000000000
# 130       f4a_seek_other  0.000000000
# 131          f4b_bipedal  0.000000000
# 132         f4b_abn_hair  0.000000000
# 133       f4b_skin_flaky  0.000000000

# AUC          SE     lower     upper level Model nvar
# 1 0.8503952 0.003113447 0.8442929 0.8564974  0.95    LR    2
# 2 0.8978208 0.002530504 0.8928612 0.9027805  0.95    LR    5
# 3 0.9058597 0.002417219 0.9011221 0.9105974  0.95    LR   10
# 4 0.8375491 0.003133531 0.8314075 0.8436907  0.95    RF    2
# 5 0.9116285 0.002238055 0.9072420 0.9160150  0.95    RF    5
# 6 0.9210763 0.002134209 0.9168934 0.9252593  0.95    RF   10

# nvar      intc intc_LCI intc_UCI slope slope_LCI slope_UCI
# <dbl>     <dbl>    <dbl>    <dbl> <dbl>     <dbl>     <dbl>
#   1     2 -0.0154     -0.381    0.351  1.03     0.712      1.41
# 2     5 -0.0148     -0.423    0.394  1.05     0.756      1.41
# 3    10  0.000332   -0.418    0.420  1.01     0.727      1.37

####################
#shigella AFe>=0.5 only MUAC Pakistan only
################### var screening, AUC ####
names<-append(x=names, values=c("f4b_muac")) #add MUAC to list
names <- names[!names %in% c("f4b_haz")] #take HAZ off list

shigella0.5MUAC.Pak <- CPR.funct(data=Pak,outcome="shigella_afe0.5",iter=100,nvars_opts=c(2,5,10))
shigella0.5MUAC.Pak[["df_imps"]]
shigella0.5MUAC.Pak[["AUC_df"]]
shigella0.5MUAC.Pak[["calib"]]

# names      var_red
# 1          f4a_drh_blood 13.921425043
# 2               base_age  7.880024858
# 3               f4b_muac  6.641009001
# 4               f4b_resp  5.977212952
# 5               f4b_temp  5.298615234
# 6          f4a_ppl_house  4.431700958
# 7           f4a_drh_days  3.126995069
# 8       f4a_yng_children  3.000581441
# 9          f4b_recommend  2.942218913
# 10       f4a_water_avail  2.740298968
# 11              f4b_eyes  2.465313985
# 12         f4a_drh_vomit  2.257324821
# 13         f4a_prim_schl  2.208528098
# 14         f4a_slp_rooms  2.197977619
# 15        f4a_max_stools  2.057455153
# 16         f4a_breastfed  1.803174141
# 17        f4a_trt_method  1.788808909
# 18        f4a_offr_drink  1.643002424
# 19         f4a_share_fac  1.514647364
# 20          f4a_ms_water  1.392661752
# 21             f4b_mouth  1.367261186
# 22            f4b_mental  1.346102361
# 23        f4a_house_tele  1.255975542
# 24     f4a_drh_bellypain  1.225559441
# 25             f3_drh_iv  1.154693424
# 26        f4a_disp_feces  1.153854597
# 27      f4a_relationship  1.145437968
# 28         f4a_drh_cough  1.127685560
# 29      f4a_house_fridge  1.115700565
# 30             f3_gender  1.096401007
# 31  f4a_drh_lethrgy_miss  1.075620947
# 32          f4a_cur_skin  1.074026117
# 33         f4a_wash_cook  1.049975017
# 34        f4a_drh_thirst  1.033639545
# 35      f4a_cur_drymouth  1.032928234
# 36              f4b_skin  1.028258508
# 37        f4a_wash_nurse  1.010699775
# 38        f4a_wash_child  1.003259518
# 39         f3_drh_turgor  0.989806168
# 40       f4a_hometrt_ors  0.968696829
# 41       f4a_water_house  0.964307362
# 42     f4a_drh_lessdrink  0.957363472
# 43         f4a_fac_waste  0.955227298
# 44        f4b_under_nutr  0.951649599
# 45      f4a_cur_restless  0.945844678
# 46       f4a_cur_thirsty  0.942774054
# 47      f4a_seek_outside  0.942144982
# 48      f4a_hometrt_none  0.936618427
# 49       f4a_house_phone  0.929805749
# 50       f4a_store_water  0.904111555
# 51          f4a_wash_def  0.902635527
# 52      f4a_water_bought  0.899037611
# 53        f4a_drh_strain  0.865413184
# 54         f4a_trt_water  0.858118454
# 55          f4a_wash_use  0.834844274
# 56      f4a_drh_restless  0.831999996
# 57          f4a_wash_eat  0.816596527
# 58       f4a_wash_animal  0.784964149
# 59            f4a_ani_no  0.780610274
# 60       f4a_house_scoot  0.764296731
# 61             f4a_floor  0.745064986
# 62     f4a_hometrt_maize  0.727495485
# 63         f4a_drh_consc  0.706756740
# 64          f4a_ani_fowl  0.701918991
# 65     f4a_hometrt_othr1  0.694547323
# 66       f4a_fuel_natgas  0.681012393
# 67          f4a_dad_live  0.664186949
# 68         f4a_fuel_wood  0.628719413
# 69        f4a_water_yard  0.576589192
# 70      f4a_seek_privdoc  0.565184640
# 71          f4a_ani_goat  0.563069150
# 72        f4a_house_bike  0.528212628
# 73    f4a_cur_fastbreath  0.525624819
# 74      f4a_water_pubtap  0.521390564
# 75        f4a_seek_pharm  0.491965667
# 76       f4a_house_radio  0.471897381
# 77        f4a_hometrt_ab  0.461133943
# 78           f4a_ani_cat  0.455450224
# 79      f4a_hometrt_zinc  0.442212441
# 80          f4a_seek_doc  0.423691888
# 81         f4a_wash_othr  0.409593660
# 82        f4a_water_othr  0.401206477
# 83     f4a_hometrt_othr2  0.372781679
# 84           f4a_ani_cow  0.360411779
# 85         f4a_house_car  0.354137464
# 86          f4b_abn_hair  0.344904942
# 87           f4a_ani_dog  0.304123051
# 88         f4a_fuel_kero  0.274993723
# 89      f4a_house_agland  0.261807988
# 90        f4a_fuel_grass  0.261020557
# 91   f4a_hometrt_othrliq  0.258371150
# 92        f4a_house_boat  0.255249719
# 93      f4a_hometrt_milk  0.250806011
# 94           f3_drh_hosp  0.249151512
# 95             f4b_admit  0.217882509
# 96        f4b_skin_flaky  0.209779034
# 97        f4a_house_elec  0.195469573
# 98        f4a_seek_other  0.193844587
# 99         f4a_ani_other  0.189257776
# 100        f4a_fuel_dung  0.179569677
# 101      f4b_chest_indrw  0.166902377
# 102     f4a_hometrt_herb  0.153102220
# 103       f4a_house_cart  0.149562985
# 104     f4a_drh_prolapse  0.145711069
# 105       f4a_house_none  0.089057000
# 106       f4a_water_pond  0.087205358
# 107           f4b_rectal  0.083629422
# 108       f4a_seek_remdy  0.066355746
# 109       f4a_water_bore  0.047937210
# 110         f4a_drh_conv  0.038343999
# 111      f4a_ani_rodents  0.036869892
# 112      f4a_seek_healer  0.036280293
# 113   f4a_water_covpwell  0.027053093
# 114        f4a_ani_sheep  0.024803449
# 115        f4a_fuel_elec  0.010468612
# 116     f4a_fuel_propane  0.005042283
# 117      f4a_fuel_biogas  0.003619444
# 118          f4b_bipedal  0.003071558
# 119    f4a_water_covwell  0.001197833
# 120                 site  0.000000000
# 121        f4a_fuel_coal  0.000000000
# 122    f4a_fuel_charcoal  0.000000000
# 123        f4a_fuel_crop  0.000000000
# 124       f4a_fuel_other  0.000000000
# 125  f4a_water_prospring  0.000000000
# 126       f4a_water_well  0.000000000
# 127   f4a_water_unspring  0.000000000
# 128    f4a_water_pubwell  0.000000000
# 129      f4a_water_river  0.000000000
# 130   f4a_water_deepwell  0.000000000
# 131       f4a_water_rain  0.000000000
# 132  f4a_water_shallwell  0.000000000
# 133      f4a_seek_friend  0.000000000

# AUC          SE     lower     upper level Model nvar
# 1 0.7386650 0.004571672 0.7297046 0.7476253  0.95    LR    2
# 2 0.7109034 0.004878651 0.7013414 0.7204653  0.95    LR    5
# 3 0.7175361 0.004907864 0.7079169 0.7271553  0.95    LR   10
# 4 0.7535636 0.004443313 0.7448548 0.7622723  0.95    RF    2
# 5 0.7348203 0.004703145 0.7256023 0.7440383  0.95    RF    5
# 6 0.7379715 0.004696006 0.7287675 0.7471755  0.95    RF   10

# nvar    intc intc_LCI intc_UCI slope slope_LCI slope_UCI
# <dbl>   <dbl>    <dbl>    <dbl> <dbl>     <dbl>     <dbl>
#   1     2 -0.0182   -0.428    0.371 0.961     0.586      1.37
# 2     5 -0.0192   -0.430    0.372 0.931     0.564      1.33
# 3    10 -0.0191   -0.435    0.377 0.834     0.490      1.21

####################
#shigella AFe>=0.5 only MUAC by country - Africa, val in Asia (not including Bdesh)
################### var screening, AUC ####
names<-append(x=names, values=c("f4b_muac")) #add MUAC to list
names <- names[!names %in% c("f4b_haz")] #take HAZ off list

shigella0.5MUAC.Afr <- CPR.funct(data=Afr,outcome="shigella_afe0.5",iter=100,nvars_opts=c(2,5,10))
shigella0.5MUAC.Afr[["df_imps"]]
shigella0.5MUAC.Afr[["AUC_df"]]
shigella0.5MUAC.Afr[["calib"]]

# names      var_red
# 1               base_age 26.854273731
# 2          f4a_drh_blood 21.956262549
# 3               f4b_resp 14.743321818
# 4               f4b_muac 14.624948791
# 5               f4b_temp 13.672218699
# 6          f4a_ppl_house 11.981312631
# 7          f4a_slp_rooms  8.222368867
# 8           f4a_drh_days  7.593206578
# 9       f4a_yng_children  6.856329631
# 10         f4a_breastfed  6.432585917
# 11         f4a_share_fac  5.094109627
# 12         f4a_prim_schl  5.031292041
# 13         f4b_recommend  4.514271822
# 14        f4a_offr_drink  4.500681148
# 15         f4a_drh_vomit  4.350247921
# 16          f4a_dad_live  3.925262403
# 17              f4b_eyes  3.742335480
# 18                  site  3.577097407
# 19        f4a_max_stools  3.433710408
# 20          f4a_ms_water  3.178913416
# 21             f4b_mouth  3.177585927
# 22         f4a_drh_cough  2.829809131
# 23        f4a_disp_feces  2.800333808
# 24      f4a_relationship  2.679759799
# 25              f4b_skin  2.606468246
# 26      f4a_cur_drymouth  2.522849571
# 27       f4a_water_avail  2.511207358
# 28        f4a_trt_method  2.444280807
# 29     f4a_drh_bellypain  2.416892747
# 30        f4a_house_bike  2.391774479
# 31         f4a_wash_cook  2.335854128
# 32  f4a_drh_lethrgy_miss  2.316605028
# 33             f3_gender  2.274249701
# 34        f4a_wash_nurse  2.202880530
# 35      f4a_water_pubtap  2.192636141
# 36          f4a_wash_use  2.170546965
# 37            f4b_mental  2.146752860
# 38       f4a_ani_rodents  2.128242031
# 39          f4a_ani_fowl  2.109408364
# 40     f4a_fuel_charcoal  2.109263691
# 41        f4a_wash_child  2.102485298
# 42      f4a_house_agland  2.100858286
# 43          f4a_wash_def  2.081446632
# 44       f4a_hometrt_ors  2.050855008
# 45       f4a_cur_thirsty  2.039341405
# 46             f4b_admit  2.035341328
# 47           f4a_ani_cat  2.028895600
# 48      f4a_drh_restless  2.008692745
# 49      f4a_hometrt_none  2.006219010
# 50        f4a_drh_thirst  1.989903327
# 51          f4a_ani_goat  1.983512343
# 52         f4a_ani_sheep  1.965060277
# 53           f3_drh_hosp  1.938471028
# 54        f4a_house_elec  1.902794723
# 55        f4a_house_tele  1.878873822
# 56        f4a_drh_strain  1.871798455
# 57        f4a_house_cart  1.865632480
# 58        f4b_under_nutr  1.865103299
# 59       f4a_house_scoot  1.862757660
# 60         f4a_ani_other  1.860258580
# 61           f4a_ani_dog  1.841934002
# 62    f4a_cur_fastbreath  1.817515063
# 63     f4a_hometrt_maize  1.808667381
# 64      f4a_cur_restless  1.805537672
# 65           f4a_ani_cow  1.789770592
# 66      f4a_house_fridge  1.760416758
# 67       f4a_store_water  1.727865719
# 68      f4a_seek_outside  1.705161801
# 69       f4a_house_radio  1.687351427
# 70          f4a_cur_skin  1.644121643
# 71      f4a_hometrt_herb  1.614692140
# 72         f4a_house_car  1.594936191
# 73         f4a_fuel_wood  1.593573791
# 74          f4a_wash_eat  1.566049316
# 75         f3_drh_turgor  1.564479665
# 76       f4a_house_phone  1.554073665
# 77     f4a_water_pubwell  1.526956325
# 78             f3_drh_iv  1.522272291
# 79         f4a_trt_water  1.512272956
# 80         f4a_fac_waste  1.483835497
# 81     f4a_hometrt_othr1  1.471740354
# 82             f4a_floor  1.444051373
# 83     f4a_drh_lessdrink  1.439941867
# 84       f4a_seek_healer  1.408916570
# 85        f4a_water_bore  1.279378336
# 86        f4a_seek_pharm  1.248598689
# 87         f4a_wash_othr  1.200642920
# 88         f4a_fuel_elec  1.185269723
# 89        f4a_hometrt_ab  1.141552006
# 90        f4a_water_yard  1.079074414
# 91        f4a_water_well  1.050117296
# 92          f4b_abn_hair  1.002766013
# 93      f4a_water_bought  0.998152642
# 94            f4a_ani_no  0.987108260
# 95    f4a_water_deepwell  0.955903988
# 96        f4b_skin_flaky  0.931009223
# 97          f4a_drh_conv  0.926307794
# 98     f4a_water_covwell  0.880948147
# 99        f4a_fuel_other  0.855616777
# 100       f4a_water_rain  0.805734293
# 101      f4a_wash_animal  0.718769868
# 102      f4b_chest_indrw  0.716935186
# 103      f4a_fuel_biogas  0.715973088
# 104   f4a_water_covpwell  0.654292379
# 105       f4a_seek_remdy  0.612055973
# 106     f4a_drh_prolapse  0.584523266
# 107     f4a_hometrt_milk  0.582995922
# 108      f4a_water_river  0.575534265
# 109    f4a_hometrt_othr2  0.551970885
# 110      f4a_fuel_natgas  0.509737053
# 111        f4a_drh_consc  0.450090698
# 112          f4b_bipedal  0.448246698
# 113       f4a_house_none  0.428294076
# 114       f4a_seek_other  0.390027282
# 115       f4a_water_pond  0.377378012
# 116  f4a_hometrt_othrliq  0.330624569
# 117       f4a_water_othr  0.303518370
# 118   f4a_water_unspring  0.276996970
# 119        f4a_fuel_crop  0.275365008
# 120           f4b_rectal  0.257568073
# 121      f4a_seek_friend  0.245365337
# 122     f4a_seek_privdoc  0.199386776
# 123        f4a_fuel_kero  0.195828314
# 124        f4a_fuel_coal  0.173718926
# 125  f4a_water_shallwell  0.155157137
# 126      f4a_water_house  0.133233067
# 127         f4a_seek_doc  0.129610448
# 128       f4a_fuel_grass  0.107568091
# 129  f4a_water_prospring  0.104023343
# 130     f4a_hometrt_zinc  0.094172170
# 131       f4a_house_boat  0.073173049
# 132     f4a_fuel_propane  0.004994769
# 133        f4a_fuel_dung  0.000000000

# AUC          SE     lower     upper level Model nvar
# 1 0.7416872 0.002664002 0.7364659 0.7469086  0.95    LR    2
# 2 0.7201292 0.002874937 0.7144945 0.7257640  0.95    LR    5
# 3 0.7284253 0.002928582 0.7226854 0.7341652  0.95    LR   10
# 4 0.7664200 0.002604577 0.7613151 0.7715248  0.95    RF    2
# 5 0.7461229 0.002705956 0.7408194 0.7514265  0.95    RF    5
# 6 0.7571733 0.002660589 0.7519586 0.7623880  0.95    RF   10

# nvar    intc intc_LCI intc_UCI slope slope_LCI slope_UCI
# <dbl>   <dbl>    <dbl>    <dbl> <dbl>     <dbl>     <dbl>
#   1     2 0.00308   -0.240    0.236 1.01      0.741      1.30
# 2     5 0.00342   -0.241    0.237 0.996     0.726      1.27
# 3    10 0.00394   -0.242    0.240 0.973     0.716      1.24

shigella0.5MUAC_Afr_2var <- glm(shigella_afe0.5~base_age+f4a_drh_blood,
                            data=Afr,family="binomial",control=glm.control(maxit=50))
summary(shigella0.5MUAC_Afr_2var)
# Estimate Std. Error z value Pr(>|z|)    
# (Intercept)    -2.312031   0.103476 -22.344  < 2e-16 ***
#   base_age        0.023508   0.004072   5.773 7.78e-09 ***
#   f4a_drh_blood1  1.803718   0.124945  14.436  < 2e-16 ***

round(exp(coef(shigella0.5MUAC_Afr_2var)),2)
# (Intercept)       base_age f4a_drh_blood1 
# 0.10           1.02           6.07 
round(exp(confint(shigella0.5MUAC_Afr_2var)),2)
# 2.5 % 97.5 %
#   (Intercept)     0.08   0.12
# base_age        1.02   1.03
# f4a_drh_blood1  4.75   7.76

Afr_2var_AfrFit<-Afr %>% select(shigella_afe0.5,base_age,f4a_drh_blood)
Afr_2var_AfrFit$Afr_pred_glm <- as.numeric(predict(shigella0.5MUAC_Afr_2var,newdata=Afr_2var_AfrFit,type="response"))
Afr_2var_AfrFit <- roc(response=Afr_2var_AfrFit$shigella_afe0.5,predictor=Afr_2var_AfrFit$Afr_pred_glm)
paste(round(Afr_2var_AfrFit$auc,2)," (",
      round(ci.auc(Afr_2var_AfrFit)[1],2),", ",
      round(ci.auc(Afr_2var_AfrFit)[3],2),")",sep="")
# "0.74 (0.71, 0.76)"

Asia_2var_AfrFit<-Asia.noBdesh %>% select(shigella_afe0.5,base_age,f4a_drh_blood)
Asia_2var_AfrFit$Afr_pred_glm <- as.numeric(predict(shigella0.5MUAC_Afr_2var,newdata=Asia_2var_AfrFit,type="response"))
Asia_2var_AfrFit <- roc(response=Asia_2var_AfrFit$shigella_afe0.5,predictor=Asia_2var_AfrFit$Afr_pred_glm)
paste(round(Asia_2var_AfrFit$auc,2)," (",
      round(ci.auc(Asia_2var_AfrFit)[1],2),", ",
      round(ci.auc(Asia_2var_AfrFit)[3],2),")",sep="")
# "0.77 (0.74, 0.8)"

####################
#shigella AFe>=0.5 only MUAC by country - Asia (not including Bdesh), val in Africa
################### var screening, AUC ####
names<-append(x=names, values=c("f4b_muac")) #add MUAC to list
names <- names[!names %in% c("f4b_haz")] #take HAZ off list

shigella0.5MUAC.Asia.noBdesh <- CPR.funct(data=Asia.noBdesh,outcome="shigella_afe0.5",iter=100,nvars_opts=c(2,5,10))
shigella0.5MUAC.Asia.noBdesh[["df_imps"]]
shigella0.5MUAC.Asia.noBdesh[["AUC_df"]]
shigella0.5MUAC.Asia.noBdesh[["calib"]]

# names      var_red
# 1          f4a_drh_blood 29.095330103
# 2               base_age 16.453124695
# 3               f4b_muac 11.630932970
# 4               f4b_resp 10.385011302
# 5               f4b_temp  9.020516363
# 6          f4a_ppl_house  7.874656744
# 7          f4a_share_fac  6.354448717
# 8           f4a_drh_days  5.311638480
# 9       f4a_yng_children  4.618553092
# 10         f4a_prim_schl  4.391482715
# 11              f4b_eyes  4.283330465
# 12         f4b_recommend  4.149788132
# 13         f4a_drh_vomit  3.886504395
# 14         f4a_slp_rooms  3.743568867
# 15       f4a_water_avail  3.575185139
# 16        f4a_disp_feces  3.439922286
# 17        f4a_trt_method  3.393119477
# 18        f4a_max_stools  3.214753499
# 19        f4a_offr_drink  3.160253216
# 20         f4a_breastfed  2.922849893
# 21            f4b_mental  2.551085248
# 22             f4b_mouth  2.456903154
# 23          f4a_ms_water  2.382639542
# 24     f4a_drh_bellypain  2.353316608
# 25         f4a_drh_cough  2.262631709
# 26             f3_gender  2.042194321
# 27       f4a_cur_thirsty  1.978340825
# 28          f4a_cur_skin  1.959714943
# 29        f4a_wash_nurse  1.945323941
# 30         f4a_wash_cook  1.924957087
# 31        f4a_house_tele  1.899726573
# 32        f4a_wash_child  1.846567094
# 33       f4a_hometrt_ors  1.832518890
# 34          f4a_wash_eat  1.809399586
# 35        f4a_drh_thirst  1.791580131
# 36      f4a_hometrt_none  1.742748324
# 37     f4a_drh_lessdrink  1.722706377
# 38       f4a_house_phone  1.721428778
# 39         f4a_trt_water  1.704452932
# 40          f4a_wash_def  1.679587705
# 41  f4a_drh_lethrgy_miss  1.676053512
# 42          f4a_wash_use  1.673317615
# 43      f4a_relationship  1.618783198
# 44     f4a_hometrt_maize  1.574840090
# 45      f4a_cur_restless  1.551631752
# 46      f4a_cur_drymouth  1.549447991
# 47           f4a_ani_cat  1.542930085
# 48      f4a_seek_outside  1.502562180
# 49      f4a_drh_restless  1.500970508
# 50      f4a_house_fridge  1.474432225
# 51          f4a_ani_fowl  1.455164092
# 52       f4a_house_radio  1.440402805
# 53              f4b_skin  1.402841105
# 54       f4a_water_house  1.393965678
# 55         f4a_fac_waste  1.383895080
# 56        f4a_drh_strain  1.377774767
# 57         f4a_fuel_wood  1.357337529
# 58        f4b_under_nutr  1.353515484
# 59           f4a_ani_dog  1.269482562
# 60        f4a_house_bike  1.264953414
# 61       f4a_fuel_natgas  1.256724797
# 62         f3_drh_turgor  1.253494764
# 63      f4a_water_pubtap  1.239893299
# 64       f4a_store_water  1.220792360
# 65       f4a_house_scoot  1.214412758
# 66            f4a_ani_no  1.212242755
# 67        f4a_water_yard  1.182701419
# 68             f4a_floor  1.179498964
# 69     f4a_hometrt_othr1  1.174429232
# 70    f4a_cur_fastbreath  1.148660525
# 71          f4a_ani_goat  1.103853735
# 72             f3_drh_iv  1.083611810
# 73         f4a_fuel_kero  1.082351403
# 74         f4a_fuel_coal  1.078000030
# 75      f4a_fuel_propane  1.029264873
# 76      f4a_water_bought  1.024489061
# 77      f4a_seek_privdoc  0.966533958
# 78       f4a_wash_animal  0.942630577
# 79           f4a_ani_cow  0.892410867
# 80       f4a_ani_rodents  0.891110620
# 81                  site  0.885623070
# 82           f3_drh_hosp  0.851702408
# 83             f4b_admit  0.811677039
# 84        f4a_hometrt_ab  0.773412086
# 85          f4a_dad_live  0.770031274
# 86         f4a_drh_consc  0.740379278
# 87          f4a_seek_doc  0.697453777
# 88        f4a_seek_pharm  0.637929380
# 89      f4a_hometrt_milk  0.572408969
# 90         f4a_wash_othr  0.571999771
# 91      f4a_hometrt_zinc  0.561839825
# 92         f4a_ani_other  0.554101416
# 93     f4a_hometrt_othr2  0.494429427
# 94        f4a_house_elec  0.470372093
# 95        f4a_water_othr  0.464010523
# 96      f4a_drh_prolapse  0.461089591
# 97          f4b_abn_hair  0.453815618
# 98    f4a_water_deepwell  0.437486767
# 99         f4a_house_car  0.414002030
# 100       f4a_house_none  0.397092908
# 101       f4a_seek_other  0.373974943
# 102       f4a_fuel_grass  0.364806031
# 103      f4b_chest_indrw  0.360302888
# 104       f4a_house_boat  0.357187969
# 105     f4a_hometrt_herb  0.274457966
# 106  f4a_hometrt_othrliq  0.251118354
# 107     f4a_house_agland  0.237603716
# 108       f4b_skin_flaky  0.225250096
# 109           f4b_rectal  0.203053159
# 110       f4a_house_cart  0.169808532
# 111       f4a_water_well  0.162458299
# 112        f4a_fuel_dung  0.160053595
# 113       f4a_water_pond  0.140588484
# 114         f4a_drh_conv  0.071871934
# 115       f4a_seek_remdy  0.063395455
# 116       f4a_water_bore  0.060791094
# 117      f4a_seek_healer  0.060019150
# 118        f4a_fuel_elec  0.055218693
# 119    f4a_fuel_charcoal  0.051462626
# 120    f4a_water_covwell  0.022706562
# 121        f4a_ani_sheep  0.022138112
# 122   f4a_water_covpwell  0.020278636
# 123  f4a_water_shallwell  0.010671819
# 124      f4a_fuel_biogas  0.003539304
# 125          f4b_bipedal  0.002922190
# 126      f4a_seek_friend  0.002812229
# 127        f4a_fuel_crop  0.000000000
# 128       f4a_fuel_other  0.000000000
# 129  f4a_water_prospring  0.000000000
# 130   f4a_water_unspring  0.000000000
# 131    f4a_water_pubwell  0.000000000
# 132      f4a_water_river  0.000000000
# 133       f4a_water_rain  0.000000000

# AUC          SE     lower     upper level Model nvar
# 1 0.7633949 0.003257704 0.7570099 0.7697798  0.95    LR    2
# 2 0.7510761 0.003325915 0.7445574 0.7575948  0.95    LR    5
# 3 0.7457770 0.003450045 0.7390150 0.7525389  0.95    LR   10
# 4 0.7625287 0.003264888 0.7561297 0.7689278  0.95    RF    2
# 5 0.7621163 0.003307360 0.7556340 0.7685986  0.95    RF    5
# 6 0.7697784 0.003211106 0.7634847 0.7760720  0.95    RF   10

# nvar     intc intc_LCI intc_UCI slope slope_LCI slope_UCI
# <dbl>    <dbl>    <dbl>    <dbl> <dbl>     <dbl>     <dbl>
#   1     2 -0.00827   -0.309    0.280 0.992     0.731      1.27
# 2     5 -0.00830   -0.309    0.280 0.977     0.718      1.25
# 3    10 -0.0115    -0.314    0.279 0.922     0.671      1.19

shigella0.5MUAC_Asia_2var <- glm(shigella_afe0.5~base_age+f4a_drh_blood,
                             data=Asia.noBdesh,family="binomial",control=glm.control(maxit=50))
summary(shigella0.5MUAC_Asia_2var)
# Estimate Std. Error z value Pr(>|z|)    
# (Intercept)    -2.183348   0.127423 -17.135  < 2e-16 ***
#   base_age        0.023378   0.004866   4.804 1.55e-06 ***
#   f4a_drh_blood1  2.436047   0.154864  15.730  < 2e-16 ***

round(exp(coef(shigella0.5MUAC_Asia_2var)),2)
# (Intercept)       base_age f4a_drh_blood1 
# 0.11           1.02          11.43 
round(exp(confint(shigella0.5MUAC_Asia_2var)),2)
# 2.5 % 97.5 %
#   (Intercept)     0.09   0.14
# base_age        1.01   1.03
# f4a_drh_blood1  8.46  15.54

Asia_2var_AsiaFit<-Asia.noBdesh %>% select(shigella_afe0.5,base_age,f4a_drh_blood)
Asia_2var_AsiaFit$Asia_pred_glm <- as.numeric(predict(shigella0.5MUAC_Asia_2var,newdata=Asia_2var_AsiaFit,type="response"))
Asia_2var_AsiaFit <- roc(response=Asia_2var_AsiaFit$shigella_afe0.5,predictor=Asia_2var_AsiaFit$Asia_pred_glm)
paste(round(Asia_2var_AsiaFit$auc,2)," (",
      round(ci.auc(Asia_2var_AsiaFit)[1],2),", ",
      round(ci.auc(Asia_2var_AsiaFit)[3],2),")",sep="")
# [1] "0.77 (0.74, 0.8)"

Afr_2var_AsiaFit<-Afr %>% select(shigella_afe0.5,base_age,f4a_drh_blood)
Afr_2var_AsiaFit$Asia_pred_glm <- as.numeric(predict(shigella0.5MUAC_Asia_2var,newdata=Afr_2var_AsiaFit,type="response"))
Afr_2var_AsiaFit <- roc(response=Afr_2var_AsiaFit$shigella_afe0.5,predictor=Afr_2var_AsiaFit$Asia_pred_glm)
paste(round(Afr_2var_AsiaFit$auc,2)," (",
      round(ci.auc(Afr_2var_AsiaFit)[1],2),", ",
      round(ci.auc(Afr_2var_AsiaFit)[3],2),")",sep="")
# [1] "0.74 (0.71, 0.76)"

