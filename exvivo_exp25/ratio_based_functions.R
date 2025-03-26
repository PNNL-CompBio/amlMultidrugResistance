##create signatures based on parameters


##load drug names
druglist<<-list(UT='none',D='Decitabine',V='Venetoclax',G='Gilteritinib',
              GV='Gilteritinib+Venetoclax',GD='Gilteritinib+Decitabine',DV='Decitabine+Venetoclax',
              GVD='Gilteritinib+Decitabine+Venetoclax')

treatlist<<-druglist
names(treatlist)=paste(names(treatlist),'-UT',sep='')


###load signature file
sigs <<- readr::read_tsv(syn$get('syn64943910')$path)|>
  mutate(marker=ifelse(logFC>0,'resistant','sensitive'))|>
  mutate(sig=adj.P.Val<0.05)|>
  subset(endsWith(contrast,'UT'))|>
  mutate(treatment=treatlist[contrast])



sigs|>subset(sig)|>
  group_by(marker,contrast)|>
  summarize(nMarkers=n())

##permanent colors fors ignature
sigcolors<<-c(Ratio='#666666',`resistant signal`="#9999DD",`sensitive signal`='#AAAA77',none='#EEEEEE')

#colors for true false
tfcolors<<-c(`TRUE`='#DD7333',`FALSE`='#4477EE')

##calculate signal from data
## this is the core functionality - it takes the signatures and a specific condition and returns the signal found in the data
getSigsFromData<-function(condition='V-UT',numProts=10,protList=c(),lfcthresh=0.0,byPval=TRUE,doPlot=FALSE,proteins){

  if(length(protList)==0)
    protList=unique(sigs$feature)

  ##take all the resistant indicators - features that are up-regulated upon drug treatment
  resistant_markers <- sigs %>%
    filter(contrast == condition)|>
    #subset(sig)|>
    arrange(t_test_pval) %>%
    subset(feature%in%protList)|>
    subset(logFC>lfcthresh)|>
    slice(1:10)|>
    pull(feature)|>
    unique()

  ##take all the 'sensitive' indicators, drugs that are  down-regulated upon drug treatment
  sensitive_markers <- sigs %>%
    filter(contrast == condition)|>
    arrange(t_test_pval) %>%
    #subset(sig)|>
    subset(feature%in%protList)|>
    subset(logFC<(-1*lfcthresh))|>
    slice(1:10)|>
    pull(feature)|>
    unique()

  if(doPlot){
    ###let's plot the proteins
    sigs<-sigs|>
      mutate(sigMarker=ifelse(feature%in%sensitive_markers,'sensitive signal',ifelse(feature%in%resistant_markers,'resistant signal','none')))
    lfcs<-sigs|>
      subset(contrast==condition)|>
      arrange(logFC)#|>
    #  subset(sigMarker!='none')

    p<-ggplot(lfcs,aes(y=-1*log10(t_test_pval),x=logFC,col=sigMarker))+
      geom_point()+
      scale_color_manual(values=sigcolors)+
      theme_bw()+
      ggtitle(paste(treatlist[[condition]],'Signature markers by p-value'))
    print(p)
  }

  ## TODO: add in some check about the 'quality' of the signature in the original data?
  ##return markers
  return(list(resistant=resistant_markers,sensitive=sensitive_markers))

}


##for each sample in an expression matrix, calculate the ratio
calcRatio<-function(resistant_markers,sensitive_markers,exprData){

  ##calculate the signature scores in patients amples
  sensitive_signal <- apply(exprData[sensitive_markers, ], 2, mean,na.rm=T)
  resistant_signal <- apply(exprData[resistant_markers, ], 2, mean,na.rm=T)
  ratio <- resistant_signal - sensitive_signal

  return(cbind(ratio,resistant_signal,sensitive_signal))

}
