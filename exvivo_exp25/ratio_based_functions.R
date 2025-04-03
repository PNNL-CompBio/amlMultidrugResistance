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

dcolors=c(none='#CCCAAA',Gilteritinib='#EEDD55',
          `Venetoclax`='#0072B2',
          `Decitabine`='#D14610',
          `Gilteritinib+Venetoclax`='#05820F',
          `Decitabine+Venetoclax`='#D172B2',
          `Gilteritinib+Decitabine`='#EDAE49',
          `Gilteritinib+Decitabine+Venetoclax`='#666666')

##permanent colors fors ignature
sigcolors<<-c(Ratio='#666666',`resistant signal`="#9999DD",`sensitive signal`='#AAAA77',none='#DDDDDD')

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
    subset(logFC>lfcthresh)
  resistant_markers = resistant_markers$feature[1:numProts]


  ##take all the 'sensitive' indicators, drugs that are  down-regulated upon drug treatment
  sensitive_markers <- sigs %>%
    filter(contrast == condition)|>
    arrange(t_test_pval) %>%
    #subset(sig)|>
    subset(feature%in%protList)|>
    subset(logFC<(-1*lfcthresh))

  sensitive_markers = sensitive_markers$feature[1:numProts]

  if(doPlot){
    ###let's plot the proteins
    sigs<-sigs|>
      mutate(sigMarker=ifelse(feature%in%sensitive_markers,'sensitive signal',ifelse(feature%in%resistant_markers,'resistant signal','none')))
    lfcs<-sigs|>
      subset(contrast==condition)|>
      arrange(logFC)#|>
    #  subset(sigMarker!='none')

 #   lfcs$significance<-apply(lfcs,1,function(x){if(x$logFC<0) return(log10(x$t_test_pval)) else return(-1*log10(x$t_test_pval))})
    lfcs<-lfcs|>
      mutate(significance=ifelse(logFC<0, log10(t_test_pval), -1*log10(t_test_pval)))|>
      arrange(significance)|>
      subset(!is.na(significance))
    lfcs$feature=factor(lfcs$feature,levels=lfcs$feature)


    p<-ggplot(lfcs,aes(y=-1*log10(t_test_pval),x=logFC,col=sigMarker,size=10))+
      geom_point()+
      scale_color_manual(values=sigcolors)+
      theme_bw()+
      ggtitle(paste("Selection of ",treatlist[[condition]],'markers by significance'))

    p2<-ggplot(lfcs[c(1:200,(nrow(lfcs)-200):nrow(lfcs)),],aes(x=feature,y=significance,fill=sigMarker))+geom_bar(stat='identity')+scale_fill_manual(values=sigcolors)+theme_bw()+coord_flip()
    print(cowplot::plot_grid(p,p2,ncol=2))

    ##also print original expression?
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
