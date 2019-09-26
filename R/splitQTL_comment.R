#' Sub-function to split comment column from QTL output
#' 
#' Takes a list of candidate markers and search for genes a determined interval
#' @param output_qtls Output from QTL annotation 

splitQTL_comment<-function(output_qtls){
  #Splitting extra_info column
  output_qtls<-as.data.frame(output_qtls)
  output_qtls$QTL_type<-as.character(output_qtls$QTL_type)
  
  
  tmp.type<-sapply(strsplit(as.character(output_qtls[grep("_Association",output_qtls$QTL_type),"QTL_type"]), "_Association"), `[`, 1)
  if(!is.null(tmp.type)){
    output_qtls[grep("_Association",output_qtls$QTL_type),"Association_type"]<-"Association"
    output_qtls[grep("_Association",output_qtls$QTL_type),"QTL_type"]<-tmp.type
  }
  
  tmp.type<-sapply(strsplit(as.character(output_qtls[grep("_QTL",output_qtls$QTL_type),"QTL_type"]), "_QTL"), `[`, 1)
  if(!is.null(tmp.type)){
    output_qtls[grep("_QTL",output_qtls$QTL_type),"QTL_type"]<-"QTL"
    output_qtls[grep("_QTL",output_qtls$QTL_type),"QTL_type"]<-tmp.type
  }
  
  tmp.type<-sapply(strsplit(as.character(output_qtls[grep("_Mendelian",output_qtls$QTL_type),"QTL_type"]), "_Mendelian"), `[`, 1)
  if(!is.null(tmp.type)){
    output_qtls[grep("_Mendelian",output_qtls$QTL_type),"QTL_type"]<-"Mendelian"
    output_qtls[grep("_Mendelian",output_qtls$QTL_type),"QTL_type"]<-tmp.type
  }
  
  
  for(i in 1:nrow(output_qtls)){
    
    tmp.split<-unlist(strsplit(as.character(output_qtls[i,"extra_info"]),";"))
    
    QTL_ID<-tmp.split[grep("QTL_ID",tmp.split)]
    QTL_ID<-unlist(strsplit(as.character(QTL_ID),"="))[2]
    if(!is.null(QTL_ID)){
      output_qtls[i,"QTL_ID"]<-QTL_ID
    }
    
    trait_ID<-tmp.split[grep("trait_ID",tmp.split)]
    trait_ID<-unlist(strsplit(as.character(trait_ID),"="))[2]
    if(!is.null(trait_ID)){
      output_qtls[i,"trait_ID"]<-trait_ID
    }
    
    breed<-tmp.split[grep("breed",tmp.split)]
    breed<-unlist(strsplit(as.character(breed),"="))[2]
    if(!is.null(breed)){
      output_qtls[i,"breed"]<-breed
    }
    
    
    Name<-tmp.split[grep("Name",tmp.split)]
    Name<-unlist(strsplit(as.character(Name),"="))[2]
    if(!is.null(Name)){
      output_qtls[i,"Name"]<-Name
    }
    
    Abbrev<-tmp.split[grep("Abbrev",tmp.split)]
    Abbrev<-unlist(strsplit(as.character(Abbrev),"="))[2]
    if(!is.null(Abbrev)){
      output_qtls[i,"Abbrev"]<-Abbrev
    }
    
    Model<-tmp.split[grep("Model",tmp.split)]
    Model<-unlist(strsplit(as.character(Model),"="))[2]
    if(!is.null(Model)){
      output_qtls[i,"Model"]<-Model
    }
    
    Test_Base<-tmp.split[grep("Test_Base",tmp.split)]
    Test_Base<-unlist(strsplit(as.character(Test_Base),"="))[2]
    if(!is.null(Test_Base)){
      output_qtls[i,"Test_Base"]<-Test_Base
    }
    
    pubmed<-tmp.split[grep("PUBMED_ID",tmp.split)]
    pubmed<-unlist(strsplit(as.character(pubmed),"="))[2]
    if(!is.null(pubmed)){
      output_qtls[i,"pubmed_id"]<-pubmed
    }
    
    p_value<-tmp.split[grep("P-value",tmp.split)]
    p_value<-unlist(strsplit(as.character(p_value),"="))[2]
    if(!is.null(p_value)){
      output_qtls[i,"p_value"]<-p_value
    }
    
    bayes<-tmp.split[grep("Bayes-value",tmp.split)]
    bayes<-unlist(strsplit(as.character(bayes),"="))[2]
    if(!is.null(bayes)){
      output_qtls[i,"bayes_value"]<-bayes
    }
    
    FlankMarkers<-tmp.split[grep("FlankMarkers",tmp.split)]
    FlankMarkers<-unlist(strsplit(as.character(FlankMarkers),"="))[2]
    if(!is.null(FlankMarkers)){
      output_qtls[i,"Flank_Markers"]<-FlankMarkers
    }
  }
  
  output_qtls[,-which(colnames(output_qtls)%in%"extra_info")]
}