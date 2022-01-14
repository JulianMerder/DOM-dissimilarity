JSD_chemtest<-function(D,method="reference",index="A",ref=1,permute=F,repe=1000,rename=F,samplevector=NULL,sqrtJSD=F){
  
  if(is.null(samplevector)){
    samples<-grep("ample",names(D),value=T)
  } else{
    samples<-names(D)[samplevector]  
  }
  
  namescheck <- c("C13","O18","S34","N15","Cl","Cl37","Br81","Br79","I")
  D[namescheck[!(namescheck %in% colnames(D))]] = 0
  D$totalc<-D$C+D$C13
  D <- D %>% mutate(
  AI=round(((1+totalc-O-O18-S34-S-0.5*(H+N+N15+P))/(totalc-O-O18-N15-S34-N-S-P)),2),
  AI=ifelse(AI<0 | is.na(AI) | is.infinite(AI),0,AI),
  AI.mod=ifelse(((C+C13)-0.5*(O+O18)-S-N-P-N15-S34)<0 | (1+(C+C13)-0.5*(O+O18)-S-S34 -0.5*(H+N+N15+P))<0 | ((1+(C+C13)-0.5*(O+O18)-S-S34 -0.5*(H+N+P+N15))/((C+C13)-0.5*(O+O18)-S-N-N15-S34-P))<0,0,((1+(C+C13)-0.5*(O+O18)-S-S34-0.5*(H+N+P+N15))/((C+C13)-0.5*(O+O18)-S-S34-N-N15-P))),
  AI.mod=ifelse(AI.mod<0 | is.na(AI.mod) | is.infinite(AI.mod),0,AI.mod),
  DBE=round(1+0.5*(2*totalc-H-Cl-Cl37-Br81-Br79-I+N+P+N15),2))   
  D<-as.data.frame(D)
  D$Aromatic=ifelse(D$AI.mod>0.5,1,0)
  D$Aromatic.O_rich=ifelse(D$AI.mod>0.5 & D$O.C>0.5,1,0)
  D$Aromatic.O_poor=ifelse(D$AI.mod>0.5 & D$O.C<=0.5,1,0)
  D$Highly.unsaturated=ifelse(D$AI.mod<=0.5 & D$H.C<1.5,1,0)
  D$Highly.unsaturated.O_rich=ifelse(D$AI.mod<=0.5 & D$H.C<1.5& D$O.C>0.5,1,0)
  D$Highly.unsaturated.O_poor=ifelse(D$AI.mod<=0.5 & D$H.C<1.5& D$O.C<=0.5,1,0)
  D$Unsaturated=ifelse(D$H.C>=1.5 & D$DBE!=0 & D$H.C<=2,1,0)
  D$Unsaturated.O_rich=ifelse(D$H.C>=1.5 & D$DBE!=0 & D$O.C>0.5 & D$H.C<=2,1,0)
  D$Unsaturated.O_poor=ifelse(D$H.C>=1.5 & D$DBE!=0 & D$O.C<=0.5 & D$H.C<=2,1,0)
  D$Unsaturated.with.N=ifelse(D$H.C>=1.5 & D$DBE!=0 & (D$N+D$N15)>0 & D$H.C<=2,1,0)
  D$Saturated=ifelse(D$DBE==0,1,0)
  D$Saturated.O_rich=ifelse(D$DBE==0& D$O.C>0.5,1,0)
  D$Saturated.O_poor=ifelse(D$DBE==0& D$O.C<=0.5,1,0)
  D$AI.mod<-round(D$AI.mod,2)
  D$rest<-rowSums(D[,c("Cl","Cl37","Br81","Br79","I")])
  D$CH<-ifelse((D$C+D$C13)>0&D$H>0&(D$O+D$O18)==0&(D$S+D$S34)==0&(D$N+D$N15)==0&D$rest==0,1,0)
  D$CHO<-ifelse((D$C+D$C13)>0&D$H>0&(D$O+D$O18)>0&(D$S+D$S34)==0&(D$N+D$N15)==0&D$rest==0,1,0)
  D$CHNO<-ifelse((D$C+D$C13)>0&D$H>0&(D$O+D$O18)>0&(D$S+D$S34)==0&(D$N+D$N15)>0&D$rest==0,1,0)
  D$CHOS<-ifelse((D$C+D$C13)>0&D$H>0&(D$O+D$O18)>0&(D$S+D$S34)>0&(D$N+D$N15)==0&D$rest==0,1,0)
  D$CHNOS<-ifelse((D$C+D$C13)>0&D$H>0&(D$O+D$O18)>0&(D$S+D$S34)>0&(D$N+D$N15)>0&D$rest==0,1,0)
  D$newcategory<-1+(D$C+D$C13)-0.5*(D$O+D$O18)-(D$S+D$S34)-0.5*(D$H+D$N+D$P+D$N15)
  D$newcategory<-round(D$newcategory/0.5)*0.5
  
  if(index=="A") categories<-c("Aromatic","Highly.unsaturated","Unsaturated","Saturated") #need to be unique for each formula
  if(index=="A1") categories<-c("Aromatic.O_rich","Aromatic.O_poor","Highly.unsaturated.O_rich","Highly.unsaturated.O_poor","Unsaturated.O_rich","Unsaturated.O_poor","Saturated.O_rich","Saturated.O_poor") #need to be unique for each formula
  if(index=="B") categories<-c("CH","CHO","CHNO","CHOS","CHNOS")
  if(index=="C") categories<-c("newcategory")
  
  #if(isos)  D$newcategory<-1+(D$C+D$C13)-0.5*(D$O+D$O18)-(D$S+D$S34)-0.5*(D$H+D$N+D$P+D$N15)
  #if(!isos) D$newcategory<-1+(D$C)-0.5*(D$O)-(D$S)-0.5*(D$H+D$N+D$P)
  
  K<-D[,c(categories,samples)]
  
  K1<-K %>%
    pivot_longer(cols=all_of(samples), names_to = "sample", values_to = "intensity")
  K1$intensity[K1$intensity==0]<-NA
  K1<-K1[!is.na(K1$intensity),]
  
  
  #K2<-K1 %>% group_by(sample) %>% summarize("Aromatic"=sum(Aromatic==1),"Highly.unsaturated"=sum(Highly.unsaturated==1),"Unsaturated"=sum(Unsaturated==1),"Saturated"=sum(Saturated==1))%>%
    #arrange(factor(sample,levels=samples))
  #K3<-as.matrix(K2[,2:5])
  #rownames(K3)<-K2$sample
  if(index %in% c("A","A1","B")){
  K2<-K1 %>% group_by(sample) %>% summarise_at(categories, sum, na.rm = TRUE)%>%
    arrange(factor(sample,levels=samples))
  K3<-as.matrix(K2[,2:5])
  rownames(K3)<-K2$sample
  #K4<-KL(K3, est.prob = "empirical",unit="log2")
  K4<-JSD(K3, est.prob = "empirical",unit="log2")
  K4m<-matrix(NA,nrow(K3),ncol=repe)
  
  }
  if(index %in% c("C")){
  K2<-K1 %>% group_by(sample,newcategory) %>% summarize(value=n())
  K6 <- K2 %>%
    pivot_wider(names_from = newcategory, values_from = value)
  p1<-sort(as.numeric(colnames(K6)[-1]))
  
  K3<- as.data.frame(K6[,as.character(c("sample",p1))])
  rownames(K3)<-K6$sample
  K3[is.na(K3)]<-0
  K3 <- K3 %>% arrange(factor(sample,levels=samples))
  rownames(K3)<-K3$sample
  K3<-K3[,-1]
  #K8<-suppressMessages(KL(as.matrix(K7), est.prob = "empirical",unit="log2"))
  K4<-suppressMessages(JSD(as.matrix(K3), est.prob = "empirical",unit="log2"))
  K4m<-matrix(NA,nrow(K3),ncol=repe)
  }
  
  if(method=="matrix"){
    ok4 <- upper.tri(K4,diag = F)
    ok4<-cbind(row = row(K4)[ok4], col = col(K4)[ok4])
    K4m<-matrix(NA,nrow(ok4),ncol=repe)
  }
  
  
  #permutations
  
  
  if(isTRUE(permute)){
    for (i in 1:repe){
      if(method=="reference"){
        weights7<-unlist(K3[ref,])+1
        if(index %in% c("C")){
        P7<-t(apply(K3,1,FUN=function(x){table(factor(sample(p1, size=sum(x), replace=TRUE,prob = weights7),levels=p1))}))
        }else{
        P7<-t(apply(K3,1,FUN=function(x){table(factor(sample(categories, size=sum(x), replace=TRUE,prob = weights7),levels=categories))}))
        }
        P7[ref,]<-unlist(K3[ref,])
        K4m[,i]<-suppressMessages( JSD(as.matrix(P7), est.prob = "empirical",unit="log2")[ref,] ) 
        
      }
      if(method=="matrix"){
        ok <- ok4
        for(j in 1:nrow(ok)){
          row=ok[j,1]
          col=ok[j,2]
          if(index %in% c("C")){
          P7<-data.frame(A=unlist(K3[row,]),B=as.numeric(unlist(table(factor(sample(p1, size=sum(K3[col,]), replace=TRUE,prob = unlist(K3[row,])+1),levels=p1)))))
          }else{
          P7<-data.frame(A=unlist(K3[row,]),B=as.numeric(unlist(table(factor(sample(categories, size=sum(K3[col,]), replace=TRUE,prob = unlist(K3[row,])+1),levels=categories)))))
          }
          #F7<-suppressMessages( KL(t(P7), est.prob = "empirical",unit="log2"))
          F7<-suppressMessages( JSD(t(P7), est.prob = "empirical",unit="log2"))
          K4m[j,i]<-F7
        }  
        
        
      }
      
    }
  } 
  
  if(method=="reference"){
    p8<-rowSums(K4[ref,]<=K4m)/repe
    if(isTRUE(sqrtJSD)) p8<-rowSums(sqrt(K4[ref,])<=sqrt(K4m))/repe
  }
  
  
  
  
  
  if(method=="matrix"){
    p8<-matrix(NA,nrow=nrow(K4),ncol=ncol(K4))
    ok <- upper.tri(K4,diag = F)
    ok<-cbind(row = row(K4)[ok], col = col(K4)[ok])
    p8x<-sapply(1:nrow(K4m),FUN=function(x){
      row=ok[x,1]
      col=ok[x,2]
      aw<-sum(K4[row,col]<=K4m[x,],na.rm=T)
      if(isTRUE(sqrtJSD)) aw<-sum(sqrt(K4[row,col])<=sqrt(K4m[x,]),na.rm=T)
      aw
      })/repe
      
    p8[ok]<-p8x
    diag(p8)<-1
    p8[lower.tri(p8)] <- t(p8)[lower.tri(p8)]
    
  }
  
  
  if(method!="matrix"){
    p8[p8==0]<-0.001
  }
  
  
  
  
  
  if(method=="reference"){
    
    smple<-rownames(K3)
    if(isTRUE(rename)){
      smple<-sub("\\_.*", "", smple)
      
    }
    if(isTRUE(sqrtJSD)) K4<-sqrt(K4)
    final<-data.frame(sample=smple,JS_divergence_to_reference=c(K4[ref,]),case=c(rep(index,nrow(K4))),pvalue=c(p8))
    final$sample<-factor(final$sample,levels=smple)
    plot1a<-ggplot(final,aes(x=sample,y=JS_divergence_to_reference,color=case))+geom_point()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ylab("JS to reference")
    
    colorr<-as.character(final$pvalue<=0.01)
    colorr[final$sample==samples[ref]]<-"reference"
    colorr<-factor(colorr,levels=c("FALSE","TRUE","reference"))
    plot1b<-ggplot(final,aes(x=sample,y=JS_divergence_to_reference,color=colorr))+geom_point()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ylab("JSD to reference")+
      labs(color = "p <= 0.01")
    
  }
  
  if(method=="matrix"){
    row.names(K4)<-rownames(K3)
    colnames(K4)<-rownames(K3)
    if(isTRUE(rename)){
      row.names(K4)<-sub("\\_.*", "", rownames(K3))
      colnames(K4)<-sub("\\_.*", "", rownames(K3))
    }
    
    if(isTRUE(sqrtJSD)) K4<-sqrt(K4)
    final<-list(JSresults=K4,p_values=p8)  
    if(!isTRUE(permute)){
      final<-list(JSresults=K4,p_values="not selected")  
      
    }
    plot1a<-ggplot(reshape2::melt(t(K4)), aes(Var1, Var2, fill=value)) +
      geom_tile(height=1, width=1) +
      scale_fill_viridis_c() +
      theme_minimal() +
      coord_equal() +
      labs(x="",y="",fill="JSd") +
      theme(axis.text.x=element_text(size=7, margin=margin(0,-3,0,0),angle=90),
            axis.text.y=element_text(size=7, margin=margin(0,-3,0,0)),
            panel.grid.major=element_blank())+ggtitle(index)
    
      
    
    K4[p8>0.05]<-NA
    
    
    
    plot1b<-ggplot(reshape2::melt(t(K4)), aes(Var1, Var2, fill=value,color="")) +
      geom_tile(height=1, width=1) +
      scale_fill_viridis_c(na.value="white") +scale_colour_manual(values=NA) +              
      guides(colour=guide_legend("p > 0.05", override.aes=list(fill="white",color="black"),order=2))+
      theme_minimal() +
      coord_equal() +
      labs(x="",y="",fill="JSd") +
      theme(axis.text.x=element_text(size=7, margin=margin(0,-3,0,0),angle=90),
            axis.text.y=element_text(size=7, margin=margin(0,-3,0,0)),
            panel.grid.major=element_blank())+ggtitle(index)
    
    
    
    
    
  }
  
  if(!isTRUE(permute)){
    plot1b=ggplot() +
      theme_void() +
      geom_text(aes(0,0,label="Note: no p-values calculated (set permute to TRUE in function call)")) +
      xlab(NULL)
    
  }
  return(list(plot1=plot1a,plot2=plot1b,JSresults=final))
  
  
}



KU<-read.csv("Crosstab.csv")

library(dplyr)
library(philentropy)
library(tidyr)
library(ggplot2)
A<-JSD_chemtest(D=KU,method="matrix",index="C",permute=T,rename=TRUE)
A$plot1
A$plot2
A1<-JSD_chemtest1(D=KU,method="matrix",permute=F,rename=TRUE)
A1$plot1
A1$plot2