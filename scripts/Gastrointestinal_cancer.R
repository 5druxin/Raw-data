if (T) {
  dir.create("scripts")
  dir.create("results")
  dir.create("files")
  dir.create("PDFs")
  dir.create("origin_datas/GEO",recursive = T)
  dir.create("origin_datas/TCGA")
}
library(reshape2)
library(ggpubr)
library(ggsci)
library(maftools)
library(tidyr)
library(pheatmap)
library(clusterProfiler)
library(ComplexHeatmap)
library(GSVA)
library(GSEABase)
library(corrplot)
library(dplyr)
library(survminer)
library(colorspace)
options(stringsAsFactors = F)
source('mg_base.R')
custom_theme <- function() {
  theme_survminer() %+replace%
    theme(
      plot.title=element_text(hjust=0.5)
    )
}
my_boxplot=function(dat,group,group_cols=ggsci::pal_aaas()(10),test_method='kruskal.test',
                    fill= "Group",label=c("p.format",'p.signif')[1],
                    xlab='',ylab='',title='',x.size=10,y.size=10,legend.position='top'){
  
  # dat=ARGs.score[tcga.subtype$Samples,'score']
  # group=tcga.subtype$Cluste
  data=data.frame(Group=group,value=dat)
  data=crbind2DataFrame(data)
  data=melt(data)
  data=data[which(!is.na(data[,1])),]
  p=ggplot(data, aes(x=Group, y=value,fill=Group)) +
    geom_boxplot()+
    scale_fill_manual(values = group_cols)+   #箱线图填充颜色
    ggpubr::stat_compare_means(aes(group=Group), label = label, method = test_method)+
    labs(x=xlab, y = ylab, fill = fill) +
    theme_bw()+
    theme(legend.position =legend.position,                 #图例位置
          plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(size = x.size),
          axis.text.y = element_text(size = y.size)) # 将图表标题居中
  return(p)
}
my_violin=function(dat,group,group_cols=ggsci::pal_aaas()(10),test_method='kruskal.test',
                   fill= "Group",label=c("p.format",'p.signif')[1],
                   xlab='',ylab='',title='',x.size=10,y.size=10,legend.position='top'){
  
  # dat = tcga.b.cell$GeneSet,
  # group = tcga.subtype$Cluster
  data=data.frame(Group=group,value=dat)
  data=crbind2DataFrame(data)
  data=melt(data)
  data=data[which(!is.na(data[,1])),]
  p=ggplot(data,aes(x=Group, y=value,fill=Group)) +
    geom_violin()+  
    scale_fill_manual(values = group_cols)+
    geom_boxplot(width=0.2,position=position_dodge(0.8),outlier.colour = NA,fill="white")+
    theme_classic(base_size = 20)+labs(x=xlab,y=ylab,title = title)+
    ggpubr::stat_compare_means(aes(group=Group), label = label, method =test_method)+
    theme(legend.position = legend.position,axis.text = element_text(color = 'black'),
          title = element_text(size = 12),
          axis.text.x = element_text(size = 12),axis.title.x = element_text(size = 15),
          axis.text.y = element_text(size = 12),axis.title.y = element_text(size = 12))
  return(p)
}

my_mutiboxplot=function(dat,group,group_cols=ggsci::pal_aaas()(10),
                        test_method=c('t.test','wilcox.test','anova','kruskal.test')[4],bw=T,
                        xlab='',ylab='score',title='',size=10,angle = 45, hjust = 1,legend.position='top',fill='group'){
  # dat=tcga.est[tcga.subtype.cli$Samples,]
  # group=tcga.subtype.cli$Cluster
  dat.bind=cbind(dat,Cluster=group)
  dat.bind=crbind2DataFrame(dat.bind)
  dat.melt=melt(dat.bind)
  #data=data[which(!is.na(data[,1])),]
  colnames(dat.melt)=c('Group','type','value')
  p=dat.melt %>%
    ggplot(aes(x=type, y=value,fill=Group)) +
    geom_boxplot()+
    scale_fill_manual(values =group_cols)+   #箱线图填充颜色
    ggpubr::stat_compare_means(aes(group=Group), label = "p.signif", method = test_method)+
    labs(x="", y = ylab, fill =fill,title =title) +
    #theme_light()+
    theme_bw()+
    #theme_classic()
    theme(legend.position = legend.position,                 #图例位置
          plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(size = size,angle = angle, hjust = hjust)) # 将图表标题居中
  return(p)
}

my_mutiboxplot_seg=function(dat,group,group_cols=ggsci::pal_aaas()(10),
                            test_method=c('t.test','wilcox.test','anova','kruskal.test')[4],
                            fill='group',xlab='',ylab='score',title='',xsize=10,xangle = 45,xhjust = 1,ysize=10,
                            legend.position='top',nrow,ncol){
  # dat=tcga.est[tcga.subtype.cli$Samples,]
  # group=tcga.subtype.cli$Cluster
  dat.bind=cbind(dat,Cluster=group)
  dat.bind=crbind2DataFrame(dat.bind)
  dat.melt=melt(dat.bind)
  #data=data[which(!is.na(data[,1])),]
  colnames(dat.melt)=c('Group','type','value')
  p=dat.melt %>%
    ggplot(aes(x=Group, y=value,fill=Group)) +
    geom_boxplot()+facet_wrap(~type,scales = 'free',nrow = nrow,ncol = ncol)+
    scale_fill_manual(values =group_cols)+   #箱线图填充颜色
    ggpubr::stat_compare_means(aes(group=Group), label = "p.signif", method = test_method)+
    labs(x="", y = ylab, fill = fill,title =title) +
    #theme_light()+
    #theme_classic()+
    theme(legend.position = legend.position,                 #图例位置
          plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(size = xsize,angle = xangle, hjust = xhjust),
          axis.text.y = element_text(size= ysize)) # 将图表标题居中
  return(p)
}


mg_violin_1=function(data,xangle=0,ylab='value',xlab='',leg.title='Group',test_method='anova',cmp_test_method='t.test',
                     legend.pos='r',melt=F,jitter=T,ylim=NULL,show_compare=NULL,point_size=NULL,group_col){
  library(ggplot2)
  if(is.null(ylim)){
    
  }
  if(melt){
    data_m=data
    colnames(data_m)=c('Group','value')
  }else{
    data_m=reshape2::melt(data)
    colnames(data_m)=c('Group','value')
  }
  if(!is.null(ylim)){
    data_m$value[data_m$value>ylim[2]]<-NA
  }
  data_m=data_m[which(!is.na(data_m[,1])),]
  if(xangle==0){
    tx=element_text(colour="black",family="Times")
  }else{
    tx=element_text(angle=xangle,hjust = 1,colour="black",family="Times")
  }
  
  pos='right'
  if(is.null(legend.pos)){
    pos='none'
  }else if(legend.pos=='tr'){
    pos=c(1,1)
  }else if(legend.pos=='br'){
    pos=c(1,0)
  }else if(legend.pos=='tl'){
    pos=c(0,1)
  }else if(legend.pos=='bl'){
    pos=c(0,0)
  }else if(legend.pos=='t'){
    pos='top'
  }else if(legend.pos=='r'){
    pos='right'
  }else if(legend.pos=='b'){
    pos='bottom'
  }
  uni.group=unique(data_m[,1])
  ct=length(uni.group)
  
  p1<-ggplot(data_m,aes(x=Group,y=value))+geom_violin(alpha=0.7)
  # if(ct<=4){
  #   p1=p1+ggsci::scale_fill_lancet()
  # }else if(ct<=10){
  #   p1=p1+ggsci::scale_fill_npg(name=leg.title)
  # }else if(ct<=20){
  #   p1=p1+ggsci::scale_fill_d3(palette = "category20",name=leg.title)
  # }else if(ct<=30){
  #   cbPalette=c(ggsci::pal_npg("nrc", alpha = 0.6)(10),ggsci::pal_d3("category20", alpha = 0.6)(20))
  #   p1=p1+scale_fill_manual(values=cbPalette[1:ct])
  # }else if(ct<=38){
  #   cbPalette=c(ggsci::pal_lancet()(10)
  #               ,ggsci::pal_npg("nrc", alpha = 0.6)(10)
  #               ,ggsci::pal_d3("category20", alpha = 0.6)(20)
  #               ,ggsci::pal_nejm("default", alpha = 0.6)(8))
  #   p1=p1+scale_fill_manual(values=cbPalette[1:ct])
  # }
  p1=p1+scale_fill_manual(values=group_col)
  if(jitter){
    if(is.null(point_size)){
      p1<-p1+geom_jitter(alpha=0.3,col='black',show.legend=FALSE,width = 0.2)
    }else{
      p1<-p1+geom_jitter(alpha=0.3,col='black',show.legend=FALSE,width = 0.2,size=point_size)
    }
  }
  
  p1=p1+theme_bw()+geom_boxplot(width=0.2,aes(fill=Group),outlier.shape = NA)
  p1=p1+theme(axis.text.x=tx, 
              axis.text.y=element_text(family="Times",face="plain"), 
              axis.title.y=element_text(family="Times",face="plain"), 
              legend.text=element_text(face="plain", family="Times", colour="black" 
              ),
              legend.title=element_text(face="plain", family="Times", colour="black"
              ),
              legend.justification=pos, legend.position=pos
              ,legend.background = element_rect(fill = NA, colour = NA)
  )+ylab(ylab)+xlab(xlab)
  til=''
  if(test_method=='anova'){
    if(length(unique(data_m[,1]))<3){
      x1=data_m[,2][which(data_m[,1]==unique(data_m[,1])[1])]
      x2=data_m[,2][which(data_m[,1]==unique(data_m[,1])[2])]
      pv=t.test(x1,x2)$p.value 
      til=paste0('t-tests p=',signif(pv,2))
    }else{
      fit <- aov(value~Group, data = data_m)
      pv=summary(fit)[[1]][5][[1]]
      fv=summary(fit)[[1]][4][[1]]
      til=paste0('ANOVA tests p=',signif(pv,2))
    }
  }else{
    if(length(unique(data_m[,1]))<3){
      x1=data_m[,2][which(data_m[,1]==unique(data_m[,1])[1])]
      x2=data_m[,2][which(data_m[,1]==unique(data_m[,1])[2])]
      pv=wilcox.test(x1,x2)$p.value 
      til=paste0('wilcox.tests p=',signif(pv,2))
    }else{
      fit=kruskal.test(value~Group, data = data_m)
      pv=fit$p.value
      til=paste0('Kruskal-Wallis test p=',signif(pv,2))
    }
  }
  p1=p1+ggtitle(til) 
  if(!is.null(ylim)){
    p1=p1+ylim(ylim)
  }
  if(is.null(show_compare)){
    if(length(uni.group)>5){
      show_compare=F
    }else{
      show_compare=T
    }
  }
  if(show_compare){
    comps=list()
    for(i in 1:(length(uni.group)-1)){
      for(j in (i+1):length(uni.group)){
        comps=c(comps,list(c(uni.group[i],uni.group[j])))
      }
    }
    p1=p1+ggpubr::stat_compare_means(comparisons = comps,method = cmp_test_method,label= "p.signif")
  }
  return(p1)
}
coxFun <- function(dat){
  library(survival)
  colnames(dat)=c('time','status','gene')
  fmla=as.formula("Surv(time,status)~gene")
  cox=coxph(fmla,data=dat)
  p=summary(cox)[[7]][5]
  result=c(p,summary(cox)[[8]][1],summary(cox)[[8]][3],summary(cox)[[8]][4])
  return(result)
}
get_riskscore.lasso<-function(dat,os,os.time,labels=c('A','B')){
  library(glmnet)
  set.seed(2021)
  fit1=glmnet(as.matrix(dat)
              ,cbind(time=os.time,
                     status=os)
              ,family="cox"
              ,nlambda=100
              , alpha=1) 
  
  cv.fit<-cv.glmnet(as.matrix(dat)
                    ,cbind(time=os.time,
                           status=os)
                    ,family="cox"
                    ,nfolds = 10
                    ,nlambda=100
                    , alpha=1)
  sig.coef <- coefficients(cv.fit,s=cv.fit$lambda.min)[which(coefficients(cv.fit,s=cv.fit$lambda.min)[,1]!=0),1]
  #print(cv.fit$lambda.min)
  #length(names(sig.coef))
  #10
  mg_plot_lasso <- function(fit,cv_fit,lambda=NULL,show_text=T,figLabels=c('A','B')){
    if(is.null(lambda)){
      lmda=cv_fit$lambda.min
    }else{
      lmda=lambda
    }
    fit.coef=fit$beta[(apply(fit$beta,1,function(x){
      return(sum(x!=0))
    })>0),]
    
    fit.coef=as.matrix(fit.coef)
    colnames(fit.coef)=fit$lambda
    #fit$lambda==cv_fit$lambda
    library(ggplot2)
    dat=data.table::melt(t(as.matrix(fit.coef)))
    dat_z=dat[which(dat$value==0),]
    dat=dat[which(dat$value!=0),]
    dat.sv=rbind()
    for (u in unique(dat_z[,2])) {
      t.z=dat_z[which(dat_z[,2]==u),1]
      t.zx=max(t.z)
      dat.sv=rbind(dat.sv,c(t.zx,u,0))
      t.zn=min(t.z)
      if(t.zx!=t.zn){
        dat.sv=rbind(dat.sv,c(t.zn,u,0))
      }
    }
    colnames(dat.sv)=colnames(dat_z)
    #dat_z=dat_z[dat_z[,2]%in%names(which(fit.coef[,which(fit$lambda==lmda)]!=0)),]
    dat=crbind2DataFrame(rbind(dat,dat.sv))
    mn=min(-log(dat$Var1))
    mx=max(-log(dat$Var1))
    if(show_text){
      mx=(mx-mn)*0.1+mx
    }
    p=ggplot(dat, aes(x=-log(Var1), y=value,colour=Var2))+geom_line()+theme_bw()+theme(legend.position = "none")
    p=p+coord_cartesian(xlim=c(mn, mx))+xlab('-ln(lambda)')+ylab('Coefficients')
    if(show_text){
      fl=fit.coef[which(fit.coef[,which(fit$lambda==lmda)]!=0),ncol(fit.coef)]
      for_label=data.frame(Var1=rep(min(dat$Var1),length(fl)),Var2=names(fl),value=fl)
      p=p+ggrepel::geom_label_repel(
        aes(label = Var2,color=Var2),
        data = for_label,hjust = 0
      )
    }
    p=p+geom_vline(aes(xintercept=-log(lmda)), colour="#BB0000", linetype="dashed")
    p=p+annotate('text',x=-log(lmda),y=min(dat[,3]),label=paste0('lambda=',round(lmda,4)))
    tgc=data.frame(lambda=cv_fit$lambda,cvm=cv_fit$cvm,cvup=cv_fit$cvup,cvlo=cv_fit$cvlo,cvsd=cv_fit$cvsd
                   ,col=ifelse(cv_fit$lambda>=cv_fit$lambda.min&cv_fit$lambda<=cv_fit$lambda.1se,ifelse(cv_fit$lambda==lmda,'A','C'),'B'))
    p1=ggplot(tgc, aes(x=log(lambda), y=cvm)) + xlab('ln(lambda)')+ ylab('Parial Likelihood Deviance')+
      geom_errorbar(aes(ymin=cvm-cvsd, ymax=cvm+cvsd)) +
      geom_point(aes(colour=col))
    p1=p1+theme_bw()+theme(legend.position = "none")
    gal=ggpubr::ggarrange(p,p1, ncol = 2, nrow = 1
                          #,align = "hv"
                          ,labels = figLabels)
    return(gal)
  }
  lasso.pdf <- mg_plot_lasso(fit1,
                             cv.fit,
                             show_text=T,
                             figLabels=labels)
  return(list(lasso.gene=names(sig.coef),lambda.min=cv.fit$lambda.min,plot=lasso.pdf))
}
get_riskscore<-function(dat,os,os.time,step=T,direction=c("both", "backward", "forward")[1]){
  tcga_dat1 <- cbind(time=os.time,
                     status=os,
                     dat)
  tcga_dat1=crbind2DataFrame(tcga_dat1)
  colnames(tcga_dat1)=gsub('-','__',colnames(tcga_dat1))
  gene111=gsub('-','__',colnames(dat))
  fmla <- as.formula(paste0("Surv(time, status) ~"
                            ,paste0(gene111,collapse = '+')))
  cox <- coxph(fmla, data =as.data.frame(tcga_dat1))
  if(step==T){
    cox1 <- step(cox,direction =direction)
  }else{
    cox1=cox
  }
  lan <- coef(cox1)
  #round(lan, 3)
  genes <- names(cox1$coefficients)
  mult_results=paste0(round(lan, 3), '*', names(lan),collapse = '+')
  risk.tcga <- as.numeric(lan%*%as.matrix(t(tcga_dat1[,genes])))
  
  data_gene_score_final<-tcga_dat1
  data_gene_score_final$Samples<-rownames(data_gene_score_final)
  data_gene_score_final$riskscore=risk.tcga
  data_gene_score_final$riskscorez=mosaic::zscore(risk.tcga)
  #最佳截断
  optimalCutoff <- survminer::surv_cutpoint(data.frame(time=data_gene_score_final$time/365,
                                                       event=data_gene_score_final$status,
                                                       risk=data_gene_score_final$riskscore), 
                                            time = "time", event = "event",variables = c("risk"))
  optimalCutoff=optimalCutoff$cutpoint$cutpoint[1]
  #print(optimalCutoff)
  #optimalCutoff=median(data_gene_score_final$riskscore)
  #optimalCutoff=0
  data_gene_score_final$Risk=ifelse(data_gene_score_final$riskscore>optimalCutoff,'High','Low')
  table(data_gene_score_final$Risk)
  data_gene_score_final$cutoff=optimalCutoff
  return(list(result=data_gene_score_final,module.gene=genes,model=mult_results))
}
bioForest=function(rt=null,col){
  #读取输入文件
  gene <- rownames(rt)
  hr <- sprintf("%.3f",rt[,2])
  hrLow  <- sprintf("%.3f",rt[,3])
  hrHigh <- sprintf("%.3f",rt[,4])
  Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")
  pVal <- ifelse(rt[,1]<0.001, "<0.001", sprintf("%.3f", rt[,1]))
  
  #输出图形
  n <- nrow(rt)
  nRow <- n+1
  ylim <- c(1,nRow)
  layout(matrix(c(1,2),nc=2),width=c(3,2.5))
  
  #绘制森林图左边的临床信息
  xlim = c(0,3)
  par(mar=c(4,2.5,2,1))
  plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
  text.cex=0.8
  text(0,n:1,gene,adj=0,cex=text.cex)
  text(1.5-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5-0.5*0.2,n+1,'pvalue',cex=text.cex,font=2,adj=1)
  text(3.1,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3.1,n+1,'Hazard ratio',cex=text.cex,font=2,adj=1)
  
  #绘制右边森林图
  par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
  xlim = c(0,max(as.numeric(hrLow),as.numeric(hrHigh)))
  plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio")
  arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="darkblue",lwd=2.5)
  abline(v=1,col="black",lty=2,lwd=2)
  boxcolor = ifelse(as.numeric(hr) > 1, col[2], col[1])
  points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=1.5)
  axis(1)
}
getGeneFC=function(gene.exp,group,ulab=ulab,dlab=dlab){
  degs_C1_C3=mg_limma_DEG(gene.exp, 
                          group,
                          ulab=ulab,
                          dlab=dlab)
  
  ## 未过滤任何基因
  degs_C1_C3_sig<-degs_C1_C3$DEG[which(degs_C1_C3$DEG$adj.P.Val <= 1),]
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(stringr)
  degs_C1_C3_sig_gene<-rownames(degs_C1_C3_sig)
  
  degs_gene_entz=bitr(degs_C1_C3_sig_gene,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db")
  degs_gene_entz <- dplyr::distinct(degs_gene_entz,SYMBOL,.keep_all=TRUE)
  
  gene_df <- data.frame(logFC=degs_C1_C3_sig$logFC,
                        SYMBOL = rownames(degs_C1_C3_sig))
  gene_df <- merge(gene_df,degs_gene_entz,by="SYMBOL")
  head(gene_df)
  
  geneList<-gene_df$logFC
  names(geneList)=gene_df$ENTREZID 
  head(geneList)
  
  geneList=sort(geneList,decreasing = T)
  head(geneList) ## 降序排列
  return(geneList)
}
mg_nomogram=function(clinical_riskscore,
                     os,
                     status,
                     title='Nomogram',
                     quick=T,
                     mks = c(1,3,5)){
  #clinical_riskscore=dat1[,3:5]
  #os=dat1[,1]
  #status=dat1[,2]
  #sum(is.na(norm.stat.al[,3]))
  norm.stat.al=data.frame(clinical_riskscore,time=os,status=status)
  norm.stat.al=as.data.frame(norm.stat.al)
  library(rms)
  env <- globalenv()
  env$MG_Grobal_DDSet <- rms::datadist(norm.stat.al) 
  options(datadist='MG_Grobal_DDSet')
  fmla <- as.formula(paste0("Surv(time, status) ~",paste0(colnames(clinical_riskscore),collapse = '+')))
  cox2 <- cph(fmla, data = norm.stat.al,surv = T,x = T,y = T)
  #summary(cox2)
  #surv=Survival(cox2)
  
  fp <- predict(cox2)
  cindex=getC_index(fp,norm.stat.al$time,norm.stat.al$status)
  #cindex.orig=1-rcorr.cens(fp,Surv(norm.stat.al$time,norm.stat.al$status))
  cut.time=c()
  if(quantile(os[!is.na(os)])['75%']<12){
    cut.time=mks
  }else if(quantile(os[!is.na(os)])['75%']<365){
    cut.time=c(12*mks[1],12*mks[2],12*mks[3])
  }else{
    cut.time=c(365*mks[1],365*mks[2],365*mks[3])
  }
  cut.time=cut.time[which(cut.time<quantile(os,seq(0,1,0.01))['100%'])]
  print(cut.time)
  #regplot(cox2)
  #  print(regplot(cox3#对观测2的六个指标在列线图上进行计分展示
  #,observation=pbc[2,] #也可以不展示
  #预测3年和5年的死亡风险，此处单位是day
  #              ,title=title
  #              ,failtime = cut.time
  #              ,prfail = TRUE #cox回归中需要TRUE
  #              ,showP = T #是否展示统计学差异
  #              ,droplines = F#观测2示例计分是否画线
  #,colors = mg_colors[1:3] #用前面自己定义的颜色
  #,rank="decreasing") #根据统计学差异的显著性进行变量的排序
  #,interval="confidence"
  #,rank="decreasing"
  #,clickable=T
  #              ,points=TRUE)) #展示观测的可信区间
  
  #  plot(nom)
  surv=Survival(cox2)
  survs=list()
  cal_all=list()
  for(i in 1:length(cut.time)){
    f1<-cph(formula = fmla,data=norm.stat.al,x=T,y=T,surv = T,na.action=na.delete,time.inc = cut.time[i]) 
    cal1<-calibrate(f1, cmethod="KM", method="boot",u=cut.time[i],m=floor(sum(f1$n)/3)) 
    cal_all=c(cal_all,list(cal1))
    #    surv0 <- function(x)surv(cut.time[i],lp=x) 
    #    survs=c(survs,list(surv0))
  }
  #f2<-cph(formula = fmla,data=norm.stat.al,x=T,y=T,surv = T,na.action=na.delete,time.inc = cut.time[2]) 
  #cal3<-calibrate(f2, cmethod="KM", method="boot",u=cut.time[2],m=100,B=200) 
  #f3<-cph(formula = fmla,data=norm.stat.al,x=T,y=T,surv = T,na.action=na.delete,time.inc = cut.time[3]) 
  #cal5<-calibrate(f3, cmethod="KM", method="boot",u=cut.time[3],m=100,B=200) 
  if(length(cut.time)==1){
    surv1 <- function(x)surv(cut.time[1],lp=x) 
    survs=list(surv1)
  }else if(length(cut.time)==2){
    surv1 <- function(x)surv(cut.time[1],lp=x) 
    surv2 <- function(x)surv(cut.time[2],lp=x) 
    survs=list(surv1,surv2)
  }else if(length(cut.time)==3){
    surv1 <- function(x)surv(cut.time[1],lp=x) 
    surv2 <- function(x)surv(cut.time[2],lp=x) 
    surv3 <- function(x)surv(cut.time[3],lp=x) 
    survs=list(surv1,surv2,surv3)
  }
  nom=nomogram(cox2,fun=survs,lp= F
               ,funlabel=c(paste0(mks[1], '-Year Survival'),
                           paste0(mks[2], '-Year Survival'),
                           paste0(mks[3], '-Year Survival'))[1:length(cut.time)]
               ,maxscale=100
               ,fun.at=seq(0,1,0.2)
  )
  
  if(!quick){
    cal_all=list()
    for(i in 1:length(cut.time)){
      cal1=get_best_calibrate(cox2,cut.time[i])
      cal_all=c(cal_all,list(cal1))
    }
    #cal3=get_best_calibrate(cox2,cut.time[2])
    #cal5=get_best_calibrate(cox2,cut.time[3])
  }
  lay2 <- customLayout::lay_new(matrix(1:2))
  lay1 <- customLayout::lay_new(matrix(1:1))
  cl <- customLayout::lay_bind_col(lay2, lay1, widths = c(1, 1.5)) 
  #customLayout::lay_show(cl)
  customLayout::lay_set(cl) 
  
  plot(cal_all[[1]],lwd = 2,lty = 1,errbar.col = mg_colors[1], bty = "l",xlim = c(0,1),ylim= c(0,1),xlab = "Nomogram-prediced OS (%)"
       ,ylab = "Observed OS (%)",col = mg_colors[1],cex.lab=1.2,cex.axis=1, cex.main=1.2, cex.sub=0.6,subtitles = F)
  #lines(cal_all[[1]][,c('mean.predicted',"KM")],type = 'b', lwd = 1, col = mg_colors[1], pch = 16)
  mtext("")
  if(length(cal_all)>1){
    for(i in 2:length(cal_all)){
      plot(cal_all[[i]],lwd = 2,lty = 1,errbar.col = mg_colors[i],xlim = c(0,1),ylim= c(0,1),col = mg_colors[i],add = T)
      #lines(cal_all[[i]][,c('mean.predicted',"KM")],type = 'b', lwd = 1, col = mg_colors[i], pch = 16)
    }
    #plot(cal3,lwd = 2,lty = 0,errbar.col = mg_colors[3],xlim = c(0,1),ylim= c(0,1),col = mg_colors[3],add = T)
    #lines(cal3[,c('mean.predicted',"KM")],type = 'b', lwd = 1, col = mg_colors[3], pch = 16)
  }
  abline(0,1, lwd = 2, lty = 3, col = 'black')
  legend("topleft", legend = c(paste0(mks[1], '-Year'),
                               paste0(mks[2], '-Year'),
                               paste0(mks[3], '-Year'))[1:length(cut.time)],col =mg_colors[1:length(cut.time)],lwd = 2,cex = 1.2,bty = "n")
  
  fp <- predict(cox2)
  cindex=getC_index(fp,norm.stat.al$time,norm.stat.al$status)
  dca_dat=data.frame(Nomogram=fp,status=norm.stat.al$status)
  fp.al=cbind(fp)
  for(i in 1:ncol(clinical_riskscore)){
    fmla1 <- as.formula(paste0("Surv(time, status) ~",colnames(clinical_riskscore)[i]))
    cox1 <- cph(fmla1, data = norm.stat.al,surv = T,x = T,y = T)
    fp1 <- predict(cox1)
    fp.al=cbind(fp.al,fp1)
  }
  colnames(fp.al)=c('Nomogram',colnames(clinical_riskscore))
  fp.al=crbind2DataFrame(fp.al)
  fp.al$status=norm.stat.al$status
  mg_plotDCA(fp.al$status
             ,c('Nomogram',colnames(clinical_riskscore))
             ,c('Nomogram',colnames(clinical_riskscore)),fp.al)
  #plot(cal1,xlim=c(0,1),ylim=c(0,1))
  #plot(cal3,xlim=c(0,1),ylim=c(0,1))
  #plot(cal5,xlim=c(0,1),ylim=c(0,1))
  plot(nom)
  options(datadist=NULL)
  return(list(Mod=cox2,Cindex=cindex,CutTime=cut.time))
}



#筛选编码蛋白基因
genecode=read.delim('GeneTag.genecode.v32.txt',sep='\t',header = T)
table(genecode$TYPE)
mrna_genecode=genecode[which(genecode$TYPE=='protein_coding'),]
head(mrna_genecode)



#####数据准备########
#####STAD########
stad_cli<-read.delim('origin_datas/TCGA/STAD/Merge_STAD_clinical.txt',sep='\t',header = T)
colnames(stad_cli)[1:30]
stad_cli=data.frame(Samples=paste0(stad_cli$A0_Samples,'-01'),
                    Age=stad_cli$A17_Age,
                    Gender=stad_cli$A18_Sex,
                    T.stage=stad_cli$A3_T,
                    N.stage=stad_cli$A4_N,
                    M.stage=stad_cli$A5_M,
                    Stage=stad_cli$A6_Stage,#Grade=stad_cli$A7_Grade,
                    OS.time=stad_cli$A1_OS,Status=stad_cli$A2_Event)
rownames(stad_cli)=stad_cli$Samples
stad_cli=stad_cli %>% drop_na(OS.time)
stad_cli$OS.time
stad_cli=stad_cli[which(stad_cli$OS.time>0),]
table(stad_cli$Status)
stad_cli$OS=ifelse(stad_cli$Status=='Alive',0,1)

table(stad_cli$Stage)
stad_cli$Stage <- gsub('[ABC]', '', stad_cli$Stage)
stad_cli$Stage <- gsub('Stage ', '', stad_cli$Stage)
stad_cli$Stage[stad_cli$Stage=='']<-NA

# table(stad_cli$Grade)
# stad_cli$Grade[stad_cli$Grade=='GX']<-NA

table(stad_cli$T.stage)
stad_cli$T.stage <- gsub('[abc]', '', stad_cli$T.stage)
stad_cli$T.stage[stad_cli$T.stage=='TX']<-NA

table(stad_cli$N.stage)
stad_cli$N.stage <- gsub('[abc]', '', stad_cli$N.stage)
stad_cli$N.stage[stad_cli$N.stage=='NX' | stad_cli$N.stage=='']<-NA

table(stad_cli$M.stage)
stad_cli$M.stage[stad_cli$M.stage=='MX']<-NA
stad_cli=crbind2DataFrame(stad_cli)

fivenum(stad_cli$Age)
stad_cli$Age1=ifelse(stad_cli$Age>67,'>67','<=67')
head(stad_cli)


##########胃癌表达谱
stad_data<-read.delim('origin_datas/TCGA/STAD/Merge_TCGA-STAD_TPM.txt',sep='\t',header = T,row.names = 1,check.names = F)
stad_data[1:4,1:4]
table(substr(colnames(stad_data),14,15))

stad.sample_T=colnames(stad_data)[which(as.numeric(substr(colnames(stad_data),14,15))==1)]#肿瘤样本
stad.sample_N=colnames(stad_data)[which(as.numeric(substr(colnames(stad_data),14,15))==11)]#正常样本
stad_type=data.frame(Samples=c(stad.sample_T,stad.sample_N),
                     Type=rep(c('Tumor','Normal'),c(length(stad.sample_T),length(stad.sample_N))))
rownames(stad_type)=stad_type$Samples
table(stad_type$Type)
# Normal  Tumor 
# 32    375 

stad_tpm_log=log2(stad_data[intersect(mrna_genecode$SYMBOL,rownames(stad_data)),c(stad.sample_T,stad.sample_N)]+1)
dim(stad_tpm_log)

stad_tpm_log_T=stad_tpm_log[,intersect(stad_cli$Samples,stad.sample_T)]
dim(stad_tpm_log_T)
stad_cli=stad_cli[intersect(stad_cli$Samples,stad.sample_T),]
dim(stad_cli)
#350  12
stad_cli$dataset='STAD'


#####COAD#####
coad_cli<-read.delim('origin_datas/TCGA/COAD+READ/Merge_COAD_clinical.txt',sep='\t',header = T)
colnames(coad_cli)[1:30]
coad_cli=data.frame(Samples=paste0(coad_cli$A0_Samples,'-01'),
                    Age=coad_cli$A17_Age,
                    Gender=coad_cli$A18_Sex,
                    T.stage=coad_cli$A3_T,
                    N.stage=coad_cli$A4_N,
                    M.stage=coad_cli$A5_M,
                    Stage=coad_cli$A6_Stage,
                    OS.time=coad_cli$A1_OS,Status=coad_cli$A2_Event)
rownames(coad_cli)=coad_cli$Samples
coad_cli=coad_cli %>% drop_na(OS.time)
coad_cli$OS.time
coad_cli=coad_cli[which(coad_cli$OS.time>0),]
table(coad_cli$Status)
coad_cli$OS=ifelse(coad_cli$Status=='Alive',0,1)

table(coad_cli$Stage)
coad_cli$Stage <- gsub('[ABC]', '', coad_cli$Stage)
coad_cli$Stage <- gsub('Stage ', '', coad_cli$Stage)
coad_cli$Stage[coad_cli$Stage=='']<-NA

table(coad_cli$T.stage)
coad_cli$T.stage <- gsub('[ab]', '', coad_cli$T.stage)
coad_cli$T.stage[coad_cli$T.stage=='Tis']<-NA

table(coad_cli$N.stage)
coad_cli$N.stage <- gsub('[abc]', '', coad_cli$N.stage)
coad_cli$N.stage[coad_cli$N.stage=='NX']<-NA

table(coad_cli$M.stage)
coad_cli$M.stage <- gsub('[ab]', '', coad_cli$M.stage)
coad_cli$M.stage[coad_cli$M.stage=='MX'|coad_cli$M.stage=='']<-NA
coad_cli=crbind2DataFrame(coad_cli)

fivenum(coad_cli$Age)
coad_cli$Age1=ifelse(coad_cli$Age>68,'>68','<=68')
head(coad_cli)


##########COAD表达谱
coad_data<-read.delim('origin_datas/TCGA/COAD+READ/COAD_TPM.txt',sep='\t',header = T,row.names = 1,check.names = F)
coad_data[1:4,1:4]
table(substr(colnames(coad_data),14,15))

coad.sample_T=colnames(coad_data)[which(as.numeric(substr(colnames(coad_data),14,15))==1)]#肿瘤样本
coad.sample_N=colnames(coad_data)[which(as.numeric(substr(colnames(coad_data),14,15))==11)]#正常样本
coad_type=data.frame(Samples=c(coad.sample_T,coad.sample_N),
                     Type=rep(c('Tumor','Normal'),c(length(coad.sample_T),length(coad.sample_N))))
rownames(coad_type)=coad_type$Samples
table(coad_type$Type)
# Normal  Tumor 
#  41    456 

coad_tpm_log=log2(coad_data[intersect(mrna_genecode$SYMBOL,rownames(coad_data)),c(coad.sample_T,coad.sample_N)]+1)
dim(coad_tpm_log)

coad_tpm_log_T=coad_tpm_log[,intersect(coad_cli$Samples,coad.sample_T)]
dim(coad_tpm_log_T)
coad_cli=coad_cli[intersect(coad_cli$Samples,coad.sample_T),]
dim(coad_cli)
# 432  11
coad_cli$dataset='COAD'

#####READ#####
read_cli<-read.delim('origin_datas/TCGA/COAD+READ/Merge_READ_clinical.txt',sep='\t',header = T)
colnames(read_cli)[1:30]
read_cli=data.frame(Samples=paste0(read_cli$A0_Samples,'-01'),
                    Age=read_cli$A17_Age,
                    Gender=read_cli$A18_Sex,
                    T.stage=read_cli$A3_T,
                    N.stage=read_cli$A4_N,
                    M.stage=read_cli$A5_M,
                    Stage=read_cli$A6_Stage,
                    OS.time=read_cli$A1_OS,Status=read_cli$A2_Event)
rownames(read_cli)=read_cli$Samples
read_cli=read_cli %>% drop_na(OS.time)
read_cli$OS.time
read_cli=read_cli[which(read_cli$OS.time>0),]
table(read_cli$Status)
read_cli$OS=ifelse(read_cli$Status=='Alive',0,1)

table(read_cli$Stage)
read_cli$Stage <- gsub('[ABC]', '', read_cli$Stage)
read_cli$Stage <- gsub('Stage ', '', read_cli$Stage)
read_cli$Stage[read_cli$Stage=='']<-NA


table(read_cli$T.stage)
read_cli$T.stage <- gsub('[ab]', '', read_cli$T.stage)
read_cli$T.stage[read_cli$T.stage=='']<-NA

table(read_cli$N.stage)
read_cli$N.stage <- gsub('[abc]', '', read_cli$N.stage)
read_cli$N.stage[read_cli$N.stage=='NX' | read_cli$N.stage=='']<-NA

table(read_cli$M.stage)
read_cli$M.stage <- gsub('[a]', '', read_cli$M.stage)
read_cli$M.stage[read_cli$M.stage=='MX'|read_cli$M.stage=='']<-NA
read_cli=crbind2DataFrame(read_cli)

fivenum(read_cli$Age)
read_cli$Age1=ifelse(read_cli$Age>66,'>66','<=66')
head(read_cli)


##########READ表达谱
read_data<-read.delim('origin_datas/TCGA/COAD+READ/READ_TPM.txt',sep='\t',header = T,row.names = 1,check.names = F)
read_data[1:4,1:4]
table(substr(colnames(read_data),14,15))

read.sample_T=colnames(read_data)[which(as.numeric(substr(colnames(read_data),14,15))==1)]#肿瘤样本
read.sample_N=colnames(read_data)[which(as.numeric(substr(colnames(read_data),14,15))==11)]#正常样本
read_type=data.frame(Samples=c(read.sample_T,read.sample_N),
                     Type=rep(c('Tumor','Normal'),c(length(read.sample_T),length(read.sample_N))))
rownames(read_type)=read_type$Samples
table(read_type$Type)
# Normal  Tumor 
# 10    166 

read_tpm_log=log2(read_data[intersect(mrna_genecode$SYMBOL,rownames(read_data)),c(read.sample_T,read.sample_N)]+1)
dim(read_tpm_log)

read_tpm_log_T=read_tpm_log[,intersect(read_cli$Samples,read.sample_T)]
dim(read_tpm_log_T)
read_cli=read_cli[intersect(read_cli$Samples,read.sample_T),]
dim(read_cli)
#157  11
read_cli$dataset='READ'

######合并后PCA########
com.genes=Reduce(intersect,list(rownames(stad_tpm_log),rownames(coad_tpm_log),rownames(read_tpm_log)))
length(com.genes)

tcga_exp=cbind(stad_tpm_log[com.genes,],coad_tpm_log[com.genes,],read_tpm_log[com.genes,])
dim(tcga_exp)

tcga_cli=rbind.data.frame(stad_cli,coad_cli,read_cli)
fivenum(tcga_cli$Age)
tcga_cli$Age1=ifelse(tcga_cli$Age>67,'>67','<=67')
head(tcga_cli)
table(tcga_cli$dataset)
# COAD READ STAD 
# 432  157  350 
writeMatrix(tcga_cli,'results/TCGA_clinical.txt')

stad_type$dataset='STAD'
coad_type$dataset='COAD'
read_type$dataset='READ'
tcga_type=rbind(stad_type,coad_type,read_type)
table(tcga_type$Type)
table(tcga_type$dataset)
# COAD READ STAD 
# 497  176  407 

library(ggbiplot)
before.pca <- prcomp(t(tcga_exp[apply(tcga_exp,1,function(x){return(sd(x)>0.5)}),]), scale=T)
pca_before <- ggbiplot(before.pca, scale=1, groups = tcga_type$dataset,
                       ellipse = TRUE,ellipse.prob=0.5, circle = F,var.axes=F) +
  scale_color_manual(values = ggsci::pal_lancet('lanonc')(9)) + 
  # xlim(-5, 5) + ylim(-5,5) +
  theme_light() +
  theme(legend.direction = 'horizontal', legend.position = 'top') +
  xlab('PCA1') + ylab('PCA2')
pca_before
#去除批次效应
library(sva)
library(limma)
## 使用 limma 的 removeBatchEffect 函数
tcga_exp <- removeBatchEffect(tcga_exp,batch =  tcga_type$dataset)
dim(tcga_exp)
#boxplot(tcga_exp[, 1:10], las=2)
tcga_tpm_log <- normalizeBetweenArrays(tcga_exp)
dim(tcga_tpm_log)
# boxplot(tcga_tpm_log_T[, 1:10], las=2)

after.pca <- prcomp(t(tcga_tpm_log), scale=T)
pca_after <- ggbiplot(after.pca, scale=1, groups =  tcga_type$dataset,
                      ellipse = TRUE,ellipse.prob=0.5, circle = F,var.axes=F) +
  scale_color_manual(values = ggsci::pal_lancet('lanonc')(9)) + 
  # xlim(-5, 5) + ylim(-5,5) +
  theme_light() +
  theme(legend.direction = 'horizontal', legend.position = 'top') +
  xlab('PCA1') + ylab('PCA2')
pca_after

tcga_tpm_log_T=tcga_tpm_log[,tcga_cli$Samples]
dim(tcga_tpm_log_T)
# 19004   939

pca.plot=mg_merge_plot(pca_before,pca_after,labels = c('A','B'))
savePDF('results/PCA.plot.pdf',pca.plot,height = 5,width = 10)


############01.鉴定胆汁酸相关分子亚型##############
dir.create('results/01.Cluster')
######胆汁酸相关基因
bile.acid.genes=read.gmt('origin_datas/HALLMARK_BILE_ACID_METABOLISM.v2022.1.Hs.gmt')
bile.acid.genes=bile.acid.genes$gene
length(bile.acid.genes)
#112
setdiff(bile.acid.genes,rownames(tcga_tpm_log_T))

####胆汁酸预后相关基因
ba.gene.cox=cox_batch(dat = tcga_tpm_log_T[intersect(bile.acid.genes,rownames(tcga_tpm_log_T)),tcga_cli$Samples],
          time = tcga_cli$OS.time,event = tcga_cli$OS)
table(ba.gene.cox$p.value<0.05)
# FALSE  TRUE 
# 96    16 
ba.gene.cox.fit=ba.gene.cox[ba.gene.cox$p.value<0.05,]
data.frame(Gene=rownames(ba.gene.cox.fit),round(ba.gene.cox.fit,3))
pdf('results/01.Cluster/Fig1a.pdf',height = 7,width = 7,onefile = F)
mg_forestplot_v2(data.frame(Gene=rownames(ba.gene.cox.fit[order(ba.gene.cox.fit$HR),]),
                            round(ba.gene.cox.fit[order(ba.gene.cox.fit$HR),],3)),
                 xlog = T,colgap = 15,lineheight = 10,zero = 1,
                 boxsize = 0.3,lwd.zero=1,lwd.ci=1,lwd.xaxis=1
                 ,box_col='blue',summary_col="black",lines_col='blue',zero_col='grey'
                 ,xlab='Hazard Ratio',lty.ci = 10,graph.pos = 3)
dev.off()

write.csv(round(ba.gene.cox.fit,3),'results/01.Cluster/BILE_ACID_GENE_COX.csv')

##############一致性聚类#######################
library(ConsensusClusterPlus)
clusterAlg_name=c('hc','pam','km','kmdist')[3]
distance_name=c('pearson','spearman','euclidean','binary','maximum','canberra','minkowski')[3]
consen_gene=rownames(ba.gene.cox.fit)
length(consen_gene)

tcga_consen_data=as.matrix(tcga_tpm_log_T[intersect(consen_gene,rownames(tcga_tpm_log_T)),])
tcga_consen_data=t(scale(t(tcga_consen_data),scale = F))   
tcga_consen_data=as.matrix(tcga_consen_data)
dim(tcga_consen_data)
tcga_clust_subtype <- ConsensusClusterPlus(tcga_consen_data
                                           , maxK = 10, reps = 500, pItem = 0.8
                                           , pFeature = 1
                                           , title = "TCGA_subtype"
                                           , clusterAlg = clusterAlg_name
                                           , distance = distance_name
                                           , plot = "pdf"
                                           , writeTable = T
                                           , seed = 123456)
k=3
cluster.color=pal_jama()(7)[4:6]
tcga.subtype <- data.frame(Samples = names(tcga_clust_subtype[[k]]$consensusClass),
                           Cluster=tcga_clust_subtype[[k]]$consensusClass)
tcga.subtype$Cluster=paste0('C',tcga.subtype$Cluster)
table(tcga.subtype$Cluster)
# C1  C2  C3 
#401 315 223 
fig1b=ggplotKMCox(data.frame(time = tcga_cli[rownames(tcga.subtype),]$OS.time/365
                       , event = tcga_cli[rownames(tcga.subtype),]$OS
                       , tcga.subtype$Cluster)
            ,add_text = '',title = 'TCGA',show_confint = F,
            palette = cluster.color,labs = names(table(tcga.subtype$Cluster)))
fig1b


###########分子特征
tcga.molecular.sign=readMatrix(paste0(MG_Grobal_baseFolder,'/source/PMC5982584_supplement_2.txt'))
table(tcga.molecular.sign$`TCGA Study`)
head(tcga.molecular.sign)
tcga.molecular.sign$Samples=paste0(rownames(tcga.molecular.sign),'-01')
rownames(tcga.molecular.sign)=tcga.molecular.sign$Samples
tcga.molecular.sign=merge(tcga.molecular.sign,tcga.subtype,by='Samples')
dim(tcga.molecular.sign)
# 403  63
head(tcga.molecular.sign)
fig1c=plotMutiBar(table(tcga.molecular.sign$`Immune Subtype`,tcga.molecular.sign$Cluster))

fig1d=ggplotKMCox(data.frame(time = tcga.molecular.sign$`OS Time`/365
                       , event = tcga.molecular.sign$OS
                       , tcga.molecular.sign$`Immune Subtype`)
            ,add_text = '',title = 'Immune Subtype',show_confint = F,
            labs = names(table(tcga.molecular.sign$`Immune Subtype`)))


tcga.subtype.cli=merge(tcga_cli,tcga.subtype,bu='Samples')
rownames(tcga.subtype.cli)=tcga.subtype.cli$Samples
head(tcga.subtype.cli)

fig1e=plotMutiBar(dat = table(tcga.subtype.cli$Status,tcga.subtype.cli$Cluster) )

fig1bcde=mg_merge_plot(fig1b,fig1e,fig1c,fig1d,labels = LETTERS[5:8])
savePDF('results/01.Cluster/Fig1bcde.pdf',fig1bcde,height = 10,width = 10)


#####02.亚型间免疫特征########
dir.create('results/02.subtype.immu')
# tcga.ciber=immu_CIBERSORT(exp_data = tcga_tpm_log_T)
# save(tcga.ciber,file = 'results/tcga.ciber.RData')
load('results/tcga.ciber.RData')
fig2a=mg_PlotMutiBoxplot(tcga.ciber[tcga.subtype$Samples,1:22]
                   , group = tcga.subtype$Cluster
                   , legend.pos = 'top'
                   , group_cols = cluster.color
                   , test_method =  c('kruskal.test','wilcox.test','anova')[1]
                   , add = 'boxplot'
                   , ylab = 'Fraction')

# tcga.est=immu_estimate(tcga_tpm_log_T)
# save(tcga.est,file = 'results/tcga.est.RData')
load('results/tcga.est.RData')
tcga.est.use=cbind.data.frame(tcga.est[tcga.subtype$Samples,],Cluster=tcga.subtype$Cluster)
head(tcga.est.use)
tcga.est.use=melt(tcga.est.use)
library(ggridges)
fig2b=ggplot(tcga.est.use, aes(x=value, y=variable, color=Cluster, point_color=Cluster, fill=Cluster)) +
  geom_density_ridges(
    jittered_points=TRUE, scale = .95, rel_min_height = .01,
    point_shape = "|", point_size = 3, size = 0.25,
    position = position_points_jitter(height = 0)) +
  scale_y_discrete(expand = c(.01, 0)) +
  scale_x_continuous(expand = c(0, 0), name = "Score") +
  scale_fill_manual(values = cluster.color, labels = c("C1", "C2","C3")) +
  scale_color_manual(values = cluster.color, guide = "none") +
  scale_discrete_manual("point_color", values = cluster.color, guide = "none") +
  guides(fill = guide_legend(
    override.aes = list(
      fill = cluster.color,
      color = NA, point_color = NA))) +
  #ggtitle("Height in Australian athletes") +
  theme_ridges(center = TRUE)+theme(legend.position = 'top')

fig2ab=mg_merge_plot(fig2a,fig2b,widths = c(2.5,1),labels = LETTERS[1:2])
savePDF('results/02.subtype.immu/Fig2ab.pdf',fig2ab,height = 5,width = 15)

###########29种免疫特征：PMC6310928
# immune.signature=readxl::read_xlsx('origin_datas/29_immune_signature_PMC6310928.xlsx',col_names = T)
# immune.signature=crbind2DataFrame(immune.signature)
# immune.signature=t(immune.signature)
# immune.signature=melt(immune.signature)
# immune.signature=data.frame(immu_signature=immune.signature$Var1,gene=immune.signature$value)
# immune.signature=immune.signature %>% drop_na(gene)
# immune.genesets.list=split(x=immune.signature,f=immune.signature$immu_signature)
# immune.genesets.list=sapply(immune.genesets.list, function(x){subset(x,select='gene',drop=TRUE)})
# 
# tcga.tme.ssgsea=ssGSEAScore_by_muti_group_genes(gene.exp = tcga_tpm_log_T
#                                                 ,genelist = immune.genesets.list)
# save(tcga.tme.ssgsea,file='results/tcga.tme.ssgsea.RData')
load('results/tcga.tme.ssgsea.RData')
tcga.tme.ssgsea[1:5,1:5]

cli_anno_sub=data.frame(Cluster=tcga.subtype[order(tcga.subtype$Cluster),'Cluster'])
rownames(cli_anno_sub)=tcga.subtype[order(tcga.subtype$Cluster),'Samples']

cluster.color.use=cluster.color
names(cluster.color.use)=names(table(tcga.subtype$Cluster))
table(tcga.subtype$Cluster)
pdf('results/02.subtype.immu/Fig2c.pdf',height = 6,width = 8)
pheatmap(tcga.tme.ssgsea[,rownames(cli_anno_sub)],
         scale = 'row',breaks = c(-3,0,3),
         color =  colorRampPalette(c('purple', "grey21", 'Yellow'))(100),
         cluster_rows = F,cluster_cols = F,
         show_colnames = F,#filename = "PDFs/test.pdf",
         gaps_col = c(401,401+315),
         annotation_col = cli_anno_sub,
         annotation_colors = list(Cluster=cluster.color.use))
dev.off()


# tcga.immu.ssgsea=immu_ssgsea(tcga_tpm_log_T)
# save(tcga.immu.ssgsea,file = 'results/tcga.immu.ssgsea.RData')
load('results/tcga.immu.ssgsea.RData')
pdf('results/02.subtype.immu/Fig2d.pdf',height = 6,width = 8)
pheatmap(t(tcga.immu.ssgsea[rownames(cli_anno_sub),]),
         scale = 'row',breaks = c(-3,0,3),
         color =  colorRampPalette(c('DeepSkyBlue', "grey21", 'Yellow'))(100),
         #color = circlize::colorRamp2(c(-3, 0, 3), c('blue', 'white', 'red')),
         cluster_rows = F,cluster_cols = F,
         show_colnames = F,#filename = "PDFs/test.pdf",
         gaps_col = c(401,401+315),
         annotation_col = cli_anno_sub,
         annotation_colors = list(Cluster=cluster.color.use))
dev.off()


########血管生成得分
tcga.Angiogenesis=immu_AngiogenesisScore(tcga_tpm_log_T)
tcga.Angiogenesis=data.frame(Samples=names(tcga.Angiogenesis),Angiogenesis_score=as.numeric(tcga.Angiogenesis))
rownames(tcga.Angiogenesis)=tcga.Angiogenesis$Samples
fig2e=my_boxplot(dat = tcga.Angiogenesis[tcga.subtype$Samples,'Angiogenesis_score'],
           group = tcga.subtype$Cluster,
           group_cols = cluster.color,legend.position = 'none',
           ylab = 'Angiogenesis score',fill = 'Cluster',
           test_method =  c('kruskal.test','wilcox.test')[1])

############炎症相关通路
t_cell=read.gmt('origin_datas/pathway/KEGG_T_CELL_RECEPTOR_SIGNALING_PATHWAY.v2022.1.Hs.gmt')
tcga.t.cell=t(ssGSEAScore_by_genes(tcga_tpm_log_T,t_cell$gene))
tcga.t.cell=crbind2DataFrame(tcga.t.cell)
fig2f=my_boxplot(dat = tcga.t.cell$GeneSet,
           group = tcga.subtype$Cluster,
           group_cols = cluster.color,legend.position = 'none',
           ylab = 'T cell receptor signaling pathway',fill = 'Cluster',
           test_method =  c('kruskal.test','wilcox.test')[1])


b_cell=read.gmt('origin_datas/pathway/KEGG_B_CELL_RECEPTOR_SIGNALING_PATHWAY.v2022.1.Hs.gmt')
tcga.b.cell=t(ssGSEAScore_by_genes(tcga_tpm_log_T,b_cell$gene))
tcga.b.cell=crbind2DataFrame(tcga.b.cell)
fig2g=my_boxplot(dat = tcga.b.cell$GeneSet,
           group = tcga.subtype$Cluster,
           group_cols = cluster.color,legend.position = 'none',
           ylab = 'B cell receptor signaling pathway',fill = 'Cluster',
           test_method =  c('kruskal.test','wilcox.test')[1])

toll_like=read.gmt('origin_datas/pathway/KEGG_TOLL_LIKE_RECEPTOR_SIGNALING_PATHWAY.v2022.1.Hs.gmt')
tcga.toll.like=t(ssGSEAScore_by_genes(tcga_tpm_log_T,toll_like$gene))
tcga.toll.like=crbind2DataFrame(tcga.toll.like)
fig2h=my_boxplot(dat = tcga.toll.like$GeneSet,
           group = tcga.subtype$Cluster,
           group_cols = cluster.color,legend.position = 'none',
           ylab = 'TOLL-like receptor signaling pathway',fill = 'Cluster',
           test_method =  c('kruskal.test','wilcox.test')[1])

jak_stat=read.gmt('origin_datas/pathway/KEGG_JAK_STAT_SIGNALING_PATHWAY.v2022.1.Hs.gmt')
tcga.jak.stat=t(ssGSEAScore_by_genes(tcga_tpm_log_T,jak_stat$gene))
tcga.jak.stat=crbind2DataFrame(tcga.jak.stat)
fig2i=my_boxplot(dat = tcga.jak.stat$GeneSet,
           group = tcga.subtype$Cluster,
           group_cols = cluster.color,legend.position = 'none',
           ylab = 'Jak-stat signaling pathway',fill = 'Cluster',
           test_method =  c('kruskal.test','wilcox.test')[1])

nfkb=read.gmt('origin_datas/pathway/BIOCARTA_NFKB_PATHWAY.v2022.1.Hs.gmt')
tcga.nfkb=t(ssGSEAScore_by_genes(tcga_tpm_log_T,nfkb$gene))
tcga.nfkb=crbind2DataFrame(tcga.nfkb)
fig2j=my_boxplot(dat = tcga.nfkb$GeneSet,
           group = tcga.subtype$Cluster,
           group_cols = cluster.color,legend.position = 'none',
           ylab = 'Nfkb signaling pathway',fill = 'Cluster',
           test_method =  c('kruskal.test','wilcox.test')[1])

fig2ei=mg_merge_plot(fig2e,fig2f,fig2g,fig2h,fig2i,fig2j,ncol=2,nrow=3,common.legend = T,labels = LETTERS[5:10])
savePDF('results/02.subtype.immu/Fig2ei.pdf',fig2ei,height = 12,width = 8)

##########03.亚型间免疫治疗#########
dir.create('results/03.subtype.immune.treatment')
############免疫检查点
tcga.icgs=immu_ICGs(tcga_tpm_log_T)
icg.genes=c('PDCD1','CD274','CTLA4','LAG3','PDCD1LG2','BTLA','HAVCR2','TIGIT')
library(tinyarray)
pdf('results/03.subtype.immune.treatment/Fig3a.pdf',height = 4,width = 7)
draw_heatmap(t(tcga.icgs[order(tcga.subtype$Cluster),icg.genes]),
             group_list = factor(tcga.subtype$Cluster[order(tcga.subtype$Cluster)]),
             cluster_cols = F,
             show_rownames = T,
             legend = T,
             annotation_legend = T,
             color_an = cluster.color)
dev.off()

###########TIDE
# tcga_tide_dat <- t(scale(t(tcga_tpm_log_T),scale = F))
# dim(tcga_tide_dat)
# write.table(tcga_tide_dat,file = 'results/tcga_tide_dat.txt',quote = F, sep = '\t')

tcga_tide_res<-read.csv('results/tcga.tide.res.csv',row.names = 1,stringsAsFactors = F)
head(tcga_tide_res)
tide_sel=c('TIDE','Exclusion','Dysfunction','MDSC','CAF')
pdf('results/03.subtype.immune.treatment/Fig3b.pdf',height = 4,width = 8,onefile = F)
my_mutiboxplot_seg(dat = tcga_tide_res[tcga.subtype$Samples,tide_sel],
                   group = tcga.subtype$Cluster,
                   test_method = kruskal.test,
                   ylab = 'Score',group_cols =cluster.color,
                   legend.pos = 'none',nrow = 1,ncol = length(tide_sel),
                   xangle = 0,xhjust = 0.5)+theme_pubclean()
dev.off()



############04.亚型间差异激活的通路#############
dir.create('results/04.subtype.pathway')
tcga.subtype.use1=tcga.subtype
tcga.subtype.use1$Cluster1=ifelse(tcga.subtype.use1$Cluster=='C1','C1','Other')
tcga.subtype.use1$Cluster2=ifelse(tcga.subtype.use1$Cluster=='C2','C2','Other')
tcga.subtype.use1$Cluster3=ifelse(tcga.subtype.use1$Cluster=='C3','C3','Other')

tcga.geneList1=getGeneFC(gene.exp=tcga_tpm_log_T[,tcga.subtype.use1$Samples],
                         group=tcga.subtype.use1$Cluster1,
                         ulab='C1',
                         dlab='Other')
tcga.geneList2=getGeneFC(gene.exp=tcga_tpm_log_T[,tcga.subtype.use1$Samples],
                         group=tcga.subtype.use1$Cluster2,
                         ulab='C2',
                         dlab='Other')
tcga.geneList3=getGeneFC(gene.exp=tcga_tpm_log_T[,tcga.subtype.use1$Samples],
                         group=tcga.subtype.use1$Cluster3,
                         ulab='C3',
                         dlab='Other')


tcga.hallmark.gsea1<-GSEA(tcga.geneList1,TERM2GENE = h.all.gmt,seed=T)
tcga.hallmark.gsea2<-GSEA(tcga.geneList2,TERM2GENE = h.all.gmt,seed=T)
tcga.hallmark.gsea3<-GSEA(tcga.geneList3,TERM2GENE = h.all.gmt,seed=T)

library(enrichplot)
library(ggplot2)


tcga.hallmark.gsea.res1=tcga.hallmark.gsea1@result
tcga.hallmark.gsea.res2=tcga.hallmark.gsea2@result
tcga.hallmark.gsea.res3=tcga.hallmark.gsea3@result
write.csv(tcga.hallmark.gsea.res1,'results/04.subtype.pathway/c1_vs_other_GSEA_res.csv')
write.csv(tcga.hallmark.gsea.res2,'results/04.subtype.pathway/c2_vs_other_GSEA_res.csv')
write.csv(tcga.hallmark.gsea.res3,'results/04.subtype.pathway/c3_vs_other_GSEA_res.csv')


rownames(tcga.hallmark.gsea.res1)=gsub("HALLMARK_","",rownames(tcga.hallmark.gsea.res1))
rownames(tcga.hallmark.gsea.res2)=gsub("HALLMARK_","",rownames(tcga.hallmark.gsea.res2))
rownames(tcga.hallmark.gsea.res3)=gsub("HALLMARK_","",rownames(tcga.hallmark.gsea.res3))

tcga.hallmark.union=Reduce(union,list(rownames(tcga.hallmark.gsea.res1),
                                      rownames(tcga.hallmark.gsea.res2),
                                      rownames(tcga.hallmark.gsea.res3)))
length(tcga.hallmark.union)#41

tcga.hallmark.heatmap.dat=matrix(0,nrow = 3,ncol = length(tcga.hallmark.union))
rownames(tcga.hallmark.heatmap.dat)=c('C1', 'C2','C3')
colnames(tcga.hallmark.heatmap.dat)=tcga.hallmark.union

tcga.hallmark.heatmap.dat[1,match(rownames(tcga.hallmark.gsea.res1),colnames(tcga.hallmark.heatmap.dat))]=tcga.hallmark.gsea.res1$NES
tcga.hallmark.heatmap.dat[2,match(rownames(tcga.hallmark.gsea.res2),colnames(tcga.hallmark.heatmap.dat))]=tcga.hallmark.gsea.res2$NES
tcga.hallmark.heatmap.dat[3,match(rownames(tcga.hallmark.gsea.res3),colnames(tcga.hallmark.heatmap.dat))]=tcga.hallmark.gsea.res3$NES
range(tcga.hallmark.heatmap.dat)

pdf('results/04.subtype.pathway/Fig4.pdf',height = 10,width = 8)
pheatmap(t(tcga.hallmark.heatmap.dat),
         #scale = 'row',
         border="white", # 设置边框为白色
         color =  colorRampPalette(c("navy", "white", "red"))(200),
         display_numbers = F, # 热图上显示数值
         cluster_cols = F, # 去掉横向、纵向聚类
         cluster_rows = T,
         show_rownames = T, #去掉横、纵坐标id
         show_colnames = T,
         #legend_breaks=c(-5,0,5),
         fontsize_row = 10, # 分别设置横向和纵向字体大小
         fontsize_col = 12)
dev.off()

##########05.亚型间差异基因的鉴定#######
dir.create('results/05.subtype.DEGs')
p_fit=0.05
fc_fit=log2(1.5)
tcga.degs.c1=mg_limma_DEG(exp = tcga_tpm_log_T,group = as.character(tcga.subtype.use1$Cluster1),ulab = 'C1',dlab = 'Other')
tcga.degs.c1$Summary
#"222|92" 
write.csv(round(tcga.degs.c1$DEG[which(tcga.degs.c1$DEG$adj.P.Val<p_fit & abs(tcga.degs.c1$DEG$logFC)>fc_fit),],3),
            'results/05.subtype.DEGs/tcga.c1.degs.csv')
tcga.degs.c1.res=rownames(tcga.degs.c1$DEG[which(tcga.degs.c1$DEG$adj.P.Val<p_fit & abs(tcga.degs.c1$DEG$logFC)>fc_fit),])
length(tcga.degs.c1.res)
#314
tcga.degs.c1$DEG$type = factor(ifelse(tcga.degs.c1$DEG$adj.P.Val<p_fit & abs(tcga.degs.c1$DEG$logFC) > fc_fit, 
                                      ifelse(tcga.degs.c1$DEG$logFC> fc_fit ,'Up','Down'),'No Signif'),levels=c('Up','Down','No Signif'))
head(tcga.degs.c1$DEG)
fig5a=ggplot(tcga.degs.c1$DEG,aes(x=logFC,y=-log10(adj.P.Val),color=type))+
  geom_point()+
  scale_color_manual(values=c("#DC143C","#00008B","#808080"))+#确定点的颜色
  theme_bw()+#修改图片背景
  theme(legend.title = element_blank())+
  ylab('-log10 (adj.PVal)')+#修改y轴名称
  xlab('log2 (FoldChange)')+#修改x轴名称
  ggtitle('C1 vs Other')+
  geom_vline(xintercept=c(-fc_fit,fc_fit),lty=3,col="black",lwd=0.5) +
  geom_hline(yintercept = -log10(p_fit),lty=3,col="black",lwd=0.5)+coord_flip()
fig5a
ggsave('fig5a.pdf',fig5a,height = 5,width = 5.5)

tcga.degs.c2=mg_limma_DEG(exp = tcga_tpm_log_T,group = as.character(tcga.subtype.use1$Cluster2),ulab = 'C2',dlab = 'Other')
tcga.degs.c2$Summary
#"28|71"
write.csv(round(tcga.degs.c2$DEG[which(tcga.degs.c2$DEG$adj.P.Val<p_fit & abs(tcga.degs.c2$DEG$logFC)>fc_fit),],3),
            'results/05.subtype.DEGs/tcga.c2.degs.csv')
tcga.degs.c2.res=rownames(tcga.degs.c2$DEG[which(tcga.degs.c2$DEG$adj.P.Val<p_fit & abs(tcga.degs.c2$DEG$logFC)>fc_fit),])
length(tcga.degs.c2.res)
#99
tcga.degs.c2$DEG$type = factor(ifelse(tcga.degs.c2$DEG$adj.P.Val<p_fit & abs(tcga.degs.c2$DEG$logFC) > fc_fit, 
                                      ifelse(tcga.degs.c2$DEG$logFC> fc_fit ,'Up','Down'),'No Signif'),levels=c('Up','Down','No Signif'))
head(tcga.degs.c2$DEG)
fig5b=ggplot(tcga.degs.c2$DEG,aes(x=logFC,y=-log10(adj.P.Val),color=type))+
  geom_point()+
  scale_color_manual(values=c("#Dc243C","#00008B","#808080"))+#确定点的颜色
  theme_bw()+#修改图片背景
  theme(legend.title = element_blank())+
  ylab('-log10 (adj.PVal)')+#修改y轴名称
  xlab('log2 (FoldChange)')+#修改x轴名称
  ggtitle('C2 vs Other')+
  geom_vline(xintercept=c(-fc_fit,fc_fit),lty=3,col="black",lwd=0.5) +
  geom_hline(yintercept = -log10(p_fit),lty=3,col="black",lwd=0.5)+coord_flip()
ggsave('fig5b.pdf',fig5b,height = 5,width = 5.5)

tcga.degs.c3=mg_limma_DEG(exp = tcga_tpm_log_T,group = as.character(tcga.subtype.use1$Cluster3),ulab = 'C3',dlab = 'Other')
tcga.degs.c3$Summary
#"389|192"
write.csv(round(tcga.degs.c3$DEG[which(tcga.degs.c3$DEG$adj.P.Val<p_fit & abs(tcga.degs.c3$DEG$logFC)>fc_fit),],3),
            'results/05.subtype.DEGs/tcga.c3.degs.csv')
tcga.degs.c3.res=rownames(tcga.degs.c3$DEG[which(tcga.degs.c3$DEG$adj.P.Val<p_fit & abs(tcga.degs.c3$DEG$logFC)>fc_fit),])
length(tcga.degs.c3.res)
#581
tcga.degs.c3$DEG$type = factor(ifelse(tcga.degs.c3$DEG$adj.P.Val<p_fit & abs(tcga.degs.c3$DEG$logFC) > fc_fit, 
                                      ifelse(tcga.degs.c3$DEG$logFC> fc_fit ,'Up','Down'),'No Signif'),levels=c('Up','Down','No Signif'))
head(tcga.degs.c3$DEG)
fig5c=ggplot(tcga.degs.c3$DEG,aes(x=logFC,y=-log10(adj.P.Val),color=type))+
  geom_point()+
  scale_color_manual(values=c("#Dc343C","#00008B","#808080"))+#确定点的颜色
  theme_bw()+#修改图片背景
  theme(legend.title = element_blank())+
  ylab('-log10 (adj.PVal)')+#修改y轴名称
  xlab('log2 (FoldChange)')+#修改x轴名称
  ggtitle('C3 vs Other')+
  geom_vline(xintercept=c(-fc_fit,fc_fit),lty=3,col="black",lwd=0.5) +
  geom_hline(yintercept = -log10(p_fit),lty=3,col="black",lwd=0.5)+coord_flip()
ggsave('fig5c.pdf',fig5c,height = 5,width = 5.5)

tcga.sub.degs=unique(c(tcga.degs.c1.res,tcga.degs.c2.res,tcga.degs.c3.res))
length(tcga.sub.degs)
#748
writeMatrix(tcga.sub.degs,'results/05.subtype.DEGs/tcga.sub.degs.txt',header = F)

########差异基因富集分析
library(topGO)
library(clusterProfiler)
library(Rgraphviz)
library(org.Hs.eg.db)
DEG.entrez_id = mapIds(x = org.Hs.eg.db,
                       keys = tcga.sub.degs,
                       keytype = "SYMBOL",
                       column = "ENTREZID")
DEG.entrez_id = na.omit(DEG.entrez_id)

#####GO分析
erich.go.BP = enrichGO(gene = tcga.sub.degs,
                       OrgDb = org.Hs.eg.db,
                       keyType = "SYMBOL",
                       ont = "BP",
                       pvalueCutoff = 0.05)
table(erich.go.BP@result$p.adjust<0.05)
erich.go.BP@result=erich.go.BP@result[erich.go.BP@result$p.adjust<0.05,]
write.csv(erich.go.BP@result,'results/05.subtype.DEGs/BP_enrich_res.csv')
fig5d=enrichplot::dotplot(erich.go.BP)+ggtitle('BP')
ggsave('fig5d.pdf',fig5d,height = 6,width = 5.5)

erich.go.MF = enrichGO(gene = tcga.sub.degs,
                       OrgDb = org.Hs.eg.db,
                       keyType = "SYMBOL",
                       ont = "MF",
                       pvalueCutoff = 0.05)
table(erich.go.MF@result$p.adjust<0.05)
erich.go.MF@result=erich.go.MF@result[erich.go.MF@result$p.adjust<0.05,]
write.csv(erich.go.MF@result,'results/05.subtype.DEGs/MF_enrich_res.csv')
fig5e=enrichplot::dotplot(erich.go.MF)+ggtitle('MF')
ggsave('fig5e.pdf',fig5e,height = 6,width = 5.5)

erich.go.CC = enrichGO(gene = tcga.sub.degs,
                       OrgDb = org.Hs.eg.db,
                       keyType = "SYMBOL",
                       ont = "CC",
                       pvalueCutoff = 0.05)
table(erich.go.CC@result$p.adjust<0.05)
erich.go.CC@result=erich.go.CC@result[erich.go.CC@result$p.adjust<0.05,]
write.csv(erich.go.CC@result,'results/05.subtype.DEGs/CC_enrich_res.csv')
fig5f=enrichplot::dotplot(erich.go.CC)+ggtitle('CC')
ggsave('fig5f.pdf',fig5f,height = 6,width = 5.5)


erich.KEGG = enrichKEGG(gene = DEG.entrez_id,
                        organism = "hsa",
                        keyType = "kegg",
                        pvalueCutoff = 0.05)
table(erich.KEGG@result$p.adjust<0.05)
erich.KEGG@result=erich.KEGG@result[erich.KEGG@result$p.adjust<0.05,]
write.csv(erich.KEGG@result,'results/05.subtype.DEGs/KEGG_enrich_res.csv')
# enrichplot::dotplot(erich.KEGG)
# barplot(erich.KEGG)+ggtitle('KEGG')
fig5g=enrichplot::cnetplot(erich.KEGG,circular=T,colorEdge = T,showCategory = 10)
ggsave('fig5g.pdf',fig5g,height = 12,width = 16)

pdf('results/05.subtype.DEGs/Fig5.pdf',height = 25,width = 22,onefile = F)
mg_merge_plot(mg_merge_plot(fig5a,fig5b,fig5c,ncol=3,labels = LETTERS[1:3],common.legend = T),
              mg_merge_plot(fig5d,fig5e,fig5f,ncol=3,labels = LETTERS[4:6],widths = c(1.2,1,0.8)),
              fig5g,ncol=1,nrow=3,heights = c(1,1,2.7),labels = c('','','G'))
dev.off()



#######06.风险模型############
dir.create('results/06.module')
pre.genes=read.csv('results/06.module/748_DEGS_PPI_0.4.csv')
pre.genes=crbind2DataFrame(pre.genes)
pre.genes=pre.genes$name
length(pre.genes)
setdiff(pre.genes,rownames(tcga_tpm_log_T))


tcga_model_data=t(tcga_tpm_log_T[intersect(pre.genes,rownames(tcga_tpm_log_T)),tcga_cli$Samples])
colnames(tcga_model_data)=gsub('-','__',colnames(tcga_model_data))
tcga_model_data=merge(data.frame(Samples=tcga_cli$Samples,OS=tcga_cli$OS,OS.time=tcga_cli$OS.time),
                      data.frame(Samples=rownames(tcga_model_data),tcga_model_data),by='Samples')
rownames(tcga_model_data)=tcga_model_data$Samples
tcga_model_data=tcga_model_data[,-1]
tcga_model_data=crbind2DataFrame(tcga_model_data)
dim(tcga_model_data)

tcga.cox=cox_batch(dat = tcga_tpm_log_T[intersect(pre.genes,rownames(tcga_tpm_log_T)),tcga_cli$Samples],
                   time = tcga_cli$OS.time,event = tcga_cli$OS)
tcga.cox=na.omit(tcga.cox)
head(tcga.cox)

rownames(tcga.cox)=gsub('-','__',rownames(tcga.cox))
p_cutoff=0.05
table(tcga.cox$p.value<p_cutoff)
# FALSE  TRUE 
# 393   265 
tcga.cox.fit=tcga.cox[tcga.cox$p.value<p_cutoff,]
write.csv(round(tcga.cox.fit,3),'results/06.module/tcga.cox.csv')

tcga.lasso=get_riskscore.lasso(dat = tcga_model_data[,rownames(tcga.cox.fit)],
                               os = tcga_model_data$OS,
                               os.time = tcga_model_data$OS.time)
length(tcga.lasso$lasso.gene)#20
pdf('results/06.module/Fig6a.pdf',height = 5,width = 10)
tcga.lasso$plot
dev.off()

tcga.module.risk=get_riskscore(dat = tcga_model_data[tcga_cli$Samples,tcga.lasso$lasso.gene],
                               os = tcga_model_data[tcga_cli$Samples,]$OS,
                               os.time = tcga_model_data[tcga_cli$Samples,]$OS.time,
                               step = T,direction = c("both", "backward", "forward")[1])
length(tcga.module.risk$module.gene)
#9
tcga.module.risk$model
#"-0.125*AXIN2+-0.075*ATOH1+0.136*CHST13+0.095*PNMA2+-0.102*GYG2+0.083*MAGEA3+0.122*SNCG+0.259*HEYL+-0.07*RASSF10"
cox <- coxph(as.formula(paste0("Surv(OS.time, OS) ~"
                               ,paste0(tcga.module.risk$module.gene,collapse = '+'))), 
             data =as.data.frame(tcga_model_data))
coef(cox)
module.gene.cox=data.frame(Gene=names(coef(cox)),coef=as.numeric(coef(cox)))
module.gene.cox
writeMatrix(module.gene.cox,'results/06.module/module.gene.coef.txt',row = F)

##########关键基因单因素森林图
module.gene.sig.cox=data.frame(Genes=rownames(tcga.cox.fit[tcga.module.risk$module.gene,]),
                               round(tcga.cox.fit[tcga.module.risk$module.gene,],3))
pdf('results/06.module/Fig6b.pdf',height = 5,width = 6,onefile = F)
mg_forestplot_v2(module.gene.sig.cox[order(module.gene.sig.cox$HR),],xlog = T,colgap = 15,lineheight = 10,zero = 1,
                 boxsize = 0.3,lwd.zero=1,lwd.ci=1,lwd.xaxis=1
                 ,box_col='blue',summary_col="black",lines_col='blue',zero_col='grey'
                 ,xlab='Hazard Ratio',lty.ci = 7,graph.pos = 2)
dev.off()




fig6c=ggplotTimeROC(time = tcga.module.risk$result$time/365
             ,status = tcga.module.risk$result$status
             ,score = tcga.module.risk$result$riskscore
             ,mks = c(1,3,5))
tcga.module.risk$result$Risktype=ifelse(tcga.module.risk$result$riskscore>median(tcga.module.risk$result$riskscore),'High','Low')
tcga.risktype.cli=data.frame(tcga_cli,
                             Cluster=tcga.subtype[tcga_cli$Samples,]$Cluster,
                             Riskscore=tcga.module.risk$result[tcga_cli$Samples,]$riskscore,
                             Risktype=tcga.module.risk$result[tcga_cli$Samples,]$Risktype)


fig6d=ggsurvplot(fit=survfit( Surv(OS.time/365, OS) ~ Risktype,
                      data = tcga.risktype.cli),
         data=tcga.risktype.cli,
         conf.int = F,pval = T,fun = "pct",risk.table = T, size = 0.7,
         title='TCGA GI cancer',ggtheme=custom_theme(),
         linetype = c("solid", "dashed","strata")[1],
         palette = ggsci::pal_lancet()(8),
         #legend = c('top', 'bottom', 'left', 'right', 'none')[1],
         legend = c(0.8,0.75), # 指定图例位置
         legend.title = "",
         legend.labs = c("High","Low"))
fig6d=fig6d$plot
fig6d


chisq.test(table(tcga.risktype.cli$Status,tcga.risktype.cli$Risktype))$p.value
tcga.os.per=prop.table(table(tcga.risktype.cli$Status,tcga.risktype.cli$Risktype),margin=2)
tcga.os.per=reshape2::melt(tcga.os.per)
colnames(tcga.os.per)<-c("Status","Risktype","Percentage")
tcga.os.per$Percentage<-round(tcga.os.per$Percentage,digits=2)

fig6e=ggplot(tcga.os.per,aes(x=Risktype,y=Percentage,fill=Status))+
  geom_bar(position = "fill",stat="identity")+
  scale_fill_manual(values = c('red','blue'))+
  labs(x="",y = "Percentage",fill = "Status",title= 'Chi-Squared Test pvalue = 1.4e-07')+ 
  theme_bw()+coord_flip()
fig6e

fig6cde=mg_merge_plot(fig6c,fig6d,fig6e,ncol=3,labels = LETTERS[4:7])
savePDF('results/06.module/Fig6cde.pdf',fig6cde,height = 5,width = 20)

########07.GEO+泛癌验证##############
dir.create('results/07.Model verification')
# GSE66229[OK]##########
load('origin_datas/GEO/GSE66229_cli_exp.RData')
GSE66229_exp[1:5,1:5]
head(GSE66229_cli)


GSE66229.module.risk=get_riskscore(dat = t(GSE66229_exp[intersect(tcga.module.risk$module.gene,rownames(GSE66229_exp)),GSE66229_cli$Samples]),
                                   os = GSE66229_cli$OS,
                                   os.time = GSE66229_cli$OS.time,
                                   step = F,direction = c("both", "backward", "forward")[1])

fig7a=ggplotTimeROC(time = GSE66229.module.risk$result$time/365
              ,status = GSE66229.module.risk$result$status
              ,score = GSE66229.module.risk$result$riskscore
              ,mks = c(1,3,5))


GSE66229.module.risk$result$Risktype=ifelse(GSE66229.module.risk$result[GSE66229_cli$Samples,]$riskscore>median(GSE66229.module.risk$result[GSE66229_cli$Samples,]$riskscore),'High','Low')
GSE66229.risktype.cli=data.frame(GSE66229_cli,
                                 Riskscore=GSE66229.module.risk$result[GSE66229_cli$Samples,]$riskscore,
                                 Risktype=GSE66229.module.risk$result[GSE66229_cli$Samples,]$Risktype)


fig7b=ggsurvplot(fit=survfit( Surv(OS.time/365, OS) ~ Risktype,
                        data = GSE66229.risktype.cli),
           data=GSE66229.risktype.cli,
           conf.int = F,pval = T,fun = "pct",risk.table = T, size = 0.7,
           title='GSE66229',ggtheme=custom_theme(),
           linetype = c("solid", "dashed","strata")[1],
           palette = ggsci::pal_lancet()(8),
           #legend = c('top', 'bottom', 'left', 'right', 'none')[1],
           legend = c(0.8,0.85), # 指定图例位置
           legend.title = "",
           legend.labs = c("High","Low"))
fig7b=fig7b$plot
fig7b


########胰腺癌paad[ok]#######
paad_cli<-read.delim('origin_datas/TCGA/Merge_PAAD_clinical.txt',sep='\t',header = T)
colnames(paad_cli)[1:30]
paad_cli=data.frame(Samples=paste0(paad_cli$A0_Samples,'-01'),
                    Age=paad_cli$A17_Age,
                    Gender=paad_cli$A18_Sex,
                    T.stage=paad_cli$A3_T,
                    N.stage=paad_cli$A4_N,
                    M.stage=paad_cli$A5_M,
                    Stage=paad_cli$A6_Stage,#Grade=paad_cli$A7_Grade,
                    OS.time=paad_cli$A1_OS,Status=paad_cli$A2_Event)
rownames(paad_cli)=paad_cli$Samples
paad_cli=paad_cli %>% drop_na(OS.time)
paad_cli$OS.time
paad_cli=paad_cli[which(paad_cli$OS.time>0),]
table(paad_cli$Status)
paad_cli$OS=ifelse(paad_cli$Status=='Alive',0,1)


##########胰腺癌表达谱
paad_data<-read.delim('origin_datas/TCGA/PAAD_TPM.txt',sep='\t',header = T,row.names = 1,check.names = F)
paad_data[1:4,1:4]
table(substr(colnames(paad_data),14,15))


paad_tpm_log_T=log2(paad_data[,intersect(paad_cli$Samples,colnames(paad_data))]+1)
dim(paad_tpm_log_T)
paad_cli=paad_cli[intersect(paad_cli$Samples,colnames(paad_data)),]
dim(paad_cli)
#176  

paad.module.risk=get_riskscore(dat = t(paad_tpm_log_T[tcga.module.risk$module.gene,paad_cli$Samples]),
                               os = paad_cli$OS,
                               os.time = paad_cli$OS.time,
                               step = F,direction = c("both", "backward", "forward")[1])

fig7c=ggplotTimeROC(time = paad.module.risk$result$time/365
              ,status = paad.module.risk$result$status
              ,score = paad.module.risk$result$riskscore
              ,mks = c(1,3,5))


paad.module.risk$result$Risktype=ifelse(paad.module.risk$result$riskscore>median(paad.module.risk$result$riskscore),'High','Low')
paad.risktype.cli=data.frame(paad_cli,
                             #Cluster=paad.subtype[paad_cli$Samples,]$Cluster,
                             Riskscore=paad.module.risk$result[paad_cli$Samples,]$riskscore,
                             Risktype=paad.module.risk$result[paad_cli$Samples,]$Risktype)


fig7d=ggsurvplot(fit=survfit( Surv(OS.time/365, OS) ~ Risktype,
                        data = paad.risktype.cli),
           data=paad.risktype.cli,
           conf.int = F,pval = T,fun = "pct",risk.table = T, size = 0.7,
           title='TCGA PAAD',ggtheme=custom_theme(),
           linetype = c("solid", "dashed","strata")[1],
           palette = ggsci::pal_lancet()(8),
           #legend = c('top', 'bottom', 'left', 'right', 'none')[1],
           legend = c(0.8,0.75), # 指定图例位置
           legend.title = "",
           legend.labs = c("High","Low"))
fig7d=fig7d$plot
fig7d

####肝癌LIHC[ok]###########
lihc_cli<-read.delim('origin_datas/TCGA/Merge_LIHC_clinical.txt',sep='\t',header = T)
colnames(lihc_cli)[1:30]
lihc_cli=data.frame(Samples=paste0(lihc_cli$A0_Samples,'-01'),
                    Age=lihc_cli$A17_Age,
                    Gender=lihc_cli$A18_Sex,
                    T.stage=lihc_cli$A3_T,
                    N.stage=lihc_cli$A4_N,
                    M.stage=lihc_cli$A5_M,
                    Stage=lihc_cli$A6_Stage,#Grade=lihc_cli$A7_Grade,
                    OS.time=lihc_cli$A1_OS,Status=lihc_cli$A2_Event)
rownames(lihc_cli)=lihc_cli$Samples
lihc_cli=lihc_cli %>% drop_na(OS.time)
lihc_cli$OS.time
lihc_cli=lihc_cli[which(lihc_cli$OS.time>0),]
table(lihc_cli$Status)
lihc_cli$OS=ifelse(lihc_cli$Status=='Alive',0,1)



##########肝癌表达谱
lihc_data<-read.delim('origin_datas/TCGA/LIHC_TPM.txt',sep='\t',header = T,row.names = 1,check.names = F)
lihc_data[1:4,1:4]
table(substr(colnames(lihc_data),14,15))

lihc_tpm_log_T=log2(lihc_data[,intersect(lihc_cli$Samples,colnames(lihc_data))]+1)
dim(lihc_tpm_log_T)
lihc_cli=lihc_cli[intersect(lihc_cli$Samples,colnames(lihc_data)),]
dim(lihc_cli)
#365  11

lihc.module.risk=get_riskscore(dat = t(lihc_tpm_log_T[tcga.module.risk$module.gene,lihc_cli$Samples]),
                               os = lihc_cli$OS,
                               os.time = lihc_cli$OS.time,
                               step = F,direction = c("both", "backward", "forward")[1])

fig7e=ggplotTimeROC(time = lihc.module.risk$result$time/365
              ,status = lihc.module.risk$result$status
              ,score = lihc.module.risk$result$riskscore
              ,mks = c(1,3,5))

lihc.module.risk$result$Risktype=ifelse(lihc.module.risk$result$riskscore>median(lihc.module.risk$result$riskscore),'High','Low')
lihc.risktype.cli=data.frame(lihc_cli,
                             #Cluster=lihc.subtype[lihc_cli$Samples,]$Cluster,
                             Riskscore=lihc.module.risk$result[lihc_cli$Samples,]$riskscore,
                             Risktype=lihc.module.risk$result[lihc_cli$Samples,]$Risktype)

fig7f=ggsurvplot(fit=survfit( Surv(OS.time/365, OS) ~ Risktype,
                        data = lihc.risktype.cli),
           data=lihc.risktype.cli,
           conf.int = F,pval = T,fun = "pct",risk.table = T, size = 0.7,
           title='TCGA LIHC',ggtheme=custom_theme(),
           linetype = c("solid", "dashed","strata")[1],
           palette = ggsci::pal_lancet()(8),
           #legend = c('top', 'bottom', 'left', 'right', 'none')[1],
           legend = c(0.8,0.75), # 指定图例位置
           legend.title = "",
           legend.labs = c("High","Low"))
fig7f=fig7f$plot
fig7f


########胆管癌chol[ok]#######
chol_cli<-read.delim('origin_datas/TCGA/Merge_CHOL_clinical.txt',sep='\t',header = T)
colnames(chol_cli)[1:30]
chol_cli=data.frame(Samples=paste0(chol_cli$A0_Samples,'-01'),
                    Age=chol_cli$A17_Age,
                    Gender=chol_cli$A18_Sex,
                    T.stage=chol_cli$A3_T,
                    N.stage=chol_cli$A4_N,
                    M.stage=chol_cli$A5_M,
                    Stage=chol_cli$A6_Stage,#Grade=chol_cli$A7_Grade,
                    OS.time=chol_cli$A1_OS,Status=chol_cli$A2_Event)
rownames(chol_cli)=chol_cli$Samples
chol_cli=chol_cli %>% drop_na(OS.time)
chol_cli$OS.time
chol_cli=chol_cli[which(chol_cli$OS.time>0),]
table(chol_cli$Status)
chol_cli$OS=ifelse(chol_cli$Status=='Alive',0,1)

##########胆管癌表达谱
chol_data<-read.delim('origin_datas/TCGA/CHOL_TPM.txt',sep='\t',header = T,row.names = 1,check.names = F)
chol_data[1:4,1:4]
table(substr(colnames(chol_data),14,15))

chol_tpm_log_T=log2(chol_data[,intersect(chol_cli$Samples,colnames(chol_data))]+1)
dim(chol_tpm_log_T)
chol_cli=chol_cli[intersect(chol_cli$Samples,colnames(chol_data)),]
dim(chol_cli)
#36

chol.module.risk=get_riskscore(dat = t(chol_tpm_log_T[tcga.module.risk$module.gene,chol_cli$Samples]),
                               os = chol_cli$OS,
                               os.time = chol_cli$OS.time,
                               step = F,direction = c("both", "backward", "forward")[1])

fig7g=ggplotTimeROC(time = chol.module.risk$result$time/365
              ,status = chol.module.risk$result$status
              ,score = chol.module.risk$result$riskscore
              ,mks = c(1,3,5))

chol.module.risk$result$Risktype=ifelse(chol.module.risk$result$riskscore>median(chol.module.risk$result$riskscore),'High','Low')
chol.risktype.cli=data.frame(chol_cli,
                             #Cluster=chol.subtype[chol_cli$Samples,]$Cluster,
                             Riskscore=chol.module.risk$result[chol_cli$Samples,]$riskscore,
                             Risktype=chol.module.risk$result[chol_cli$Samples,]$Risktype)


fig7h=ggsurvplot(fit=survfit( Surv(OS.time/365, OS) ~ Risktype,
                        data = chol.risktype.cli),
           data=chol.risktype.cli,
           conf.int = F,pval = T,fun = "pct",risk.table = T, size = 0.7,
           title='TCGA CHOL',ggtheme=custom_theme(),
           linetype = c("solid", "dashed","strata")[1],
           palette = ggsci::pal_lancet()(8),
           #legend = c('top', 'bottom', 'left', 'right', 'none')[1],
           legend = c(0.8,0.75), # 指定图例位置
           legend.title = "",
           legend.labs = c("High","Low"))
fig7h=fig7h$plot
fig7h

fig7=mg_merge_plot(fig7a,fig7c,fig7e,fig7g,fig7b,fig7d,fig7f,fig7h,ncol=4,nrow=2,labels = LETTERS[1:8])
savePDF('results/07.Model verification/Fig7.pdf',fig7,height = 8,width = 16)

############08.RiskScore 与临床特征间的关系#########
dir.create('results/08.risktype.cli')
my_barplot=function(dat,ist=F,margin=T,xlb='',ylb='',lineCol=c('black','NA')[2],
                    lineW=0.5,legTitle='Group',showValue=F,xangle=0,
                    legend.position='right'){
  library(ggplot2)
  if(ist){
    dat=t(dat)
  }
  lbc=colnames(dat)
  lbr=row.names(dat)
  bk_dat=dat
  if(margin){
    dat=dat%*%diag(1/c(apply(t(dat), 1, sum)))
  }
  row.names(dat)=paste0('R',1:(nrow(dat)))
  colnames(dat)=paste0('C',1:(ncol(dat)))
  row.names(bk_dat)=paste0('R',1:(nrow(bk_dat)))
  colnames(bk_dat)=paste0('C',1:(ncol(bk_dat)))
  tp.dat=as.data.frame(cbind(bg=row.names(dat),dat))
  tp.dat[,1]=as.character(tp.dat[,1])
  for(i in 2:ncol(tp.dat)){
    tp.dat[,i]=as.numeric(as.character(tp.dat[,i]))
  }
  mt.df=reshape2::melt(tp.dat)
  colnames(mt.df)=c('bg','variable','value')
  
  ##画图
  pg=ggplot(mt.df, aes(x=variable, y=value, fill=bg))+geom_bar(stat = "identity", width=lineW, col=lineCol) 
  if(showValue){ ###是否显示比例
    pg=pg+geom_text(data=mt.df,aes(label=sprintf("%0.2f", round(value, digits = 2))),position=position_stack(vjust=0.5))
  }
  pg=pg+scale_x_discrete(breaks = paste0('C',1:(ncol(dat))),label = lbc)
  pg=pg+labs(x=xlb, y=ylb)+theme(legend.position =legend.position)
  pg=pg+scale_fill_discrete(breaks = paste0('R',1:nrow(dat)),label = lbr,name=legTitle,)#
  
  if(xangle>0){
    pg=pg+theme(axis.text.x = element_text(angle = xangle, hjust = 1))
  }
  if(xangle>0){ 
    pg=pg+theme(axis.text.x = element_text(angle = xangle, hjust = 0.5))
  }
  pval=round(chisq.test(bk_dat[,c(1,2)])$p.value,5)
  p=pg+ggtitle(paste0('p value = ',pval))
  return(p)
}

bar.cli1=my_barplot(table(tcga.risktype.cli$Gender,tcga.risktype.cli$Risktype))
bar.cli2=my_barplot(table(tcga.risktype.cli$Age1,tcga.risktype.cli$Risktype))
bar.cli3=my_barplot(table(tcga.risktype.cli$T.stage,tcga.risktype.cli$Risktype))
bar.cli4=my_barplot(table(tcga.risktype.cli$N.stage,tcga.risktype.cli$Risktype))
bar.cli5=my_barplot(table(tcga.risktype.cli$M.stage,tcga.risktype.cli$Risktype))
bar.cli6=my_barplot(table(tcga.risktype.cli$Stage,tcga.risktype.cli$Risktype))
# bar.cli7=my_barplot(table(tcga.risktype.cli$Status,tcga.risktype.cli$Risktype))



tcga_cox_datas=tcga.risktype.cli
colnames(tcga_cox_datas)
table(tcga_cox_datas$Age1)
table(tcga_cox_datas$T.stage)
tcga_cox_datas$T.stage[tcga_cox_datas$T.stage=='T1'|tcga_cox_datas$T.stage=='T2']<-'T1+T2'
tcga_cox_datas$T.stage[tcga_cox_datas$T.stage=='T3'|tcga_cox_datas$T.stage=='T4']<-'T3+T4'

table(tcga_cox_datas$N.stage)
tcga_cox_datas$N.stage[tcga_cox_datas$N.stage=='N1'|tcga_cox_datas$N.stage=='N2'|tcga_cox_datas$N.stage=='N3']<-'N1+N2+N3'

table(tcga_cox_datas$M.stage)

table(tcga_cox_datas$Stage)
tcga_cox_datas$Stage[tcga_cox_datas$Stage=='I'|tcga_cox_datas$Stage=='II']<-'I+II'
tcga_cox_datas$Stage[tcga_cox_datas$Stage=='III'|tcga_cox_datas$Stage=='IV']<-'III+IV'

p1=my_violin(dat = tcga.risktype.cli$Riskscore,
          group = tcga.risktype.cli$Gender,title = 'Gender',
          legend.position = 'none',ylab = 'Riskscore',
          group_cols = qualitative_hcl(4, palette = "Dark 3"),
          test_method =  c('kruskal.test','wilcox.test')[2])

p2=my_violin(dat = tcga.risktype.cli$Riskscore,
          group = tcga.risktype.cli$Age1,title = 'Age',
          legend.position = 'none',ylab = 'Riskscore',
          group_cols = qualitative_hcl(4, palette = "Dark 3"),
          test_method =  c('kruskal.test','wilcox.test')[2])

p3=my_violin(dat = tcga.risktype.cli$Riskscore,
          group = tcga.risktype.cli$T.stage,title = 'T.stage',
          legend.position = 'none',ylab = 'Riskscore',
          group_cols = qualitative_hcl(4, palette = "Dark 3"),
          test_method =  c('kruskal.test','wilcox.test')[1])

p4=my_violin(dat = tcga.risktype.cli$Riskscore,
          group = tcga.risktype.cli$N.stage,title = 'N.stage',
          legend.position = 'none',ylab = 'Riskscore',
          group_cols = qualitative_hcl(4, palette = "Dark 3"),
          test_method =  c('kruskal.test','wilcox.test')[1])

p5=my_violin(dat = tcga.risktype.cli$Riskscore,
          group = tcga.risktype.cli$M.stage,title = 'M.stage',
          legend.position = 'none',ylab = 'Riskscore',
          group_cols = qualitative_hcl(4, palette = "Dark 3"),
          test_method =  c('kruskal.test','wilcox.test')[2])


p6=my_violin(dat = tcga.risktype.cli$Riskscore,
          group = tcga.risktype.cli$Stage,title = 'Stage',
          legend.position = 'none',ylab = 'Riskscore',
          group_cols = qualitative_hcl(4, palette = "Dark 3"),
          test_method =  c('kruskal.test','wilcox.test')[1])

p7=my_violin(dat = tcga.risktype.cli$Riskscore,
          group = tcga.risktype.cli$Cluster,title = 'Cluster',
          legend.position = 'none',ylab = 'Riskscore',
          group_cols = qualitative_hcl(4, palette = "Dark 3"),
          test_method =  c('kruskal.test','wilcox.test')[1])

pdf('results/08.risktype.cli/Fig8a.pdf',height = 8,width = 18)
mg_merge_plot(bar.cli1,bar.cli2,bar.cli3,bar.cli4,bar.cli5,bar.cli6,bar.cli7,
                    p1,p2,p3,p4,p5,p6,p7,ncol=6,nrow=2,labels = 'A')
dev.off()

######风险模型与临床特征在患者预后诊断比较####
#Age
tcga_cox_datas=crbind2DataFrame(tcga_cox_datas)
Age_sig_cox <- summary(coxph(formula=Surv(OS.time, OS)~Age,
                             data=tcga_cox_datas))
Age_sig_cox_dat <- data.frame(Names=rownames(Age_sig_cox[[8]]),
                              HR = round(Age_sig_cox[[7]][,2],3),
                              lower.95 = round(Age_sig_cox[[8]][,3],3),
                              upper.95 = round(Age_sig_cox[[8]][,4],3),
                              p.value=round(Age_sig_cox[[7]][,5],3))
Age_sig_cox_dat

#Gender
Gender_sig_cox <- summary(coxph(formula=Surv(OS.time, OS)~Gender,
                                data=tcga_cox_datas))
Gender_sig_cox_dat <- data.frame(Names=rownames(Gender_sig_cox[[8]]),
                                 HR = round(Gender_sig_cox[[7]][,2],3),
                                 lower.95 = round(Gender_sig_cox[[8]][,3],3),
                                 upper.95 = round(Gender_sig_cox[[8]][,4],3),
                                 p.value=round(Gender_sig_cox[[7]][,5],3))
Gender_sig_cox_dat

#T.stage
T.stage_sig_cox <- summary(coxph(formula=Surv(OS.time, OS)~T.stage,
                                 data=tcga_cox_datas))
T.stage_sig_cox_dat <- data.frame(Names=rownames(T.stage_sig_cox[[8]]),
                                  HR = round(T.stage_sig_cox[[7]][,2],3),
                                  lower.95 = round(T.stage_sig_cox[[8]][,3],3),
                                  upper.95 = round(T.stage_sig_cox[[8]][,4],3),
                                  p.value=round(T.stage_sig_cox[[7]][,5],3))
T.stage_sig_cox_dat

#N.stage
N.stage_sig_cox <- summary(coxph(formula=Surv(OS.time, OS)~N.stage,
                                 data=tcga_cox_datas))
N.stage_sig_cox_dat <- data.frame(Names=rownames(N.stage_sig_cox[[8]]),
                                  HR = round(N.stage_sig_cox[[7]][,2],3),
                                  lower.95 = round(N.stage_sig_cox[[8]][,3],3),
                                  upper.95 = round(N.stage_sig_cox[[8]][,4],3),
                                  p.value=round(N.stage_sig_cox[[7]][,5],3))
N.stage_sig_cox_dat

#M.stage
M.stage_sig_cox <- summary(coxph(formula=Surv(OS.time, OS)~M.stage,
                                 data=tcga_cox_datas))
M.stage_sig_cox_dat <- data.frame(Names=rownames(M.stage_sig_cox[[8]]),
                                  HR = round(M.stage_sig_cox[[7]][,2],3),
                                  lower.95 = round(M.stage_sig_cox[[8]][,3],3),
                                  upper.95 = round(M.stage_sig_cox[[8]][,4],3),
                                  p.value=round(M.stage_sig_cox[[7]][,5],3))
M.stage_sig_cox_dat

#Stage
Stage_sig_cox <- summary(coxph(formula=Surv(OS.time, OS)~Stage,
                               data=tcga_cox_datas))
Stage_sig_cox_dat <- data.frame(Names=rownames(Stage_sig_cox[[8]]),
                                HR = round(Stage_sig_cox[[7]][,2],3),
                                lower.95 = round(Stage_sig_cox[[8]][,3],3),
                                upper.95 = round(Stage_sig_cox[[8]][,4],3),
                                p.value=round(Stage_sig_cox[[7]][,5],3))
Stage_sig_cox_dat


#riskscore
riskscore_sig_cox<-summary(coxph(formula=Surv(OS.time, OS)~Riskscore,
                                 data=tcga_cox_datas))
riskscore_sig_cox_dat <- data.frame(Names=rownames(riskscore_sig_cox[[8]]),
                                    HR = round(riskscore_sig_cox[[7]][,2],3),
                                    lower.95 = round(riskscore_sig_cox[[8]][,3],3),
                                    upper.95 = round(riskscore_sig_cox[[8]][,4],3),
                                    p.value=round(riskscore_sig_cox[[7]][,5],3))
riskscore_sig_cox_dat

sig_cox_dat <- rbind(Age_sig_cox_dat,
                     Gender_sig_cox_dat,
                     T.stage_sig_cox_dat,
                     N.stage_sig_cox_dat,
                     M.stage_sig_cox_dat,
                     Stage_sig_cox_dat,
                     riskscore_sig_cox_dat)
data.sig <- data.frame(Names=sig_cox_dat$Names,
                       p.value=sig_cox_dat$p.value,
                       sig_cox_dat$HR,
                       sig_cox_dat$lower.95,
                       sig_cox_dat$upper.95)
data.sig <- crbind2DataFrame(data.sig)
data.sig$Names
rownames(data.sig) <- c("Age",
                        "Gender",
                        "T.stage",
                        "N.stage",
                        "M.stage",
                        "Stage",
                        "RiskScore")
data.sig$Names <- rownames(data.sig)
data.sig
data.sig$p.value=ifelse(data.sig$p.value<0.001,'<0.001',data.sig$p.value)
pdf('results/08.risktype.cli/Fig8b.pdf',height = 6,width = 8,onefile = F)
mg_forestplot_v2(data.sig,xlog = T,colgap = 15,lineheight = 10,zero = 1,
                 boxsize = 0.3,lwd.zero=1,lwd.ci=1,lwd.xaxis=1
                 ,box_col='blue',summary_col="black",lines_col='blue',zero_col='grey'
                 ,xlab='Hazard Ratio',lty.ci = 10,graph.pos = 3)
dev.off()


#多因素
#T.stage +  N.stage +M.stage
muti_sig_cox <- summary(coxph(formula=Surv(OS.time, OS)~Age+Stage+T.stage +  N.stage +M.stage+Riskscore, 
                              data=tcga_cox_datas))
muti_cox_dat <- data.frame(Names=rownames(muti_sig_cox[[8]]),
                           HR = round(muti_sig_cox[[7]][,2],3),
                           lower.95 = round(muti_sig_cox[[8]][,3],3),
                           upper.95 = round(muti_sig_cox[[8]][,4],3),
                           p.value=round(muti_sig_cox[[7]][,5],3))
data.muti <- data.frame(Names=muti_cox_dat$Names,
                        p.value=muti_cox_dat$p.value,
                        muti_cox_dat$HR,
                        muti_cox_dat$lower.95,
                        muti_cox_dat$upper.95)
data.muti <- crbind2DataFrame(data.muti)
data.muti
rownames(data.muti) <- c("Age",
                         #"Gender",
                         "Stage",
                         "T.stage",
                         "N.stage",
                         "M.stage",
                         "RiskScore")
data.muti$Names <- rownames(data.muti)
data.muti
data.muti$p.value=ifelse(data.muti$p.value<0.001,'<0.001',data.muti$p.value)
pdf('results/08.risktype.cli/Fig8c.pdf',height = 6,width = 8,onefile = F)
mg_forestplot_v2(data.muti,xlog = T,colgap = 15,lineheight = 10,zero = 1,
                 boxsize = 0.3,lwd.zero=1,lwd.ci=1,lwd.xaxis=1
                 ,box_col='blue',summary_col="black",lines_col='blue',zero_col='grey'
                 ,xlab='Hazard Ratio',lty.ci = 10,graph.pos = 3)
dev.off()


pdf('results/08.risktype.cli/Fig8d.pdf', width = 12, height = 10)
nom.plot=mg_nomogram(data.frame(RiskScore=tcga_cox_datas$Riskscore,
                                Age=tcga_cox_datas$Age,
                                N.stage=tcga_cox_datas$N.stage,
                                M.stage=tcga_cox_datas$M.stage),
                     os = tcga_cox_datas$OS.time,
                     status = tcga_cox_datas$OS,
                     mks = c(1,3,5))
dev.off()
mg_nomogram_buti(nom.plot$Mod,cut.time = c(1*365,3*365,5*365))


######09.分组间免疫特征#########
dir.create('results/09.risktype.immu')
pdf('results/09.risktype.immu/Fig9a.pdf',height = 6,width = 15,onefile = F)
my_mutiboxplot(dat =tcga.ciber[tcga.risktype.cli$Samples,1:22],
               group = tcga.risktype.cli$Risktype,
               group_cols = qualitative_hcl(2, palette = "Dynamic"),
               test_method ='wilcox.test',bw = T,
               ylab = 'Fraction',fill = 'Risk type')
dev.off()

mg_PlotMutiBoxplot(tcga.est[tcga.risktype.cli$Samples,]
                   , group = tcga.risktype.cli$Risktype
                   , legend.pos = 'top'
                   #, group_cols = cluster.color
                   , test_method =  c('kruskal.test','wilcox.test')[2]
                   , add = 'boxplot'
                   , ylab = 'Score')
est_cor_RS=cbind.data.frame(Riskscore=tcga.risktype.cli$Riskscore,
                            tcga.est[tcga.risktype.cli$Samples,])
pdf('results/09.risktype.immu/Fig9b.pdf',height = 6,width = 6,onefile = F)
corrplot(cor(est_cor_RS),
         tl.col = 'black',
         # cl.align.text = 'l',
         tl.srt=90, diag = F,
         col=colorRampPalette(c('blue', 'white','red'))(100),
         method = c("circle", "square", "ellipse", "number", "shade", "color", "pie")[7],
         order = c("original", "AOE", "FPC", "hclust", "alphabet")[1],
         type="upper",)
dev.off()



mg_PlotMutiBoxplot(tcga.immu.ssgsea[tcga.risktype.cli$Samples,]
                   , group = tcga.risktype.cli$Risktype
                   , legend.pos = 'top'
                   #, group_cols = cluster.color
                   , test_method =  c('kruskal.test','wilcox.test')[2]
                   , add = 'boxplot'
                   , ylab = 'Score')
library(tidyverse)
library(ggcor)
library(vegan)
cr=psych::corr.test(x=tcga.risktype.cli$Riskscore,
                    y=tcga.immu.ssgsea[tcga.risktype.cli$Samples,]
                    ,method = 'spearman')
df=t(rbind(cr$r,cr$p))
colnames(df)=c('r','p.value')
df=data.frame(Riskscore='Riskscore',tcga.immu.ssgsea=rownames(df),df)
df
df <- df %>%
  mutate(lty = cut(r, breaks = c(-1, 0, 1),
                   labels = c("r <= 0", "r > 0")),
         col = cut(p.value, breaks = c(0, 0.01, 0.05, 1),
                   labels = c("< 0.01", "< 0.05", ">= 0.05"),
                   right = FALSE, include.lowest = TRUE))
head(df)

corrmat.color=colorRampPalette(c('blue', 'white','orange'))(100)
quickcor(tcga.immu.ssgsea[tcga.risktype.cli$Samples,], cor.test = TRUE,type = "lower") + #upper
  geom_square(data = get_data(p.value < 0.05, type = "lower")) + 
  anno_link(df, mapping = aes(colour = col,
                              size = abs(r),
                              linetype = lty),diag.label = TRUE) +
  scale_fill_gradient2n(colours = corrmat.color) +
  remove_x_axis()+
  scale_size_area(max_size = 2) +
  scale_linetype_manual(values = c("dotted", "solid")) +
  guides(
    fill = guide_colourbar(title = "corr", order = 1),
    colour = guide_legend(title = "spearman's p", order = 2),
    size = guide_legend(title = "spearman's r", order = 3),
    linetype = "none")


load('results/tcga.tme.ssgsea.RData')
tcga.tme.ssgsea=t(tcga.tme.ssgsea)
mg_PlotMutiBoxplot(tcga.tme.ssgsea[tcga.risktype.cli$Samples,]
                   , group = tcga.risktype.cli$Risktype
                   , legend.pos = 'top'
                   #, group_cols = cluster.color
                   , test_method =  c('kruskal.test','wilcox.test')[2]
                   , add = 'boxplot'
                   , ylab = 'Score')


cr=psych::corr.test(x=tcga.risktype.cli$Riskscore,
                    y=tcga.tme.ssgsea[tcga.risktype.cli$Samples,]
                    ,method = 'spearman')
df=t(rbind(cr$r,cr$p))
colnames(df)=c('r','p.value')
df=data.frame(Riskscore='Riskscore',tcga.tme.ssgsea=rownames(df),df)
df
df <- df %>%
  mutate(lty = cut(r, breaks = c(-1, 0, 1),
                   labels = c("r <= 0", "r > 0")),
         col = cut(p.value, breaks = c(0, 0.01, 0.05, 1),
                   labels = c("< 0.01", "< 0.05", ">= 0.05"),
                   right = FALSE, include.lowest = TRUE))
head(df)

corrmat.color=colorRampPalette(c('blue', 'white','orange'))(100)
quickcor(tcga.tme.ssgsea[tcga.risktype.cli$Samples,], cor.test = TRUE,type = "upper") + #upper
  geom_square(data = get_data(p.value < 0.05, type = "upper")) + 
  anno_link(df, mapping = aes(colour = col,
                              size = abs(r),
                              linetype = lty),diag.label = TRUE) +
  scale_fill_gradient2n(colours = corrmat.color) +
  remove_x_axis()+
  scale_size_area(max_size = 2) +
  scale_linetype_manual(values = c("dotted", "solid")) +
  guides(
    fill = guide_colourbar(title = "corr", order = 1),
    colour = guide_legend(title = "spearman's p", order = 2),
    size = guide_legend(title = "spearman's r", order = 3),
    linetype = "none")




################免疫治疗
pdf('results/09.risktype.immu/Fig9e.pdf',height = 9,width = 14)
my_mutiboxplot_seg(dat = tcga.icgs[tcga.risktype.cli$Samples,icg.genes],
                   group = tcga.risktype.cli$Risktype,
                   group_cols = qualitative_hcl(2, palette = "Dynamic"),
                   test_method = 'wilcox.test',fill = 'Risktype',
                   ylab = 'log2(TPM+1)',legend.position = 'none',
                   ncol = 4,nrow=2,xangle = 0,xhjust = 0.5,
                   xsize = 12,ysize = 12)
dev.off()


####与TIDE评分相关性
TIDE_cor_RS=cbind.data.frame(tcga_tide_res[tcga.risktype.cli$Samples,tide_sel],
                            Riskscore=tcga.risktype.cli$Riskscore)
ggcorplot <- function(plot_df,a,b,method="spearman"){
  corr_eqn <- function(x,y,digits=3) {
    test <- cor.test(x,y,method=method,exact=FALSE)
    paste(#paste0("Method = ",method),
      paste0("r = ",round(test$estimate,digits)),
      paste0("p.value= ",round(test$p.value,digits)),
      sep = ", ")
  }
  plot_df <- plot_df[,c(a,b)]
  names(plot_df) <- c("var1","var2")
  require(ggplot2)
  p=ggplot(plot_df,aes(var1,var2))+
    geom_point(col="black")+
    geom_smooth(method=lm, se=T,na.rm=T, fullrange=T,size=2,col="red")+
    geom_rug(col="#006fbc")+
    theme_minimal()+
    xlab(a)+
    ylab(b)+
    ## 依靠函数来生成title
    labs(title = corr_eqn(plot_df$var1,plot_df$var2),subtitle =paste0("Method = ",method) )+
    theme(plot.title = element_text(hjust = 0.5),
          plot.margin = margin(1, 1, 1, 1, "cm"))
  return(p)
}

fig9f=list()
for (i in 1:length(colnames(TIDE_cor_RS))) {
  print(i)
  #inputgene = yourgene
  #i=3
  fig9f[[i]] = ggcorplot(plot_df = TIDE_cor_RS,a = 'Riskscore', b = colnames(TIDE_cor_RS)[i], method = 'spearman')
}
fig9f=mg_merge_plot(fig9f[1:5],ncol = 1,nrow = 5)
fig9f
savePDF('results/09.risktype.immu/Fig9f.pdf',fig9f,height = 18,width = 5)



#######10.风险分组间通路/突变差异###########
dir.create('results/10.risktype.pathway')
###########通路活性差异
library(progeny)
tcga.pathway.activ=progeny(as.matrix(tcga_tpm_log_T),scale = T)
dim(tcga.pathway.activ)
range(tcga.pathway.activ)


pathway_cor_RS=cbind.data.frame(Riskscore=tcga.risktype.cli$Riskscore,
                                tcga.pathway.activ[tcga.risktype.cli$Samples,])
cor_res <- Hmisc::rcorr(as.matrix(pathway_cor_RS),type = 'spearman')
cor_res$P[is.na(cor_res$P)] <- 0

pdf('results/10.risktype.pathway/Fig10a.pdf',height = 7,width = 7)
corrplot(as.matrix(cor_res$r),
         p.mat = as.matrix(cor_res$P),
         mar = c(0,0,1,1),diag = F,
         col=diverging_hcl(100,palette = 'Green-Orange'),
         tl.srt = 90,tl.cex = 1,tl.col = 'black',tl.offset = 0.5,
         cl.pos = c("b","r","n")[1],cl.align.text = 'l',cl.length = 5,
         cl.ratio = 0.1,cl.cex = 0.8,
         addgrid.col = 'white',
         method = c("circle", "square", "ellipse", "number", "shade", "color", "pie")[6],
         insig = 'label_sig',
         sig.level=c(0.001,0.01,0.05),
         pch.cex=1,is.corr=T,xpd=T)
dev.off()


############fgsea
tcga.geneList=getGeneFC(gene.exp=tcga_tpm_log_T,group=tcga.risktype.cli$Risktype
                        ,ulab='High',dlab = 'Low')

gmt2list <- function(annofile){
  
  if (!file.exists(annofile)) {
    stop("There is no such a gmt file!")
  }
  
  if (tools::file_ext(annofile) == "xz") {
    annofile <- xzfile(annofile)
    x <- scan(annofile, what="", sep="\n", quiet=TRUE)
    close(annofile)
  } else if (tools::file_ext(annofile) == "gmt") {
    x <- scan(annofile, what="", sep="\n", quiet=TRUE)
  } else {
    stop ("Only gmt and gmt.xz are accepted for gmt2list")
  }
  
  y <- strsplit(x, "\t")
  names(y) <- sapply(y, `[[`, 1)
  
  annoList <- lapply(y, `[`, c(-1,-2))
  return(annoList)
}

library(fgsea)
fgseaRes <- fgsea(pathways = kegmt,
                  stats =tcga.geneList ,
                  minSize=10,
                  maxSize=500,
                  nperm=1000)
head(fgseaRes)

nrow(fgseaRes[padj<0.05&ES > 0])
#20
nrow(fgseaRes[padj<0.05&ES < 0])
#9

topPathwaysUp <- fgseaRes[ES > 0][head(order(padj), n=10), pathway]
topPathwaysDown <- fgseaRes[ES < 0][head(order(padj), n=9), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))

library(forcats)
library(ggplot2)
library(ggstance)
fgseaRes$sign<-ifelse(fgseaRes$ES>0,"Activated","Suppressed") ## 定义分组
fgseaRes$Description<- gsub('HALLMARK_','',fgseaRes$pathway) ##删除KEGG_的前缀
fgseaRes$Description<- gsub('_',' ',fgseaRes$Description) ## 将下划线替换为空

pdf('results/10.risktype.pathway/Fig10b.pdf',height = 7,width = 10)
fgseaRes %>% 
  group_by(sign) %>%  ## 很重要，先按照sign分组
  arrange(pval) %>% ## 按照p值排序
  slice(1:10) %>%  ## 设置各组要显示的数目，这个可以自己定义
  ## 开始ggplot2 作图
  ggplot( aes(NES, fct_reorder(Description, NES),fill=pval)) + 
  geom_col(orientation='y') + 
  scale_fill_continuous(low="#f87669", high="#2874C5", guide=guide_colorbar(reverse=TRUE)) +   
  ylab(NULL)+  ggtitle("GSEA HALLMARK")+ ##设置标题
  theme_bw(base_size = 12)+  ## 设定主题
  theme(plot.title = element_text(hjust = 0.5)) ##标题居中
dev.off()

head(fgseaRes)
class(fgseaRes)
aa=data.frame(fgseaRes)
write.xlsx(aa,'results/10.risktype.pathway/risktype_GESA.xlsx',overwrite = T)



#########体细胞突变#########
merge_mafs_use = function(mafs, verbose = TRUE, ...){
  
  if(all(unlist(lapply(mafs, is, "MAF")))){
    
    if(verbose){
      cat(paste0("Merging ", length(mafs) ," MAF objects\n"))
    }
    mafs_dat = lapply(mafs, function(m) {
      data.table::rbindlist(l = list(m@data, m@maf.silent), use.names = TRUE, fill = TRUE)
    })
    mafs_dat = data.table::rbindlist(l = mafs_dat, fill = TRUE, use.names = TRUE)
    
    mafs_clin = lapply(mafs, function(m) {
      m@clinical.data
    })
    mafs_clin = data.table::rbindlist(l = mafs_clin, fill = TRUE, use.names = TRUE)
    maf = read.maf(maf = mafs_dat, clinicalData = mafs_clin, verbose = verbose, ...)
  }else if(all(unlist(lapply(mafs, is, "data.frame")))){
    if(verbose){
      cat(paste0("Merging ", length(mafs) ," MAF data.frames\n"))
    }
    mafs_dat = data.table::rbindlist(l = mafs, fill = TRUE, use.names = TRUE)
    maf = read.maf(maf = mafs_dat, verbose = verbose, ...)
  }else{
    if(verbose){
      cat(paste0("Merging ", length(mafs) ," MAF files\n"))
    }
    maf = lapply(mafs, function(x) {
      x = data.table::fread(x, stringsAsFactors = FALSE, fill = TRUE, showProgress = TRUE, header = TRUE, skip = "Hugo_Symbol")
      
      required.fields = c(
        'Hugo_Symbol',
        'Chromosome',
        'Start_Position',
        'End_Position',
        'Reference_Allele',
        'Tumor_Seq_Allele2',
        'Variant_Classification',
        'Variant_Type',
        'Tumor_Sample_Barcode'
      )
      
      #Change column names to standard names; i.e, camel case
      for (i in 1:length(required.fields)) {
        colId = suppressWarnings(grep(
          pattern = paste("^", required.fields[i], "$", sep = ""),
          x = colnames(x),
          ignore.case = TRUE
        ))
        if (length(colId) > 0) {
          colnames(x)[colId] = required.fields[i]
        }
      }
      x
    })
    #names(maf) = gsub(pattern = "\\.maf$", replacement = "", x = basename(path = unlist(mafs)), ignore.case = TRUE)
    #names(maf) = basename(mafs)
    maf = data.table::rbindlist(l = maf, fill = TRUE, idcol = "Source_MAF", use.names = TRUE)
    
    maf = read.maf(maf = maf, verbose = verbose)
  }
  
  maf
}

# stad.maf=getTCGAMAFByCode('STAD')
# coad.maf=getTCGAMAFByCode('COAD')
# read.maf=getTCGAMAFByCode('READ')
# tcga.maf=merge_mafs_use(mafs = list(stad.maf,
#                                     coad.maf,
#                                     read.maf),
#                         verbose = TRUE)
# save(tcga.maf,file='results/tcga.maf.RData')
load('results/tcga.maf.RData')

tcga.risktype.use=tcga.risktype.cli[,c('Samples','Risktype')]
table(tcga.risktype.use$Risktype)
colnames(tcga.risktype.use)[1]='Tumor_Sample_Barcode'
tcga.risktype.use$Tumor_Sample_Barcode=substr(tcga.risktype.use$Tumor_Sample_Barcode,1,12)
tcga.risktype.use.high=tcga.risktype.use[which(tcga.risktype.use$Risktype=='High'),]
tcga.risktype.use.low=tcga.risktype.use[which(tcga.risktype.use$Risktype=='Low'),]

write.table(tcga.risktype.use.high,file='results/tcga.risktype.use.high.txt')
write.table(tcga.risktype.use.low,file='results/tcga.risktype.use.low.txt')

tcga.maf.high=subsetMaf(tcga.maf,tsb=intersect(tcga.maf@data$Tumor_Sample_Barcode,tcga.risktype.use.high$Tumor_Sample_Barcode))
tcga.maf.high<-read.maf(tcga.maf.high@data,isTCGA=T,clinicalData = 'results/tcga.risktype.use.high.txt')
tcga.maf.high@clinical.data

tcga.maf.low=subsetMaf(tcga.maf,tsb=intersect(tcga.maf@data$Tumor_Sample_Barcode,tcga.risktype.use.low$Tumor_Sample_Barcode))
tcga.maf.low<-read.maf(tcga.maf.low@data,isTCGA=T,clinicalData = 'results/tcga.risktype.use.low.txt')
tcga.maf.low@clinical.data


tcga.mut.dat <- tcga.maf
tcga.mut.dat <- as.data.frame(tcga.mut.dat@data)
tcga.mut.dat <- tcga.mut.dat[, c("Hugo_Symbol", "Tumor_Sample_Barcode", "Variant_Classification")]

tcga.mut.dat$Variant_Classification <- 1
tcga.mut.dat <- reshape2::dcast(data = tcga.mut.dat, Hugo_Symbol ~ Tumor_Sample_Barcode)
class(tcga.mut.dat)
rownames(tcga.mut.dat) <- tcga.mut.dat$Hugo_Symbol
tcga.mut.dat <- tcga.mut.dat[, -1]

colnames(tcga.mut.dat) <- paste0(colnames(tcga.mut.dat), '-01')
mut.samples <- intersect(colnames(tcga.mut.dat), tcga.risktype.cli$Samples)


tcga.mut.dat <- tcga.mut.dat[, mut.samples]
rownames(tcga.risktype.cli)=tcga.risktype.cli$Samples
tcga_mut_cli <- tcga.risktype.cli[mut.samples, ]

tcga.mut.dat.freq <- as.data.frame(rowSums(tcga.mut.dat))
colnames(tcga.mut.dat.freq) <- 'Freq'
tcga.mut.dat.freq$Genes <- rownames(tcga.mut.dat.freq)
library(dplyr)
str(tcga.mut.dat.freq)
head(tcga.mut.dat.freq)
tcga.mut.dat.freq <- dplyr::arrange(tcga.mut.dat.freq, desc(Freq))
head(tcga.mut.dat.freq)
dim(tcga.mut.dat.freq)

mut.genes <- rownames(tcga.mut.dat.freq)[tcga.mut.dat.freq$Freq > 3]
length(mut.genes)

tcga.mut.dat <- ifelse(tcga.mut.dat > 0, 'Mutant', 'WildType')
dim(tcga.mut.dat)

mut.res <- data.frame(High = NA,Low = NA)
mut.p <- c()

for (ge in mut.genes) {
  #  print(ge)
  tmp <- table(tcga.mut.dat[ge, ], tcga_mut_cli$Risktype)
  pvalue <- fisher.test(tmp)
  mut.p <- c(mut.p, pvalue$p.value)
  mut.res <- rbind(mut.res, tmp[1, ])
}
mut.res <- na.omit(mut.res)
rownames(mut.res) <- mut.genes
class(mut.res)
mut.res$P.value <- mut.p

table(mut.res$P.value < 0.05)
# FALSE  TRUE
#15483   812
mut.res.filtered <- mut.res[which(mut.res$P.value < 0.05), ]
mut.res.filtered
dim(mut.res.filtered)


pdf('results/10.risktype.pathway/Fig10c.pdf',height = 6,width = 16)
coOncoplot(m1=tcga.maf.high, m2=tcga.maf.low, m1Name="High", m2Name="low",genes =  rownames(mut.res.filtered)[1:15])
dev.off()

save.image(file='GI_cancer.RData')

