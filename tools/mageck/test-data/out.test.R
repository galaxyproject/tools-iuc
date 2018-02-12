pdf(file='output.pdf',width=4.5,height=4.5);
gstable=read.table('output.gene_summary.txt',header=T)
# 
#
# parameters
# Do not modify the variables beginning with "__"

# gstablename='__GENE_SUMMARY_FILE__'
startindex=3
# outputfile='__OUTPUT_FILE__'
targetgenelist=c("ACIN1","ACTR8","AHCY","ACLY","AATF","AGBL5","AHCTF1","ABT1","ADIRF","ABCF1")
# samplelabel=sub('.\\w+.\\w+$','',colnames(gstable)[startindex]);
samplelabel='HL60.final,KBM7.final_vs_HL60.initial,KBM7.initial neg.'


# You need to write some codes in front of this code:
# gstable=read.table(gstablename,header=T)
# pdf(file=outputfile,width=6,height=6)


# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")

######
# function definition

plotrankedvalues<-function(val, tglist, ...){
  
  plot(val,log='y',ylim=c(max(val),min(val)),type='l',lwd=2, ...)
  if(length(tglist)>0){
    for(i in 1:length(tglist)){
      targetgene=tglist[i];
      tx=which(names(val)==targetgene);ty=val[targetgene];
      points(tx,ty,col=colors[(i %% length(colors)) ],cex=2,pch=20)
      # text(tx+50,ty,targetgene,col=colors[i])
    }
    legend('topright',tglist,pch=20,pt.cex = 2,cex=1,col=colors)
  }
}



plotrandvalues<-function(val,targetgenelist, ...){
  # choose the one with the best distance distribution
  
  mindiffvalue=0;
  randval=val;
  for(i in 1:20){
    randval0=sample(val)
    vindex=sort(which(names(randval0) %in% targetgenelist))
    if(max(vindex)>0.9*length(val)){
      # print('pass...')
      next;
    }
    mindiffind=min(diff(vindex));
    if (mindiffind > mindiffvalue){
      mindiffvalue=mindiffind;
      randval=randval0;
      # print(paste('Diff: ',mindiffvalue))
    }
  }
  plot(randval,log='y',ylim=c(max(randval),min(randval)),pch=20,col='grey', ...)
  
  if(length(targetgenelist)>0){
    for(i in 1:length(targetgenelist)){
      targetgene=targetgenelist[i];
      tx=which(names(randval)==targetgene);ty=randval[targetgene];
      points(tx,ty,col=colors[(i %% length(colors)) ],cex=2,pch=20)
      text(tx+50,ty,targetgene,col=colors[i])
    }
  }
  
}




# set.seed(1235)



pvec=gstable[,startindex]
names(pvec)=gstable[,'id']
pvec=sort(pvec);

plotrankedvalues(pvec,targetgenelist,xlab='Genes',ylab='RRA score',main=paste('Distribution of RRA scores in \\n',samplelabel))

# plotrandvalues(pvec,targetgenelist,xlab='Genes',ylab='RRA score',main=paste('Distribution of RRA scores in \\n',samplelabel))


pvec=gstable[,startindex+1]
names(pvec)=gstable[,'id']
pvec=sort(pvec);

plotrankedvalues(pvec,targetgenelist,xlab='Genes',ylab='p value',main=paste('Distribution of p values in \\n',samplelabel))

# plotrandvalues(pvec,targetgenelist,xlab='Genes',ylab='p value',main=paste('Distribution of p values in \\n',samplelabel))



# you need to write after this code:
# dev.off()






# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(561.4907165816957,824.0396348113272,428.37415340969943,579.047491896501),c(3424.79939695118,3818.2871009576584,1992.3498917052,690.0506672205338),c(846.6456878299913,985.6508562937211,335.0024675413113,415.97581680707134),c(2432.636481525409,2122.257249136931,1067.465489792653,155.6333179800872),c(1308.1851773762019,2186.1913587343615,1482.5909580453515,997.3120339679854),c(405.68439208520414,268.16807081144486,170.34023773287015,109.85881269182627),c(640.8637498157573,559.4234589775174,711.6436598617687,632.2603542941043),c(946.5969148654764,470.6260845366416,663.0651476194316,457.74505288260946),c(246.9383256170808,177.59474888175154,28.39003962214503,0.0),c(568.8400715107754,612.7018836420428,564.0154538266146,270.64176251684285))
targetgene="ACIN1"
collabel=c("HL60.initial","KBM7.initial","HL60.final","KBM7.final")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(2484.0819660289676,2349.578527705573,2172.7843657481662,910.9126552363929),c(992.1629154257711,1005.1862786707138,743.8190381001997,200.26346063614164),c(1267.0287897733551,1156.1418152202027,251.09412821363824,42.34141739164138),c(1500.738276518092,1315.977089213779,800.5991173444897,1476.2277955464156),c(1925.5309914189038,2054.7712445618654,194.94493873872918,235.16652091844063),c(351.29916561001374,781.4168950797068,227.75120674654121,624.2498158686586),c(1719.74905340467,1006.9622261595313,356.45271970026533,222.0063506480656),c(903.9706562768137,1445.6212558974576,1482.5909580453515,1055.1023468944147),c(651.152846716469,1081.552020689867,576.0023594448536,1072.2677863775127),c(285.1549712482957,408.46792242802854,99.0496937928171,44.630142656054424))
targetgene="ACTR8"
collabel=c("HL60.initial","KBM7.initial","HL60.final","KBM7.final")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(301.3235520922712,657.1005708624807,228.38209651592223,137.32351586478285),c(1142.0897559789987,1099.311495578042,112.92926871919911,100.70391163417409),c(789.3207193831689,671.3081507730209,723.6305654800077,588.7745742702564),c(392.45555321286054,412.0198174056636,334.37157777193033,213.99581222261992),c(2009.3136376104133,2235.917888421252,2437.1271791188055,1937.9781176417478),c(1071.5359486598327,406.69197493921104,645.4002340767636,349.602784139093),c(61.7345814042702,218.44154112455442,614.4866353770946,452.5954210376801),c(651.152846716469,879.0940069646701,237.21455328725622,18.88198343140764),c(1625.6773103124485,1410.1023061211074,2146.286995434164,1986.613529510525),c(1053.8974968300413,882.6459019423052,106.6203710253891,105.85354347910344))
targetgene="AHCY"
collabel=c("HL60.initial","KBM7.initial","HL60.final","KBM7.final")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(1268.498660759171,1411.8782536099247,1136.2324746551822,603.6512884889412),c(327.78122983695846,454.642557137284,51.73296108924205,24.031615276336996),c(132.28838872343613,241.5288584791821,123.02350502929512,65.80085135187511),c(495.34652221997754,586.0626713097802,279.4841678357833,243.74924065998954),c(1009.8013672555626,1102.8633905556771,1237.174837756142,1004.7503910773278),c(877.5129785321263,715.7068379934587,538.1489732819936,594.496387431289),c(1594.8100196103135,1108.1912330221296,605.6541786057605,127.59643349102738),c(314.5523909646148,252.1845434120872,88.95545748272109,359.9020478289517),c(512.984974049769,269.94401830026237,205.67006481820619,126.45207085882086),c(761.3931706526657,475.9539270030942,559.5992254409475,596.7851126957021))
targetgene="ACLY"
collabel=c("HL60.initial","KBM7.initial","HL60.final","KBM7.final")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(659.9720726313648,809.832054900787,880.7221180558769,802.1982051767731),c(724.6463960072668,1086.8798631563195,695.2405258578626,307.26136674745163),c(836.3565909292796,1289.3378768815162,468.75109865008346,177.94838930811443),c(367.46774645398926,571.85509139924,300.30353022535627,116.72498848506541),c(518.8644579930328,632.2373060190355,627.7353205340956,308.9779106957614),c(405.68439208520414,259.28833336735727,324.27734146183434,166.5047629860492),c(2096.0360257735547,1960.6460276545372,1573.4390848362154,629.9716290296913),c(277.8056163192159,435.1071347602913,182.32714335110919,0.0),c(995.1026573974029,477.7298744919117,728.0467938656747,275.21921304566894),c(2185.6981559083283,1482.9161531626255,1741.8866532609427,1862.4501839161173))
targetgene="AATF"
collabel=c("HL60.initial","KBM7.initial","HL60.final","KBM7.final")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(640.8637498157573,602.0461987091378,307.2433176885473,192.82510352679924),c(354.23890758164566,280.5997032331675,204.4082852794442,275.79139436177223),c(779.0316224824572,932.3724316291956,778.5179754161547,905.1908420753603),c(624.6951689717818,554.0956165110648,370.96318439602834,558.4489645167836),c(1133.270530064103,1394.1187787217498,639.0913363829536,1131.2024619361487),c(423.32284391499564,412.0198174056636,224.59675789963623,426.84726181303336),c(296.91393913482335,829.3674772777797,489.5704610396565,1233.0507362025292),c(684.959879390236,546.9918265557948,394.30610586312537,566.4595029422292),c(440.96129574478715,630.461358530218,434.6830511035094,457.1728715665062),c(1108.2827233052317,1969.5257650986248,1066.2037102538911,1333.7546478367033))
targetgene="AGBL5"
collabel=c("HL60.initial","KBM7.initial","HL60.final","KBM7.final")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(196.96271209933826,301.9110730989776,423.9579250240324,34.33087896619571),c(1106.8128523194157,1056.6887558464218,1743.1484327997048,807.3478370217025),c(748.1643317803222,488.3855594248168,239.73811236478022,477.77139894622366),c(1095.053884432888,882.6459019423052,837.8216137379688,365.05167967388104),c(677.6105244611563,316.11865300951774,613.8557456077136,819.3636446598709),c(1078.8853035889126,1609.008424868669,348.88204246769334,193.96946615900578),c(1437.533824128006,1095.759600600407,320.4920028455483,161.35513114111984),c(845.1758168441753,660.6524658401157,541.3034221288985,640.8430740356532),c(551.2016196809839,740.570102836904,1103.42620664737,622.5332719203489),c(601.1772331987264,900.4053768304803,735.6174710982467,754.1349746240991))
targetgene="AHCTF1"
collabel=c("HL60.initial","KBM7.initial","HL60.final","KBM7.final")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(487.9971672908978,367.6211301852257,312.2904358435953,441.15179471561487),c(358.6485205390935,394.2603425174884,593.0363832181406,268.35303725242983),c(1743.266989177725,1980.1814500315297,837.1907239685878,281.5132075228048),c(1597.7497615819454,1465.1566782744503,1065.57282048451,992.7345834391593),c(119.05954985109253,378.2768151181308,185.48159219801417,128.7407961232339),c(986.2834314825072,745.8979453033566,328.0626800781203,302.11173490252224),c(523.2740709504807,694.3954681276485,336.89513684945433,597.9294753279087),c(1562.4728579223624,763.6574201915316,422.0652557158894,220.28980669975581),c(30.8672907021351,179.37069637056908,238.47633282601822,184.81456510135357),c(339.5401977234861,447.5387671820139,310.3977665354523,205.98527379717427))
targetgene="ABT1"
collabel=c("HL60.initial","KBM7.initial","HL60.final","KBM7.final")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(492.4067802483456,221.99343610218943,309.7668767660713,102.99263689858714),c(243.9985836454489,239.7529109903646,130.59418226186713,174.51530141149487),c(734.9354929079785,673.0840982618383,620.7955330709046,470.9052231529845),c(1074.4756906314647,950.1319065173708,1100.902647569846,743.8357109342404),c(702.5983312200275,1010.5141211371663,1291.4313579229083,1017.3383800315995),c(1647.7253750996879,760.1055252138966,685.7771793171477,608.2287390177673),c(951.0065278229242,864.8864270541301,606.9159581445226,769.0116888427839),c(435.0818118015233,435.1071347602913,275.69882921949727,339.8757017653375),c(89.66213013477338,209.56180368046682,208.8245136651112,304.4004601669353),c(1328.7633711776252,1571.7135276035012,1122.983789498181,1356.6419004808338))
targetgene="ADIRF"
collabel=c("HL60.initial","KBM7.initial","HL60.final","KBM7.final")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(216.0710349149457,289.479440677255,192.42137966120518,498.36992632594104),c(1127.391046120839,1198.764554951823,371.5940741654094,370.2013115188104),c(1111.2224652768637,1038.9292809582466,948.227323379644,922.3562815584581),c(1164.137820766238,1204.0923974182756,1686.9992433247955,2089.033985093009),c(48.505742531926586,248.63264843445216,665.5887066969557,248.8988725049189),c(501.2260061632414,387.1565525622184,436.5757204116524,314.69972385679404),c(1975.5066049366465,1797.2588586833258,1628.3264947723626,1289.6966864967521),c(213.13129294331378,376.5008676293133,404.4003421732214,482.921030791153),c(2012.2533795820452,1989.0611874756173,1064.3110409457481,431.9968936579627),c(264.57677744687226,353.4135502746856,442.25372833608145,191.6807408945927))
targetgene="ABCF1"
collabel=c("HL60.initial","KBM7.initial","HL60.final","KBM7.final")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}



# 
#
# parameters
# Do not modify the variables beginning with "__"

# gstablename='__GENE_SUMMARY_FILE__'
startindex=9
# outputfile='__OUTPUT_FILE__'
targetgenelist=c("ACRC","AGAP3","ADCK4","AHRR","ADRBK1","ADK","ADCK1","ADARB2","ACSS2","ADNP")
# samplelabel=sub('.\\w+.\\w+$','',colnames(gstable)[startindex]);
samplelabel='HL60.final,KBM7.final_vs_HL60.initial,KBM7.initial pos.'


# You need to write some codes in front of this code:
# gstable=read.table(gstablename,header=T)
# pdf(file=outputfile,width=6,height=6)


# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")

######
# function definition

plotrankedvalues<-function(val, tglist, ...){
  
  plot(val,log='y',ylim=c(max(val),min(val)),type='l',lwd=2, ...)
  if(length(tglist)>0){
    for(i in 1:length(tglist)){
      targetgene=tglist[i];
      tx=which(names(val)==targetgene);ty=val[targetgene];
      points(tx,ty,col=colors[(i %% length(colors)) ],cex=2,pch=20)
      # text(tx+50,ty,targetgene,col=colors[i])
    }
    legend('topright',tglist,pch=20,pt.cex = 2,cex=1,col=colors)
  }
}



plotrandvalues<-function(val,targetgenelist, ...){
  # choose the one with the best distance distribution
  
  mindiffvalue=0;
  randval=val;
  for(i in 1:20){
    randval0=sample(val)
    vindex=sort(which(names(randval0) %in% targetgenelist))
    if(max(vindex)>0.9*length(val)){
      # print('pass...')
      next;
    }
    mindiffind=min(diff(vindex));
    if (mindiffind > mindiffvalue){
      mindiffvalue=mindiffind;
      randval=randval0;
      # print(paste('Diff: ',mindiffvalue))
    }
  }
  plot(randval,log='y',ylim=c(max(randval),min(randval)),pch=20,col='grey', ...)
  
  if(length(targetgenelist)>0){
    for(i in 1:length(targetgenelist)){
      targetgene=targetgenelist[i];
      tx=which(names(randval)==targetgene);ty=randval[targetgene];
      points(tx,ty,col=colors[(i %% length(colors)) ],cex=2,pch=20)
      text(tx+50,ty,targetgene,col=colors[i])
    }
  }
  
}




# set.seed(1235)



pvec=gstable[,startindex]
names(pvec)=gstable[,'id']
pvec=sort(pvec);

plotrankedvalues(pvec,targetgenelist,xlab='Genes',ylab='RRA score',main=paste('Distribution of RRA scores in \\n',samplelabel))

# plotrandvalues(pvec,targetgenelist,xlab='Genes',ylab='RRA score',main=paste('Distribution of RRA scores in \\n',samplelabel))


pvec=gstable[,startindex+1]
names(pvec)=gstable[,'id']
pvec=sort(pvec);

plotrankedvalues(pvec,targetgenelist,xlab='Genes',ylab='p value',main=paste('Distribution of p values in \\n',samplelabel))

# plotrandvalues(pvec,targetgenelist,xlab='Genes',ylab='p value',main=paste('Distribution of p values in \\n',samplelabel))



# you need to write after this code:
# dev.off()






# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(461.5394895462105,502.5931393353569,445.40817718298644,889.1697652244688),c(76.43329126242978,90.5733219296933,447.30084649112945,357.0411412484354),c(258.6972935036084,685.515730683561,533.7327448963265,560.7376897811967),c(232.23961575892122,681.9638357059259,275.69882921949727,467.47213525636494),c(1393.4376945535273,1472.2604682297203,1039.706339939889,532.7008052921368),c(2395.88970688001,2441.927797124084,2462.9936596634266,2461.5240218762324),c(495.34652221997754,605.5980936867728,1159.575396122279,1617.5565806239213),c(682.0201374186041,822.2636873225097,1572.1773052974536,1333.7546478367033),c(961.2956247236359,1097.5355480892247,959.5833392285019,905.1908420753603),c(1940.2297012770634,1289.3378768815162,942.5493154552149,1103.737758763192))
targetgene="ACRC"
collabel=c("HL60.initial","KBM7.initial","HL60.final","KBM7.final")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(1387.5582106102636,1120.6228654438523,1214.4628060584262,1111.1761158725344),c(388.0459402554127,509.69692929062694,933.0859689144999,750.1297054113762),c(326.3113588511425,635.7892009966705,960.8451187672639,615.6670961271097),c(1328.7633711776252,1038.9292809582466,1346.3187678590552,1596.3858719281006),c(352.7690365958297,234.42506852391205,310.3977665354523,429.1359870774464),c(693.7791053051318,678.4119407282909,784.1959833405838,895.4637597016048),c(837.8264619150956,719.2587329710938,374.74852301231437,993.8789460713658),c(365.99787546817333,369.3970776740432,333.74068800254935,746.6966175147567),c(707.0079441774753,635.7892009966705,837.1907239685878,1465.3563505404536),c(486.5272963050818,673.0840982618383,784.8268731099647,734.6808098765882))
targetgene="AGAP3"
collabel=c("HL60.initial","KBM7.initial","HL60.final","KBM7.final")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(830.4771069860158,864.8864270541301,1349.4732167059603,740.974804353724),c(1481.6299537024847,1994.38902994207,2044.082852794442,1810.9538654668238),c(1234.6916280854039,1299.9935618144214,1357.6747837079133,2232.6514954349277),c(224.89026082984142,188.25043381465665,700.2876440129107,81.24974688666317),c(812.8386551562243,845.3510046771374,946.334654071501,999.6007592323984),c(1978.4463469082782,1751.0842239740703,2659.2003779409174,2851.1794981425537),c(565.9003295391435,776.0890526132542,878.1985589783528,445.72924524444096),c(680.5502664327881,534.5601941340722,550.7667686696135,1025.9210997731484),c(161.68580843975528,333.87812789769293,275.0679394501163,465.18340999195186),c(2523.768482645998,2445.4796921017187,2153.226782897355,1516.8526689897471))
targetgene="ADCK4"
collabel=c("HL60.initial","KBM7.initial","HL60.final","KBM7.final")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(345.4196816667499,163.38716897121142,474.42910657451245,481.2044868428432),c(415.9734889859159,372.9489726516783,212.6098522813972,349.03060282298975),c(1.469870985815957,83.46953197442323,0.0,62.9399447713588),c(351.29916561001374,150.9555365494888,288.9475143764983,416.54799812317464),c(561.4907165816957,170.49095892648148,199.3611671243962,411.97054759434855),c(251.34793857452865,221.99343610218943,1564.6066280648815,1502.5481360871656),c(736.4053638937945,893.3015868752103,1114.782222496228,459.46159683091923),c(338.07032673767014,607.3740411755903,378.5338616286004,65.22867003577186),c(1230.2820151279561,525.6804566899846,837.1907239685878,945.2435342025885))
targetgene="AHRR"
collabel=c("HL60.initial","KBM7.initial","HL60.final","KBM7.final")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(371.87735941143717,877.3180594758527,2395.4884543396593,1564.9158995424211),c(1109.7525942910477,1138.3823403320275,970.308465307979,999.0285779162951),c(1462.5216308868773,1209.420239884728,1537.4783679814984,1519.14139425416),c(586.4785233405669,987.4268037825386,743.8190381001997,1312.0117578247794),c(1018.6205931704583,717.4827854822763,1070.619938639558,1144.3626322065236),c(1269.9685317449869,1212.9721348623632,1591.1039983788835,1624.9949377332637),c(1321.4140162485455,1795.4829111945082,1478.8056194290655,1237.056005415252),c(908.3802692342615,832.9193722554148,1639.6825106212207,1268.5259778009315),c(923.078979092421,758.3295777250792,1479.4365091984464,1275.964334910274),c(680.5502664327881,634.013253507853,318.5993335374053,631.1159916618979))
targetgene="ADRBK1"
collabel=c("HL60.initial","KBM7.initial","HL60.final","KBM7.final")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(1472.810727787589,1829.225913482041,1263.0413183007631,1315.444845721399),c(208.7216799858659,65.71005708624807,292.1019632234033,350.17496545519623),c(1011.2712382413785,1166.7975001531076,652.9709113093356,860.5606994193058),c(557.0811036242477,685.515730683561,875.0441101314478,1019.6271052960126),c(363.0581334965414,825.8155823001447,736.8792506370087,349.602784139093),c(1505.14788947554,451.09066215964896,653.6018010787167,991.0180394908496),c(198.43258308515422,28.41515982108025,249.83234867487624,114.43626322065236),c(438.02155377315523,74.58979453033565,254.87946682992424,231.16125170571777),c(804.0194292413286,472.4020320254591,1336.2245315489592,1203.2973077651598),c(454.19013461713075,490.1615069136343,896.4943622904019,685.4732166917076))
targetgene="ADK"
collabel=c("HL60.initial","KBM7.initial","HL60.final","KBM7.final")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(662.9118146029966,1008.7381736483488,1101.533537339227,1694.8010582978616),c(1547.7741480642028,1965.9738701209897,1869.9572764452857,2353.9539344488194),c(1459.5818889152454,1179.2291325748304,1296.4784760779562,1222.1792911965672),c(1193.5352404825571,1355.0479339677643,1622.0175970785526,1905.9359639399652),c(868.6937526172306,701.4992580829187,720.4761166331027,603.6512884889412),c(798.1399452980647,768.9852626579842,1478.8056194290655,1756.0244591209105),c(1168.5474337236858,907.5091667857504,879.4603385171149,977.8578692204745),c(809.8989131845924,687.2916781723785,678.8373918539567,865.7103312642352),c(1246.4505959719315,753.0017352586266,1301.5255942330043,1264.5207085882087),c(826.0674940285679,797.4004224790644,977.8791425405509,2066.7189137649816))
targetgene="ADCK1"
collabel=c("HL60.initial","KBM7.initial","HL60.final","KBM7.final")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(1863.7964100146337,1585.9211075140413,1761.4442361117538,1464.211987908247),c(742.2848478370584,598.4943037315028,943.8110949939769,820.5080072920774),c(1568.3523418656262,2083.1864043829455,1810.6536381234716,1887.6261618246608),c(1018.6205931704583,513.248824268262,679.4682816233377,824.5132765048003),c(1140.6198849931827,1191.6607649965529,880.0912282864958,977.8578692204745),c(135.22813069506805,118.98848175077354,351.40560154521734,399.95473995618005),c(665.8515565746286,701.4992580829187,986.7115993118849,746.6966175147567),c(418.9132309575478,300.1351256101601,376.6411923204574,645.4205245644794),c(561.4907165816957,543.4399315781598,881.9838975946388,580.7640358448108),c(442.4311667306031,229.0972260574595,395.5678854018874,651.142337725512))
targetgene="ADARB2"
collabel=c("HL60.initial","KBM7.initial","HL60.final","KBM7.final")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(734.9354929079785,358.74139274113816,541.9343118982795,378.7840312603593),c(595.2977492554626,591.3905137762326,1061.787481868224,887.4532212761591),c(1655.0747300287676,943.0281165621008,1069.358159100796,2038.1098479598186),c(626.1650399575977,884.4218494311227,517.3296108924205,858.2719741548927),c(680.5502664327881,747.673892792174,533.1018551269456,1016.194017399393),c(662.9118146029966,777.8650001020718,864.9498738213518,787.3214909580882),c(880.4527205037583,621.5816210861304,671.8976043907657,1040.7978139918332),c(94.07174309222125,447.5387671820139,711.6436598617687,927.5059134033875),c(399.80490814194036,806.280159923152,1147.58849050404,1059.1076161071376),c(698.1887182625796,531.0082991564371,504.0809257354195,347.8862401907832))
targetgene="ACSS2"
collabel=c("HL60.initial","KBM7.initial","HL60.final","KBM7.final")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(408.62413405683606,523.9045092011671,483.89245311522745,701.494293542599),c(1805.0015705819953,1434.9655709645526,1712.2348341000356,2152.546111180471),c(3017.64513388016,2642.609863360463,1834.6274493599499,3573.2723190648703),c(1649.1952460855039,783.1928425685244,773.4708572611067,1332.0381038883936),c(959.82575373782,1397.6706736993847,1429.5962174173474,2811.126806015325),c(495.34652221997754,301.9110730989776,336.89513684945433,555.015876620164),c(1491.9190506031964,1331.9606166131366,2087.614246881731,1983.1804416139055),c(429.2023278582595,889.7496918975753,567.8007924429005,1132.9190058844583),c(427.7324568724435,573.6310388880576,655.4944703868597,899.4690289143276),c(690.8393633334998,767.2093151691668,1040.33722970927,993.3067647552625))
targetgene="ADNP"
collabel=c("HL60.initial","KBM7.initial","HL60.final","KBM7.final")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}



dev.off()
Sweave("output_summary.Rnw");
library(tools);

texi2dvi("output_summary.tex",pdf=TRUE);

