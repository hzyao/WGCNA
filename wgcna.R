

######################## WGCNA - 加权基因共表达网络分析 ########################


# 安装并加载包
library(WGCNA)

# WGCNA依赖的包比较多，bioconductor上的包需要自己安装，cran上依赖的包可以自动安装。
# 通常在R中运行下面4条语句应该就可以完成WGCNA的安装。
# source("https://bioconductor.org/biocLite.R")
# biocLite(c("AnnotationDbi", "impute","GO.db", "preprocessCore"))
# site = "https://mirrors.tuna.tsinghua.edu.cn/CRAN"
# install.packages(c("WGCNA", "stringr", "reshape2"), repos = site)



# 加载数据，这里包括需要用到的表达矩阵和临床信息

# 表达矩阵：常规表达矩阵即可，我使用的数据已经经过了log2(TPM+1)处理
exp_tpm <- readRDS("./data/exp_tpm.rds")
head(exp_tpm)[1:5, 1:5]
#              TSPAN6 TNMD     DPM1    SCYL3 C1orf112
# ACH-001113 4.331992    0 7.364660 2.792855 4.471187
# ACH-000242 6.729417    0 6.537917 2.456806 3.867896
# ACH-000461 4.017031    0 6.534497 2.226509 3.021480
# ACH-000528 4.512227    0 7.099926 2.843984 4.672991
# ACH-000792 3.280956    0 6.391115 1.752749 3.436961

# 如果有批次效应，需要先进行去除；
# 如果数据存在系统偏移，需要进行quantile normalization；
# 标准化推荐使用DESeq2中的varianceStabilizingTransformation方法，
# 或将基因标准化后的数据（如FPKM、CPM等）进行log2(x+1)转化。
# （大家可以依据自己数据的实际情况进行这些步骤）

# 顺便打个广告再！
# 想详细了解FPKM、TPM、CPM等数据类型的小伙伴们，
# 可以查看：https://mp.weixin.qq.com/s/j9woMCGVoaw52x-ppyjFMw



# 判断一下数据质量，其实就是检测缺失值
# 通过goodSamplesGenes函数生成一个质量检查结果对象gsg
gsg <- goodSamplesGenes(exp_tpm, verbose = 3)

# `verbose` 参数被设置为 `3`，用于控制函数 `goodSamplesGenes` 的详细输出程度。
# `verbose = 0`：不产生任何输出，只返回结果，通常用于静默模式。
# `verbose = 1`：产生基本的信息输出，以提供一些关于函数执行进度的信息。
# `verbose = 2`：产生更详细的输出，可能包括一些中间步骤的信息。
# `verbose = 3`：产生最详细的输出，通常包括每个步骤的详细信息，用于调试和详细分析。


# 检查质量检查结果对象中的allOK属性，如果为FALSE，表示质量不好
gsg$allOK
# [1] FALSE

# 哎呀！质量不好，咱们执行以下代码块
if (!gsg$allOK) 
{ 
  # 如果goodGenes属性中有不好的基因，打印并移除它们
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", paste(names(exp_tpm)[!gsg$goodGenes], collapse = ", "))) 
  
  # 如果goodSamples属性中有不好的样本，打印并移除它们
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", paste(rownames(exp_tpm)[!gsg$goodSamples], collapse = ", "))) 
  
  # 根据质量检查结果对象中的goodSamples和goodGenes属性，过滤数据集exp_tpm
  exp_tpm <- exp_tpm[gsg$goodSamples, gsg$goodGenes] 
}
# Removing genes: DEFB131A, DEFB110, FAM236D, USP17L5, USP17L10, FAM236C

# 我们删除了以上这几个基因，其实一般大家都不会需要这个步骤，因为咱自己的数据，咱应该心里有数叭哈哈哈哈哈哈哈

nGenes = ncol(exp_tpm)
nGenes
# [1] 19187

nSamples = nrow(exp_tpm)
nSamples
# [1] 476

dim(exp_tpm)
# [1]   476 19187




# 临床信息：这个一般建议大家自己整理，主要包括你想研究的性状或表型，其中每一列代表一个性状，而每一行代表一个样本
drug_screen_parp <- readRDS("./data/drug_screen_parp.rds")
head(drug_screen_parp)
#                  auc     ic50      ec50
# ACH-000007 0.9910031 1.419984 0.7635118
# ACH-000008 0.9635668 4.008706 2.6297078
# ACH-000011 0.9247148 2.114563 1.2638989
# ACH-000012 0.8990521 2.407329 1.5441064
# ACH-000013 0.7971138 3.468354 1.8830208
# ACH-000014 0.8965265 1.909032 1.1196910

# 因为我们今天关注的是样本的PARP抑制剂敏感性，大家可以看到，我准备的临床信息中，
# 有 auc、ec50、ic50等等，这些都是用于判断样本是否对某药物敏感的指标。
# 有想详细了解的小伙伴可以告诉我，后期我可以专门出一期！

# 这里要注意的是，用于关联分析的性状必须是数值型特征（比如年龄、体重、基因表达水平等）。
# 如果是分类变量（比如性别、地区等），则需要转换为0-1矩阵的形式
# （1表示属于此组或有此属性，0表示不属于此组或无此属性）。




# 第一步，咱们就是构建样本的系统聚类树
# 主要为了查看是否有离群样本
sampleTree <- hclust(dist(exp_tpm), method = "average")

# sampleTree是层次聚类的结果，通常是一个树状结构，其中包含了样本的聚类信息。
# 这个树可以是通过层次聚类算法生成的，用于将样本分组成簇。

# 可视化
plot(sampleTree, main = "Sample clustering to detect outliers", sub = "", xlab = "")

# 如果数据中有个别数据与其他数据有较大的差异性，为了让后续分析更精准，可以将其进行适当的删除即可。

# 可以在图中画线，剔除离群样本，注意：不需要剔除样本画线就别运行这行！
abline(h = 200, col = "red")

# 将层次聚类树中的样本分成不同的簇
clust <- cutreeStatic(sampleTree, cutHeight = 200, minSize = 10)

# cutHeight = 200：用于指定在层次聚类树中切割的高度。在树状结构中，高度表示样本之间的相似性或距离。
# 通过指定 cutHeight，你可以控制在哪个高度水平切割树，从而确定最终的簇数。

# minSize = 10：用于指定最小簇的大小。在进行切割时，如果某个簇的大小小于 minSize，
# 则可能会合并到其他簇中，以确保生成的簇都具有足够的样本数。

# 查看详情
table(clust)
# clust
#   0   1 
#   9 467 


# 剔除离群样本，保留我们想要留下的样本们！
keepSamples <- clust == 1
datExpr <- exp_tpm[keepSamples, ]

nGenes = ncol(datExpr)
nGenes
# [1] 19187

nSamples = nrow(datExpr)
nSamples
# [1] 467

dim(datExpr)
# [1]   467 19187

# 可以看到，我们已经剔除了9个离群样本。

# 啊，和临床信息再对应一下
datExpr <- datExpr %>% filter(rownames(datExpr) %in% rownames(drug_screen_parp))

nSamples = nrow(datExpr)
nSamples
# [1] 464

dim(datExpr)
# [1]   464 19187


# 剔除离群样本后，咱们可以再重新聚类一下
sampleTree_final <- hclust(dist(datExpr), method = "average")

# 我们前面已经整理好了临床信息，咱们可以再回忆一下
head(drug_screen_parp)
#                  auc     ic50      ec50
# ACH-000007 0.9910031 1.419984 0.7635118
# ACH-000008 0.9635668 4.008706 2.6297078
# ACH-000011 0.9247148 2.114563 1.2638989
# ACH-000012 0.8990521 2.407329 1.5441064
# ACH-000013 0.7971138 3.468354 1.8830208
# ACH-000014 0.8965265 1.909032 1.1196910

# 为了演示过程更便于理解，咱们就以auc作为关注的性状


# 绘制样本聚类图（上）与样本性状热图（下）

# 这一步是为了将数值转换为颜色编码，用于可视化，将数值映射到颜色，我们可以更直观地展示数据的特征或模式
traitColors <- numbers2colors(drug_screen_parp$auc, signed = FALSE)
head(traitColors)
#      [,1]     
# [1,] "#FFB4A1"
# [2,] "#FFBCAC"
# [3,] "#FFC9BB"
# [4,] "#FFD1C6"
# [5,] "#FFF3EF"
# [6,] "#FFD1C6"

# 可以看到所有的数值都转换为了颜色编码

# 可视化
plotDendroAndColors(sampleTree_final, traitColors,
                    groupLabels = names(drug_screen_parp),
                    main = "Sample dendrogram and trait heatmap")

# 我们可以保存一下处理好的表达矩阵和临床信息（也就是所谓的性状矩阵）。
saveRDS(datExpr, file = "./data/datExpr.rds")


# 挑选最佳软阈值

# 设置 power 参数选择范围，可以自行修改设定
powers <- c(c(1:10), seq(from = 12, to = 20, by = 2))

# 选择最佳软阈值，获取各个阈值下的 R^2 和平均连接度
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

# 我们看一下软阈值选择结果
sft
# $powerEstimate
# [1] 5
# 
# $fitIndices
#    Power  SFT.R.sq     slope truncated.R.sq      mean.k.    median.k.     max.k.
# 1      1 0.3875076  1.704334      0.9129145 1778.5099895 1.786517e+03 3556.81883
# 2      2 0.4961936 -1.538776      0.9538743  302.2895064 2.759378e+02 1021.10668
# 3      3 0.7696795 -1.984592      0.9984766   72.2718042 5.734930e+01  379.97377
# 4      4 0.8484691 -2.341492      0.9914478   21.9532922 1.454930e+01  186.06380
# 5      5 0.9581711 -2.036913      0.9590233    8.1076302 4.377465e+00  109.18653
# 6      6 0.9803518 -1.900979      0.9786923    3.5676380 1.521133e+00   95.60116
# 7      7 0.9768706 -1.684097      0.9884073    1.8428375 5.933449e-01   85.28786
# 8      8 0.9664818 -1.517685      0.9842268    1.0954397 2.505880e-01   77.09653
# 9      9 0.9644630 -1.410692      0.9794021    0.7300586 1.132233e-01   70.37984
# 10    10 0.9645054 -1.323333      0.9790107    0.5307152 5.341439e-02   64.74042
# 11    12 0.9640755 -1.236659      0.9738519    0.3327212 1.298621e-02   55.73650
# 12    14 0.9713689 -1.176731      0.9800533    0.2380254 3.576691e-03   48.81645
# 13    16 0.9728194 -1.151938      0.9806096    0.1827955 1.028644e-03   43.30473
# 14    18 0.9807890 -1.131848      0.9851551    0.1464499 3.097658e-04   38.80021
# 15    20 0.9828774 -1.123491      0.9829834    0.1206651 9.463552e-05   35.04578


# 可视化
sizeGrWindow(9, 5)
par(mfrow = c(1,2))
cex1 = 0.9

# 绘制软阈值和拟合指标的关系图
plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit, signed R^2", type = "n",
     main = paste("Scale independence"))
text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, cex = cex1, col = "red")

# 添加 R^2 水平线，使用 R^2 阈值为 0.90，官网建议最好是0.85或以上
abline(h = 0.90, col = "red")

# 绘制软阈值对平均连接度的影响图
plot(sft$fitIndices[, 1], sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)", ylab = "Mean Connectivity", type = "n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = powers, cex = cex1, col = "red")




# 正式构建加权共表达网络

# 一步法，这步时间可能有点久
net <- blockwiseModules(datExpr, power = sft$powerEstimate, 
                       maxBlockSize = nGenes, TOMType = "unsigned", 
                       minModuleSize = 30, reassignThreshold = 0, mergeCutHeight = 0.25, 
                       numericLabels = TRUE, pamRespectsDendro = FALSE, 
                       saveTOMs = F, verbose = 3)
# Calculating module eigengenes block-wise from all genes
# Flagging genes and samples with too many missing values...
# ..step 1
# ..Excluding 2 genes from the calculation due to too many missing samples or zero variance.
# ..step 2
# ..Working on block 1 .
# TOM calculation: adjacency..
# ..will not use multithreading.
# Fraction of slow calculations: 0.000000
# ..connectivity..
# ..matrix multiplication (system BLAS)..
# ..normalization..
# ..done.
# ....clustering..
# ....detecting modules..
# ....calculating module eigengenes..
# ....checking kME in modules..
# ..removing 1148 genes from module 1 because their KME is too low.
# ..removing 369 genes from module 2 because their KME is too low.
# ..removing 22 genes from module 3 because their KME is too low.
# ..removing 217 genes from module 4 because their KME is too low.
# ..removing 434 genes from module 5 because their KME is too low.
# ..removing 36 genes from module 6 because their KME is too low.
# ..removing 60 genes from module 7 because their KME is too low.
# ..removing 14 genes from module 8 because their KME is too low.
# ..removing 12 genes from module 9 because their KME is too low.
# ..removing 55 genes from module 10 because their KME is too low.
# ..removing 70 genes from module 11 because their KME is too low.
# ..removing 1 genes from module 12 because their KME is too low.
# ..removing 2 genes from module 13 because their KME is too low.
# ..removing 3 genes from module 14 because their KME is too low.
# ..removing 24 genes from module 16 because their KME is too low.
# ..removing 19 genes from module 17 because their KME is too low.
# ..removing 17 genes from module 18 because their KME is too low.
# ..removing 41 genes from module 19 because their KME is too low.
# ..removing 2 genes from module 23 because their KME is too low.
# ..removing 7 genes from module 24 because their KME is too low.
# ..removing 10 genes from module 25 because their KME is too low.
# ..removing 1 genes from module 28 because their KME is too low.
# ..removing 2 genes from module 32 because their KME is too low.
# ..removing 5 genes from module 34 because their KME is too low.
# ..merging modules that are too close..
# mergeCloseModules: Merging modules whose distance is less than 0.25
# Calculating new MEs...

# 查看模块情况
table(net$colors) 
#    0    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16   17   18   19   20   21   22 
# 7222 2758 2127 1889  860  536  451  429  378  312  270  244  203  170  143  105  101   84   81   70   67   61   60 
#   23   24   25   26   27   28   29   30   31   32   33 
#   60   58   58   55   52   52   50   50   47   47   37

# 根据模块中基因数目的多少，降序排列，依次编号，比如 1 为最大模块，模块中基因最多。
# 0 表示没有分入任何模块的基因。 


# 使用层次聚类树展示各个模块，模块可视化

# 将基因模块的标签转换为对应的颜色编码
moduleColors <- labels2colors(net$colors)
table(moduleColors)
# moduleColors
#      black           blue          brown           cyan      darkgreen       darkgrey darkolivegreen 
#        429           2127           1889            143             60             58             37 
# darkorange        darkred  darkturquoise          green    greenyellow           grey         grey60 
#         55             61             60            536            244           7222             84 
#  lightcyan     lightgreen    lightyellow        magenta   midnightblue         orange  paleturquoise 
#        101             81             70            312            105             58             47 
#       pink         purple            red      royalblue    saddlebrown         salmon        skyblue 
#        378            270            451             67             50            170             52 
#  steelblue            tan      turquoise         violet          white         yellow 
#         50            203           2758             47             52            860 


# 绘制层次聚类树
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

# 这里用不同的颜色来代表那些所有的模块，其中灰色默认是无法归类于任何模块的那些基因，
# 如果灰色模块里面的基因太多，那么前期对表达矩阵挑选基因的步骤可能就不太合适。


# 分步法
# 暂略


# 我们把刚刚得到的一些数据保存一下！防止服务器炸了咱们哭唧唧！

# 颜色标签
moduleLables <- net$colors
moduleColors <- labels2colors(net$colors)


# ME值，也就是获取eigengenes，每个ME代表一个模块
MEs <- net$MEs
head(MEs)[1:5, 1:5]
#                      ME7         ME5         ME12         ME17        ME13
# ACH-001113 -0.0007874971 -0.01357137  0.087065262  0.064226471 -0.02089271
# ACH-000242 -0.0012565168 -0.02072050  0.054847155 -0.005907044 -0.01726976
# ACH-000461 -0.0050758474 -0.03310460 -0.010548236 -0.003185603 -0.01910437
# ACH-000528 -0.0056703988 -0.04616061 -0.008741992 -0.016184701 -0.02564138
# ACH-000792 -0.0047607635  0.03734415 -0.007612451 -0.008838765  0.01790916


geneTree <- net$dendrograms[[1]]

save(modduleLables, moduleColors, MEs, geneTree, file = "./data/networkConstruction.RData")



# 将基因模块与性状进行关联

# 获取eigengenes，用颜色标签计算ME值
MEList <-  moduleEigengenes(datExpr, colors = moduleColors)
MEs0 <- MEList$eigengenes

# 查看用颜色标签计算的ME值
head(MEs0)[1:5, 1:5]
#                  MEblack       MEblue      MEbrown       MEcyan  MEdarkgreen
# ACH-001113 -0.0007874971 -0.007654900  0.036397596  0.003731432 -0.013282549
# ACH-000242 -0.0012565168  0.005612039  0.005233571 -0.023237757 -0.011480834
# ACH-000461 -0.0050758474 -0.016010056 -0.005049884 -0.016021361 -0.010828976
# ACH-000528 -0.0056703988  0.002251787  0.027288677  0.024569227 -0.001865334
# ACH-000792 -0.0047607635  0.003413084 -0.038058864  0.012558507 -0.006653088

# 可以看到我们的列名已经变成了颜色，不同的颜色代表不同的模块

# 排序
MEs <- orderMEs(MEs0)
head(MEs)[1:5, 1:5]
#                  MEblack     MEgreen     MEgrey60        MEtan    MEsalmon
# ACH-001113 -0.0007874971 -0.01357137  0.064226471  0.087065262 -0.02089271
# ACH-000242 -0.0012565168 -0.02072050 -0.005907044  0.054847155 -0.01726976
# ACH-000461 -0.0050758474 -0.03310460 -0.003185603 -0.010548236 -0.01910437
# ACH-000528 -0.0056703988 -0.04616061 -0.016184701 -0.008741992 -0.02564138
# ACH-000792 -0.0047607635  0.03734415 -0.008838765 -0.007612451  0.01790916


# 计算每个模块和每个性状之间的相关性
moduleTraitCor <- cor(MEs, drug_screen_parp , use = "p");
head(moduleTraitCor)
#                  auc          ic50         ec50
# MEblack   0.04480662 -0.0007293297 -0.005612265
# MEgreen  -0.01683223  0.0308188806 -0.009436863
# MEgrey60  0.03662910 -0.0068789966 -0.012882025
# MEtan    -0.03372756 -0.0117766639 -0.025682412
# MEsalmon  0.07163196  0.0002345794  0.153468410
# MEpink    0.03327284 -0.0216032637 -0.010037838


# 计算显著性
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)
head(moduleTraitPvalue)
#                auc      ic50         ec50
# MEblack  0.3355252 0.9874994 0.9040342072
# MEgreen  0.7176324 0.5078266 0.8393441457
# MEgrey60 0.4311935 0.8825162 0.7819718265
# MEtan    0.4685990 0.8002670 0.5810751688
# MEsalmon 0.1233612 0.9959792 0.0009113246
# MEpink   0.4746191 0.6425421 0.8292653579


# 可视化相关性和P值
textMatrix <- paste(signif(moduleTraitCor, 2), "\n(", signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) <- dim(moduleTraitCor)
par(mar = c(1, 15, 3, 3))

# 绘制热图
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = colnames(drug_screen_parp),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.7,
               # 这里大家要注意，一般会设置为zlim = c(-1, 1)，怪我怪我选的这个例子不好，性状和模块之间相关性太低啦！
               # 不过如果你的相关性也很低，咱们也可以适当调整哟！
               zlim = c(-0.25, 0.25), 
               main = paste("Module-trait relationships"))
# Warning message:
#   In greenWhiteRed(50) :
#   WGCNA::greenWhiteRed: this palette is not suitable for people
# with green-red color blindness (the most common kind of color blindness).
# Consider using the function blueWhiteRed instead.

# 哈哈哈哈哈哈哈哈哈哈哈哈！这一步如果colors = greenWhiteRed(50)的话，会有这个警告！
# 好可爱！好体贴！咱们要关爱红绿色盲的人儿！
# 所以我就乖乖听话把colors = greenWhiteRed(50)改成了colors = blueWhiteRed(50)！


table(moduleColors)
moduleColors
#      black           blue          brown           cyan      darkgreen       darkgrey darkolivegreen 
#       c429           2127           1889            143             60             58             37 
# darkorange        darkred  darkturquoise          green    greenyellow           grey         grey60 
#         55             61             60            536            244           7222             84 
#  lightcyan     lightgreen    lightyellow        magenta   midnightblue         orange  paleturquoise 
#        101             81             70            312            105             58             47 
#       pink         purple            red      royalblue    saddlebrown         salmon        skyblue 
#        378            270            451             67             50            170             52 
# csteelblue            tan      turquoise         violet          white         yellow 
#         50            203           2758             47             52            860 





# 模块内基因与表型数据关联

# 我们发现与 auc 最相关的是 saddlebrown 模块
# names (colors) of the modules

# datExpr 表示每个基因在每个样本中的表达量
# MEs 表示每个模块在每个样本中的模块特征值
# moduleColors 表示每个基因所属的模块颜色

# 获取模块名称
modNames <- substring(names(MEs), 3)
modNames
# [1]  "black"          "green"          "grey60"         "tan"            "salmon"         "pink"          
# [7]  "steelblue"      "orange"         "darkgreen"      "greenyellow"    "white"          "darkred"       
# [13] "violet"         "darkturquoise"  "royalblue"      "lightgreen"     "skyblue"        "magenta"       
# [19] "cyan"           "blue"           "brown"          "turquoise"      "yellow"         "darkolivegreen"
# [25] "lightyellow"    "paleturquoise"  "saddlebrown"    "midnightblue"   "red"            "darkgrey"      
# [31] "purple"         "darkorange"     "lightcyan"      "grey" 



# 计算模块与基因的相关性矩阵
geneModuleMembership <- as.data.frame(cor(datExpr, MEs, use = "p"))
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) <- paste("MM", modNames, sep = "")
names(MMPvalue) <- paste("p.MM", modNames, sep = "");
geneModuleMembership[1:5, 1:5]
#               MMblack     MMgreen    MMgrey60        MMtan     MMsalmon
# TSPAN6    0.061400357  0.03516486 -0.03010767  0.157223072  0.057596458
# TNMD      0.008277116 -0.10551856 -0.01892595 -0.052500790 -0.017606553
# DPM1     -0.015979075  0.10719315 -0.02328479  0.044476432 -0.134685371
# SCYL3     0.075200100 -0.08775867  0.02600854 -0.005045328 -0.007943488
# C1orf112 -0.017742573 -0.07316337 -0.09035519 -0.123376142 -0.221407856

# 计算性状与基因的相关性矩阵 
# 只有连续型性状才能进行计算，如果是离散变量，在构建样本表时就转为0-1矩阵。
geneTraitSignificance <- as.data.frame(cor(datExpr, drug_screen_parp, use = "p"))
GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) <- paste("GS.", names(drug_screen_parp), sep = "")
names(GSPvalue) <- paste("p.GS.", names(drug_screen_parp), sep = "")
head(geneTraitSignificance)
#               GS.auc       GS.ic50      GS.ec50
# TSPAN6    0.07739899 -0.0104488289  0.003197307
# TNMD      0.06665378 -0.0070637254 -0.007063725
# DPM1      0.03264356 -0.0009365489 -0.063931003
# SCYL3    -0.04207715 -0.0735641493 -0.022479278
# C1orf112  0.02379281 -0.0653475875 -0.016298842
# FGR      -0.03549283 -0.0199220785 -0.019922079


# 最后把两个相关性矩阵联合起来，指定感兴趣模块进行分析
module = "saddlebrown"
pheno = "auc"
modNames = substring(names(MEs), 3)

# 获取关注的列
module_column = match(module, modNames)
pheno_column = match(pheno, colnames(drug_screen_parp))

# 获取模块内的基因
moduleGenes <- moduleColors == module

# 可视化基因与模块、表型的相关性，绘制散点图（想看所有的咱可以批量作图，咦，后续是不是可以分享一下！）
par(mar = c(6, 6, 3, 3))
verboseScatterplot(abs(geneModuleMembership[moduleGenes, module_column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = paste("Gene significance for AUC"),
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

