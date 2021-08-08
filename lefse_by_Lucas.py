#!/usr/bin/env python
# -*- coding: UTF-8 -*-
#Editor=lucas
import argparse
import sys
import os

from sklearn.discriminant_analysis import LinearDiscriminantAnalysis as LDA
import pandas as pd
import numpy as np
import scipy.stats as stats
import math
import random
import copy
random.seed(1982)

parser = argparse.ArgumentParser(prog='abundance.py',

    epilog='''
python lefse_by_Lucas.py -s sptab.txt -g group.txt  -o out_dir
''', formatter_class=argparse.RawDescriptionHelpFormatter)


parser.add_argument('-s', '--sptab', required=True, type=str,
                    help='sptab.txt')
parser.add_argument('-g', '--group', required=True, type=str,
                    help='group.txt')  
parser.add_argument('-l', '--lda_th', required=False, type=float, default=2.0,
                    help='default=2.0') 
                    
parser.add_argument('-o', '--output', required=True, type=str,
                    help='out_dir')

def main():
    argsshell = parser.parse_args()
    data=pd.read_csv(argsshell.sptab,index_col=0)
    meta_data=pd.read_csv(argsshell.group,index_col=0)
    #设置kruskal 检验阈值
    k_alpha=0.05
    #设置循环次数
    bootstrap_num=30
    #设置抽样比例
    fract_sample=0.6666
    #设置LDA的阈值
    lda_th=argsshell.lda_th
    
    GGG=set(meta_data.iloc[:,0])
    group={}
    for i in range(len(GGG)):
        group[list(GGG)[i]]=meta_data[meta_data.iloc[:,0]==list(GGG)[i]].index.tolist()





    args=[i for i in group.values()]





    #获取所有分组做kuskal检验的P值


    DF_a=[]
    for x in range(len(group.values())):
        DF_a.append(data.loc[:,args[x]])  


    kruskal_p=[]
    for i in range(len(data.index)):
        test_arg=[]
        for j in range(len(DF_a)):
            test_arg.append(DF_a[j].iloc[i,:].tolist())
        kruskal_p.append(stats.kruskal(*test_arg)[1])



    #合并为表格
    DF1=pd.DataFrame({"index":data.index.tolist(),"kruskal_p":kruskal_p})



    #筛选出来p值小于阈值的特征(种)
    DF_filter=DF1[DF1["kruskal_p"]<k_alpha]



    #获取经过筛选后的有显著差异的种的data_
    data_=data.loc[DF_filter["index"]]




    def contast_within_classes_or_few_per_class(rand_s):
        cls_cols=data_2.iloc[:,rand_s]
        if len(set(cls_cols.columns))<ncl:
            return True
        for c in set(cls_cols.columns):
            if list(cls_cols.columns).count(c)<min_cl:
                return True
            for i in cls_cols.index:
                if (len(set(cls_cols.loc[i])) <= min_cl and min_cl > 1) or (min_cl == 1 and len(set(cls_cols.loc[i])) <= 1):
                    return True    
        return False



    #X为特征矩阵，Y为分组信息
    def lda_score(X,Y):
        global unit
        global mmm
        global cla
        lda = LDA(n_components=1,tol=0.0000000001) 
        lda=lda.fit(X,Y)
        u=lda.scalings_[:,0]
        unit=u/np.sqrt(sum(u*u))
        LD=np.dot(X,unit)
        mmm=lda.means_
        cla=lda.classes_
        return LD



    ##先判断每个组别下的不重复的丰度的个数是否大于4或者该组别的一半（改组中样本数大于8），
    #若不满足这个条件，则对该组别下的每一个数加上 一个 服从均值为0，标准差为(在（该数*0.05和 0.01这个值）中的最大值) 正态分布 的值 后 再取绝对值




    #将筛选后的每个组别下的Dataframe提取出来放到DF_b中。
    DF_b=[]
    for x in range(len(group.values())):
        DF_b.append(data_.loc[:,args[x]])




    #深拷贝DF_b
    DF_b_=copy.deepcopy(DF_b)



    for j in range(len(DF_b_)):
        for i in DF_b_[j].index:
            feats=set(DF_b_[j].loc[i].tolist())
            if len(feats) <= max(len(DF_b_[j].columns)*0.5,4):
                DF_b_[j].loc[i]=[np.abs(x+random.normalvariate(mu=0,sigma=max(x*0.05,0.01))) for x in DF_b_[j].loc[i]]




    #将拆开的DF_b合并为一整个dataframe：data_2，它的X轴为物种，Y轴为样本
    data_2=pd.concat(DF_b_,axis=1)



    ncl=len(set(group.keys())) #ncl 是有几个大的分组

    min_cl=int(float(min([list(group.keys()).count(c) for c in set(group.keys())]))*fract_sample*fract_sample*0.5)
    # min_cl 是每个组别中含量最少的样本数*fract_sample*fract_sample*0.5

    min_cl= max(min_cl,1)
    #如果比1 小 则 min_cl 取1




    rfk=int(float(len(data_2.columns))*fract_sample)



    KKKK=[]
    for k in range(bootstrap_num):
        for rtmp in range(1000):
                rand_s = [random.randint(0,len(data_2.columns)-1) for v in range(rfk)]  #从[0,样本数)中抽取0.6666*样本数个随机数构成列表rand_s
                if not contast_within_classes_or_few_per_class(rand_s=rand_s):
                        break
        #rand_s为抽取的样本的索引

        #得到抽样后的待LDA的列表
        data_3=data_2.iloc[:,rand_s]


        YYY=[]
        for i in data_3.columns:
            for x,y in group.items():
                if i in y:
                    YYY.append(x)

        RRR=[]
        for i in data_3.columns:
            for x,y in group.items():
                if i in y:
                    RRR.append([x,i])

        RRR=np.array(RRR).T

        RRR=pd.DataFrame(RRR,columns=RRR[0]).iloc[1,:]

        #准备输入的数据集的自变量和因变量
        X=np.array(data_3.T)
        Y=np.array(YYY).reshape(1,-1)[0]


        LD=lda_score(X,Y)

        LD=pd.DataFrame(LD).T

        LD.columns=[data_3.columns]

        #得到的DF_b_为按照Lefse作者预处理后的结果

        #得到分组的两两组合
        pairs=[(a,b) for a in set(group.keys()) for b in set(group.keys()) if a>b]

        #得到每个两两组别之间的effectsize 是CCC
        CCC=[]
        for i in range(len(pairs)):
            if (sum(RRR.index==pairs[i][0])*sum(RRR.index==pairs[i][1]))!=0:
                A=np.mean(LD[set(RRR[[pairs[i][0]]])].iloc[0,:].tolist())
                B=np.mean(LD[set(RRR[[pairs[i][1]]])].iloc[0,:].tolist())
                CCC.append(np.abs(A-B))
            else:
                CCC.append(np.nan)


        #计算每两个组别间的wfinal 为 unit * 对应组别的effectsize
        wfinal=[i*unit for i in CCC]

        #对刚刚求得wfinal的值进行判断，若为空值则赋值为0，否则转化为浮点数并取绝对值得到coeff


        coeff=[]
        for i in range(len(wfinal)):
            blank=[]
            for j in wfinal[i]:        
                if math.isnan(j):
                    j=0.0
                else:
                    j=np.abs(j)
                blank.append(j)
            coeff.append(np.array(blank))   

        rres=pd.DataFrame(mmm)

        rres.index=cla

        #这个是从rres中取出每一个特征在该组别下的值 当pp在means的行名中，否则返回[0]*列数                                                                 
                        #对当前的两两组别对比求循环，也就是每次取出一个组别
                    #所以得到的字典res应为keys是两两组别的名称 每一个所对应的values是均值项该组别下的值的列表,该列表下每一个值是一个特征

        res=[]
        for x in range(len(pairs)):
            for i in pairs[x]:
                if i in cla:
                    res.append((i,[float(ff) for ff in rres.loc[i]]))
                else:
                    res.append((i,[0.0]*mmm.shape[1]))

        #i 代表组别 目前有六组
        #j 为特征数

        res_2=[]
        for i in range(len(pairs)):
            res_tmp=[]
            for j,n in enumerate(data_.index):
                gm=np.abs(dict(res)[pairs[i][0]][j]-dict(res)[pairs[i][1]][j])
                res_tmp.append((n,(gm+coeff[i][j])*0.5))
            res_2.append(res_tmp)
        KKKK.append(res_2)




    #KKKK是包含多次重复抽样的结果的列表



    res_3={}
    for q in data_3.index:
        qwe=[]
        for p in range(len(pairs)):        
            qwe.append(np.mean([dict(KKKK[kk][p])[q] for kk in range(bootstrap_num)]))
        m=max(qwe)
        res_3[q]=math.copysign(1.0,m)*np.log10(1.0+np.abs(m))
        
    # res_3是所有特征的LDA score




    #res_th 是符合满足大于阈值的LDA值和它所对应的特征名
    res_th=[(k,x) for k,x in res_3.items() if np.abs(x) >lda_th]




    #找出符合lda阈值筛选条件下，每个特征所在的最大平均丰度的组别名


    cls_mean=[DF_b[i].mean(axis=1) for i in range(len(DF_b))]
    #cls_mean是列表包含每个组别的不同均值的dataframe



    cls_mean_=pd.concat(cls_mean,axis=1)
    cls_mean_.columns=list(group.keys())
    #cls_mean_为列名为组别名，每一行为不同特征所对应的均值


    after_lda_th_mean_abundance=cls_mean_.loc[dict(res_th).keys()]




    des=[]
    for i,j in enumerate(after_lda_th_mean_abundance.max(axis=1)):
        des.append(after_lda_th_mean_abundance.iloc[i][j==after_lda_th_mean_abundance.iloc[i]].index.tolist()[0])




    final_score=pd.DataFrame(res_th,columns=["feature","LDA_score"])



    final_score["group"]=des


    final_score.to_csv(os.path.join(argsshell.output),index=None)
    print('OK!')
   

if __name__ == "__main__":
    main()
