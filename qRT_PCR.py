import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import statsmodels.api as sm
from scipy import stats
from scipy.stats import f_oneway
from statsmodels.stats.multicomp import pairwise_tukeyhsd

file='Example_qRT_PCR.csv'

title=file.rstrip('.csv')

testdf=pd.read_csv(file)
groups=testdf.Group.unique()
group1=groups[0]
group2=groups[1]



    
meancontrol=(testdf['dCT'].where(testdf['Group']==group1))
meancontrol = [meancontrol_i for meancontrol_i in meancontrol if str(meancontrol_i) != 'nan']
meancontrol

controlsem=stats.sem(meancontrol)

meancontrol=sum(meancontrol)/len(meancontrol)
meancontrol

testdf['Power']=2 ** -(testdf['dCT']-meancontrol)

i=len(groups)
experimental_rqs=[]

for x in range(1,i):
    group=groups[x]
    experimental=(testdf['dCT'].where(testdf['Group']==group))
    experimental = [experimental_i for experimental_i in experimental if str(experimental_i) != 'nan']

    experimentalsem=stats.sem(experimental)

    experimental=sum(experimental)/len(experimental)
    
    experimental_ddCT=experimental-meancontrol
    control_ddCT=meancontrol-meancontrol
    control_RQ=2 ** -control_ddCT
    experimental_RQ=2** -experimental_ddCT
    
    experimental_rqs.append(experimental_RQ)
    
RQs=[control_RQ]+experimental_rqs





def create_sems():
    i=len(groups)
    experimental_rqs=[]
    experimental_dCTs=[]
    experimental_sems=[]

    for x in range(1,i):
        group=groups[x]
        experimental_dCTs=(testdf['dCT'].where(testdf['Group']==group)).values.tolist()
        experimental_dCTs = [experimental_dCTs_i for experimental_dCTs_i in experimental_dCTs if str(experimental_dCTs_i) != 'nan']
        experimental_RQs=(testdf['Power'].where(testdf['Group']==group)).values.tolist()
        experimental_RQs = [experimental_RQs_i for experimental_RQs_i in experimental_RQs if str(experimental_RQs_i) != 'nan']
        experimental_RQs_sem=stats.sem(experimental_RQs)
        experimental_sems.append(experimental_RQs_sem)
        
    return experimental_sems


def two_sample_test():
    i=len(groups)
    dCTs_ttest=[]
        
    for x in range(0,i):
        group=groups[x]
        dCTs=(testdf['dCT'].where(testdf['Group']==group)).values.tolist()
        dCTs = [dCTs_i for dCTs_i in dCTs if str(dCTs_i) != 'nan']
        dCTs_ttest.append(dCTs)
    
    
    ttest=round(stats.ttest_ind(dCTs_ttest[0],dCTs_ttest[1])[1],10)
    
    d = {'group1': group1, 'group2': group2, 'p-adj': ttest}
    tukeydf = pd.DataFrame(data=d,index=[1])
    return tukeydf


def multiple_sample_test():
    
    
    tukey=pairwise_tukeyhsd(testdf['dCT'],testdf['Group'],0.05)
    tukeydf = pd.DataFrame(data=tukey._results_table.data[1:], columns=tukey._results_table.data[0])
    return tukeydf

def assign_stars():
    if len(groups)>2:
        tukeydf['SE']=(tukeydf['upper']-tukeydf['lower'])/(2*1.96)
        tukeydf['z']=(tukeydf['meandiff']/tukeydf['SE'])
        tukeydf['p_actual']=np.exp((-0.717*tukeydf['z'])-0.416*(tukeydf['z']**2))

        tukeydf_fig=tukeydf.where(tukeydf['p_actual']<0.05)
        tukeydf_fig=tukeydf_fig.dropna()
        
    else:
        tukeydf['p_actual']=tukeydf['p-adj']
        tukeydf_fig=tukeydf.where(tukeydf['p_actual']<0.05)
        tukeydf_fig=tukeydf_fig.dropna()



    def func(x):
        if x < 0.05 and x >= 0.01:
            return "*"
        elif x < 0.01 and x >= 0.001:
            return "**"
        elif x < 0.001 and x >= 0.0001:
            return "***"
        else:
            return "****"


    tukeydf_fig['Stars'] = tukeydf_fig['p_actual'].apply(func)    
    
    i=0
    group_dict={}

    for item in groups:
        x={item: i}
        group_dict.update(x)
        i+=1
        
    tukey_test= tukeydf_fig
    tukey_test=tukey_test.replace({"group1": group_dict})
    tukey_test=tukey_test.replace({"group2": group_dict})
    tukey_test=tukey_test[['group1','group2','Stars']]
    tukey_test['ind']=[np.arange(len(groups))]* len(tukey_test)
    tukey_test['menMeans']=[RQs]* len(tukey_test)
    tukey_test
    
  
    return tukey_test

    


sems=[controlsem]+create_sems()[:len(group)]  

if len(groups)==2:
    tukeydf=two_sample_test()
    
if len(groups)>2:
    tukeydf=multiple_sample_test()
    
starsdf=assign_stars()





fig = plt.figure(figsize=(8, 7))
plt.rcParams.update({'font.size': 15})


#color='full'
color='blackandwhite'





if color=='full':
    #sns.barplot(x=groups,y=RQs,color='lightgrey',edgecolor='black')
    sns.swarmplot(x='Group',y='Power',data=testdf,size=8/np.log(len(groups)))
    
else:
    #sns.barplot(x=groups,y=RQs,color='white',edgecolor='black')
    sns.swarmplot(x='Group',y='Power',data=testdf,size=8/np.log(len(groups)),color='white',edgecolor='black',linewidth=1)
    
    
    
    
plt.errorbar(groups, RQs, yerr=sems, fmt=' ',color='black', zorder=-1, capsize=12,marker='_',linewidth=2)
plt.ylabel('RQ')
plt.xlabel('')
#ylim_max=(max(RQs))+(len(starsdf))
ylim_max=max(testdf['Power'])+(max(testdf['Power'])/10)*(len(starsdf))
ylim_min=0
plt.ylim(ylim_min,ylim_max)




def label_diff(i,j,text,X,Y):
    x = ((X[i]+X[j])/2)-(len(text)/50)
    differential= (abs(Y[i]-Y[j])/(i+j))
    y =  max([Y[i],Y[j]])+differential
    dx = abs(X[i]-X[j])

    props = {'connectionstyle':'arc','arrowstyle':'-',\
                 'shrinkA':10,'shrinkB':10,'linewidth':1}
    plt.annotate(text, xy=(x,y+(ylim_max/8)), zorder=10, size=15, bbox=dict(boxstyle='square,pad=-.1',facecolor='white', edgecolor='none', alpha=1))
    plt.annotate('', xy=(X[i],y+(ylim_max/8)), xytext=(X[j],y+(ylim_max/8)), arrowprops=props)



for n in range (0,len(starsdf)):
    label_diff(starsdf.iloc[n].values.tolist()[0],
               starsdf.iloc[n].values.tolist()[1],
               starsdf.iloc[n].values.tolist()[2],
               starsdf.iloc[n].values.tolist()[3],
               starsdf.iloc[n].values.tolist()[4])


plt.savefig(title+'_arkitect_qRTPCR.png')
