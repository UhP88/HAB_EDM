
# Script associated with manuscript: [A Prototype Coupled Modelling Approach for Predicting Harmful Algal Blooms: A Case Study in Chile]

## library
import pandas as pd
import pyEDM
import pandas as pd
import numpy as np
import pyEDM
import itertools
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.colors
import seaborn as sns
import pkg_resources
import sys


## loads the dataset
df = pd.read_csv('./BD_FITOSenoReloncaviMetriReady.csv', encoding='cp932')
area= "Metri"
sp="Pseudo-Nitzschia;seriata"

# simplify the dataset
df=pd.DataFrame(df)
df.reset_index(drop=True,inplace=True)
#print(df)

#Calculating best E and best theta for species 1
p = df.columns.get_loc(sp) 
# optimal embedding dimension (ver. rho)
library_string = '10 {}'.format(len(df))
# calculation for speices1
embed_sp1 = pyEDM.EmbedDimension(dataFrame=df, columns=df.columns[p], maxE=10, Tp=1, lib=library_string, pred=library_string, numThreads=10, showPlot=False)
best_embed_sp1_res = embed_sp1.loc[[embed_sp1['rho'].idxmax()]]
best_embed_sp1 = best_embed_sp1_res.iat[0, 0]

theta_sp1 = pyEDM.PredictNonlinear(dataFrame=df, lib=library_string, pred=library_string, columns=df.columns[p], E=int(best_embed_sp1), showPlot=False)
best_theta_sp1_res = theta_sp1.loc[[theta_sp1['rho'].idxmax()]]
best_theta_sp1 = best_theta_sp1_res.iat[0, 0]

#Creating embedding table

target =  'Pseudo-Nitzschia;seriata'
causes = ['Actinoptychus;spp.']
E=int(best_embed_sp1)
embed = pyEDM.Embed(dataFrame=df, E=E, columns=target)
embed.insert(loc=0,column='Time', value=df["index"])
for cause in causes:
    embed[cause] = df[cause].values
embed

# Compute for correlation coefficient for increasing Tp (timestep window)
from scipy import stats
corr_coeff_TEST = []
pvalue_TEST = []
rho = []
for Tp in np.arange(0,5):
    smap_out = pyEDM.SMap( dataFrame = embed, lib = library_string, pred = library_string,  theta = int(best_theta_sp1), E = int(best_embed_sp1), Tp=Tp, embedded = True, showPlot = False,verbose = False, columns=' '.join(embed.columns[1:]),target=target+'(t-0)')
    x=smap_out['predictions']
    x=x.dropna()
    OBS=x.Observations
    PRED=x.Predictions
    r,p = stats.pearsonr(OBS, PRED) 
    corr_coeff_TEST.append(r)
    pvalue_TEST.append(p)
    RHO=pyEDM.ComputeError(smap_out['predictions']['Observations'], smap_out['predictions']['Predictions'] )
    rho.append(RHO)

print({"Pearson_Correlation_Coefficient for each TP":corr_coeff_TEST})
print({"Pvalue for each TP":pvalue_TEST})
print({"ComputeErrorOutput for each TP":rho})

#Forecasting with SMap

library_string = '1 {}'.format(len(embed)-(E-1))
results = pyEDM.SMap( dataFrame=embed, \
     lib = library_string, pred = library_string, \
     theta = int(best_theta_sp1), E = int(best_embed_sp1), Tp = 1, embedded=True, \
     verbose = False, showPlot = False, columns = ' '.join(embed.columns[1:]), target=target+'(t-0)')

predictions_df = results['predictions']
predictions_df['Time'] = pd.to_datetime(predictions_df.Time)
a=pd.DataFrame({'Time': pd.date_range(start=predictions_df.Time.iloc[-3], periods=3, freq='14D',closed='right')})
b = a['Time'][0] 
# selecting new value
c = a['Time'][1]
predictions_df.loc[[len(predictions_df.index)-2],['Time']] = [b]
predictions_df.loc[[len(predictions_df.index)-1],['Time']] = [c]
predictions_df['Time'] = predictions_df['Time'].dt.strftime('%Y-%m-%d')
predictions_df.dtypes

#plotting
plt.figure(figsize=(15, 4))
plt.plot(predictions_df.Time, predictions_df.Observations, label='Observations')
plt.plot(predictions_df.Predictions, label='Predictions')
plt.xlabel("Time")
#plt.xticks(rotation = 90)
plt.xticks(np.arange(0, len(predictions_df.Time)+1, 8),rotation = 90)
plt.ylabel("Scaled abundance of plankton cells/mL")
plt.legend()
plt.title(f'SMAP-Prediction for {area}-{sp}')
plt.savefig('Area_{}_sp_{}_ForecastGraph.png'.format(area, sp),bbox_inches="tight",dpi=300)
plt.show(block=False)
plt.clf()
plt.close()

