
# Script associated with manuscript: [A Prototype Coupled Modelling Approach for Predicting Harmful Algal Blooms: A Case Study in Chile]

## library
import pandas as pd
import numpy as np
import pyEDM
import itertools
import matplotlib.pyplot as plt
import seaborn as sns

## loads the dataset
#Contains scaled phytoplankton abundance data
df = pd.read_csv('./BD_FITOSenoReloncaviMetriReady.csv', encoding='cp932')

# define area to save plots
areafito ="Metri"

## CCM with optimal embedding dimension for rho
# Extract the variables
#Change target species name as needed; Here is an example for Pseudo-Nitzschia;seriata
target_list = df.columns[0:]
n = len(target_list)
i = df.columns.get_loc("Pseudo-Nitzschia;seriata")
# run loop

for j in range(1, n):
    if j == i: #skips to the next iteration
        continue

    # define the focus species
    sp1 = target_list[i]
    sp2 = target_list[j]

    # optimal embedding dimension (ver. rho)
    library_string = '10 {}'.format(len(df))
    # calculation for speices1
    embed_sp1 = pyEDM.EmbedDimension(dataFrame=df, columns=df.columns[i], maxE=10, Tp=1, lib=library_string, pred=library_string, numThreads=10, showPlot=False)
    best_embed_sp1_res = embed_sp1.loc[[embed_sp1['rho'].idxmax()]]
    best_embed_sp1 = best_embed_sp1_res.iat[0, 0]
    # calculation for speices2
    embed_sp2 = pyEDM.EmbedDimension(dataFrame=df, columns=df.columns[j], maxE=10, Tp=1, lib=library_string, pred=library_string, numThreads=10, showPlot=False)
    best_embed_sp2_res = embed_sp2.loc[[embed_sp2['rho'].idxmax()]]
    best_embed_sp2 = best_embed_sp2_res.iat[0, 0]
    
    # CCM from sp1 to sp2
    '''
    <parameter>
    E=int(best_embed_sp1)
    Tp=0
    libSize="START FINISH STEP"
    sample=100
    seed=12345
    '''
    ccm_sp1 = pyEDM.CCM(dataFrame=df[['index', sp1, sp2]], E=int(best_embed_sp1), Tp=0, columns=sp1, target=sp2, libSizes="10 210 20", sample=100, seed=12345, showPlot=False)
    plt.plot(ccm_sp1.iloc[:, 0], ccm_sp1.iloc[:, 1], marker='o', color="red", label=ccm_sp1.columns[1])
    plt.plot(ccm_sp1.iloc[:, 0], ccm_sp1.iloc[:, 2], marker='o', color="blue", label=ccm_sp1.columns[2])
    plt.hlines(y=0, xmin=10, xmax=210, color="black", linestyles='dotted')
    plt.xlabel('Library size')
    plt.ylabel('Correlation ' + r'$\rho$')
    plt.title('From_{}_To_{}_EmbedSp1_{}'.format(sp1, sp2, int(best_embed_sp1)))
    plt.legend()
    plt.savefig('{}_ccm_sp1_{}_sp2_{}.png'.format(areafito,sp1, sp2))
    plt.show(block=False)
    plt.clf()
    plt.close()
    
    # CCM from sp2 to sp1
    '''
    <parameter>
    E=int(best_embed_sp2)
    Tp=0
    libSize="START FINISH STEP"
    sample=100
    seed=12345
    '''
    ccm_sp2 = pyEDM.CCM(dataFrame=df[['index', sp2, sp1]], E=int(best_embed_sp2), Tp=0, columns=sp2, target=sp1, libSizes="10 210 20", sample=100, seed=12345, showPlot=False)
    plt.plot(ccm_sp2.iloc[:, 0], ccm_sp2.iloc[:, 1], marker='o', color="red", label=ccm_sp2.columns[1])
    plt.plot(ccm_sp2.iloc[:, 0], ccm_sp2.iloc[:, 2], marker='o', color="blue", label=ccm_sp2.columns[2])
    plt.hlines(y=0, xmin=10, xmax=210, color="black", linestyles='dotted')
    plt.xlabel('Library size')
    plt.ylabel('Correlation ' + r'$\rho$')
    plt.title('From_{}_To_{}_EmbedSp2_{}'.format(sp2, sp1, int(best_embed_sp2)))
    plt.legend()
    plt.savefig('{}_ccm_sp2_{}_sp1_{}.png'.format(areafito,sp2, sp1))
    plt.show(block=False)
    plt.clf()
    plt.close()

## S-map with optimal embedding dimension for rho
# extract the variables
target_list = df.columns[0:]
n = len(target_list)
i = df.columns.get_loc("Pseudo-Nitzschia;seriata")
# run loop
for j in range(1, n):
    if j==i: #skips to the next iteration
        continue	


    # define the focus species
    sp1 = target_list[i]
    sp2 = target_list[j]
    
    # optimal embedding dimension (ver. rho)
    library_string = '10 {}'.format(len(df))
    # calculation for speices1
    embed_sp1 = pyEDM.EmbedDimension(dataFrame=df, columns=df.columns[i], maxE=10, Tp=1, lib=library_string, pred=library_string, numThreads=10, showPlot=False)
    best_embed_sp1_res = embed_sp1.loc[[embed_sp1['rho'].idxmax()]]
    best_embed_sp1 = best_embed_sp1_res.iat[0, 0]
    # calculation for speices2
    embed_sp2 = pyEDM.EmbedDimension(dataFrame=df, columns=df.columns[j], maxE=10, Tp=1, lib=library_string, pred=library_string, numThreads=10, showPlot=False)
    best_embed_sp2_res = embed_sp2.loc[[embed_sp2['rho'].idxmax()]]
    best_embed_sp2 = best_embed_sp2_res.iat[0, 0]

    # non-linearity test (ver. rho)
    # calculation for speices1
    theta_sp1 = pyEDM.PredictNonlinear(dataFrame=df, lib=library_string, pred=library_string, columns=df.columns[i], E=int(best_embed_sp1), showPlot=False)
    best_theta_sp1_res = theta_sp1.loc[[theta_sp1['rho'].idxmax()]]
    best_theta_sp1 = best_theta_sp1_res.iat[0, 0]
    # calculation for speices2
    theta_sp2 = pyEDM.PredictNonlinear(dataFrame=df, lib=library_string, pred=library_string, columns=df.columns[j], E=int(best_embed_sp2), showPlot=False)
    best_theta_sp2_res = theta_sp2.loc[[theta_sp2['rho'].idxmax()]]
    best_theta_sp2 = best_theta_sp2_res.iat[0, 0]

    # S-map sp1->sp2
    '''
    <parameter>
    lib="START FINISH"
    pred='START FINISH'
    E=int(best_embed_sp1)
    theta=int(best_theta_sp1)
    seed=12345
    '''
    smap_sp1 = pyEDM.SMap(dataFrame=df, lib="10 210", pred="10 210", columns=sp1, target=sp2, E=int(best_embed_sp1), theta=int(best_theta_sp1), showPlot=False)
    coef_sp1 = smap_sp1['coefficients']
    pred_sp1 = pyEDM.ComputeError(smap_sp1['predictions']['Observations'],smap_sp1['predictions']['Predictions'] )
    rho_sp1 = pred_sp1['rho']
    interaction_stg_sp1 = coef_sp1[coef_sp1.columns[-1]].values
    plt.plot(coef_sp1['index'], interaction_stg_sp1, color="blue", label=f'Interaction strength from {sp1} to {sp2}')
    plt.axhline(y=0, linestyle='dotted', color='black')
    plt.axhline(y=np.nanmean(interaction_stg_sp1), linestyle='dotted', color='red', label='average interaction strength')
    plt.xlabel('Date')
    plt.ylabel('Strength')
    plt.xticks([7/19/1991, 7/19/1993])
    plt.legend()
    plt.savefig('{}_smap_sp1_{}_sp2_{}_rho{}.png'.format(areafito,sp1, sp2,rho_sp1))
    plt.show(block=False)
    plt.clf()
    plt.close()
    
    # S-map sp2->sp1
    '''
    <parameter>
    lib="START FINISH"
    pred='START FINISH'
    E=int(best_embed_sp2)
    theta=int(best_theta_sp2)
    seed=12345
    '''
    smap_sp2 = pyEDM.SMap(dataFrame=df, lib="10 210", pred="10 210", columns=sp2, target=sp1, E=int(best_embed_sp2), theta=int(best_theta_sp2), showPlot=False)
    coef_sp2 = smap_sp2['coefficients']
    pred_sp2 = pyEDM.ComputeError(smap_sp2['predictions']['Observations'],smap_sp2['predictions']['Predictions'] )
    rho_sp2 = pred_sp2['rho']
    interaction_stg_sp2 = coef_sp2[coef_sp2.columns[-1]].values
    plt.plot(coef_sp2['index'], interaction_stg_sp2, color="blue", label=f'Interaction strength from {sp2} to {sp1}')
    plt.axhline(y=0, linestyle='dotted', color='black')
    plt.axhline(y=np.nanmean(interaction_stg_sp2), linestyle='dotted', color='red', label='average interaction strength')
    plt.xlabel('Date')
    plt.ylabel('Strength')
    plt.xticks([7/19/1991, 7/19/1993])
    plt.legend()
    plt.savefig('{}_smap_sp2_{}_sp1_{}_rho{}.png'.format(areafito,sp2, sp1,rho_sp2))
    plt.show(block=False)
    plt.clf()
    plt.close()

## function of making twin-surrogate data
from scipy.spatial.distance import pdist, squareform

def twin_surrogate(df, target='x5', e_dim=2, max_iter=100, \
                   method='original', obs_per_year=None):
    # 時間遅れ埋め込みの状態空間を作る。
    x_e = pyEDM.Embed(dataFrame=df, E=e_dim, columns=target)
    # all vs. allの距離行列を計算する。max norm.
    dist_mtx = squareform(pdist(x_e.values, metric='chebyshev'))
    # 距離行列を、距離が近い・遠いでゼロイチのマトリックスにしてしまう。どのあたりで切るかは適当...
    binary_dist_mtx = (dist_mtx > np.percentile(dist_mtx, 75.0)).astype(int)
    
    # all vs. allのペアが "twin" であるかどうかのマトリックス
    # ほかの全部との距離パターンが似ている2点は twin.
    twins_table = np.eye(len(x_e), len(x_e))
    for i, j in itertools.combinations(range(len(x_e)), 2):
        # ゼロイチパターンがすべて一致
        if np.all(binary_dist_mtx[:, i] == binary_dist_mtx[:, j]):
            if method == 'phase_lock':
                if (j-i) % obs_per_year == 0:
                    twins_table[i, j] = 1
                    twins_table[j, i] = 1
            else:
                twins_table[i, j] = 1
                twins_table[j, i] = 1

    n_iter = 0
    while True:
        failed_surrogate_generation = False
        surrogate = np.zeros((x_e.values.shape[0] + (e_dim - 1), x_e.values.shape[1]))
        # 適当な点からはじめる
        if method == 'phase_lock':
            ind = np.random.choice(np.arange(0, len(x_e), obs_per_year), 1)[0] # random sampling from same season samples
        else:
            ind = np.random.choice(len(x_e), 1)[0]
        surrogate[0, :] = x_e.values[ind, :]
        for t in np.arange(surrogate.shape[0])[1:]:
            # 次の点はオリジナル系列の t+1 か、あるいは、twin の中からランダムに選択された点
            nex = np.random.choice(np.argwhere(twins_table[ind, :]).ravel(), 1)[0] + 1
            # もし途中で系列の終端に達してしまったら全体をはじめからやり直す
            if nex >= len(x_e):
                failed_surrogate_generation = True
                break
            surrogate[t, :] = x_e.values[nex, :]
            ind = nex
        if not failed_surrogate_generation:
            break
        else:
            n_iter += 1
            if n_iter > max_iter:
                # max_iter回繰り返して全部失敗したらあきらめる
                raise Exception(f'Failed to generate surrogate in {max_iter} iterations.')

    return surrogate[:, -1]

## function of doing ccm with surrogate test
from tqdm import tqdm
#from tqdm.notebook import tqdm
#https://github.com/CosmiQ/solaris/issues/392

def ccm_surrogate_ConfidenceInterval(dataframe, E=2, \
        column='column', target='target', libSizes="10 210 20", \
        method='phase_lock', obs_per_year=24, n_surrogate=100, CI=95.0):
    # 普通のCCM計算
    original_rho = pyEDM.CCM(dataFrame = dataframe, \
        E = E, columns = column, target = target, Tp=0, \
        libSizes = libSizes, sample = 100, \
        verbose = True,  showPlot = False) [f'{column}:{target}']
    libSize_indices = original_rho.index
    original_rho = original_rho.values

    # 百本のサロゲート作ってそれぞれでCCM計算
    surrogate_rhos = np.zeros((len(original_rho), n_surrogate))
    for i in tqdm(range(n_surrogate)):
        surrogate = twin_surrogate(dataframe, target=column, \
                                   e_dim=E, method=method, obs_per_year=obs_per_year)
        dataframe['surrogate'] = surrogate[E-1:len(surrogate)]
        surrogate_rho = pyEDM.CCM(dataFrame = dataframe, \
            E = E, columns = 'surrogate', target = target, Tp=0, \
            libSizes = libSizes, sample = 100, verbose = True,  showPlot = False) [f'surrogate:{target}'].values
        surrogate_rhos[:, i] = surrogate_rho

    result_df = pd.DataFrame(np.zeros((len(original_rho), 3)), \
                             columns=['rho', f'upper{int(CI):d}', f'lower{int(CI):d}'], index=libSize_indices)
    result_df['rho'] = original_rho
    result_df[f'upper{int(CI):d}'] = np.nanpercentile(surrogate_rhos, 100.0 - (100.0 - CI) / 2.0, axis=1)[:, np.newaxis]
    result_df[f'lower{int(CI):d}'] = np.nanpercentile(surrogate_rhos, (100.0 - CI) / 2.0, axis=1)[:, np.newaxis]
    return result_df

## CCM with twin-surrogate test
# ignore warnings
import warnings
warnings.simplefilter('ignore')

# extract the variables
target_list = df.columns[0:]
n = len(target_list)
i = df.columns.get_loc("Pseudo-Nitzschia;seriata")
# run loop
for j in range(1, n):
    if j == i: #skips to the next iteration
        continue
    # define the focus species
    sp1 = target_list[i]
    sp2 = target_list[j]
    
    # optimal embedding dimension
    library_string = '10 {}'.format(len(df))
    # calculation species1
    embed_sp1 = pyEDM.EmbedDimension(dataFrame=df, columns=df.columns[i], maxE=10, Tp=1, lib=library_string, pred=library_string, numThreads=10, showPlot=False)
    best_embed_sp1_res = embed_sp1.loc[[embed_sp1['rho'].idxmax()]]
    best_embed_sp1 = best_embed_sp1_res.iat[0, 0]
    # calculation species1
    embed_sp2 = pyEDM.EmbedDimension(dataFrame=df, columns=df.columns[j], maxE=10, Tp=1, lib=library_string, pred=library_string, numThreads=10, showPlot=False)
    best_embed_sp2_res = embed_sp2.loc[[embed_sp2['rho'].idxmax()]]
    best_embed_sp2 = best_embed_sp2_res.iat[0, 0]
    
    # CCM with twin-surrogate
    '''
    <parameter>
    E=int(best_embed_sp1)
    libSizes="START FINISH STEP"
    method="phase_lock" # Or "original"
    obs_per_year=24
    n_surrogate=100
    CI=95.0 # confidence interval
    seed=12345
    '''
    ccm_sp1 = ccm_surrogate_ConfidenceInterval(df[['index', sp1, sp2]], E=int(best_embed_sp1), column=sp1, target=sp2, libSizes="10 210 20", method='phase_lock', obs_per_year=24, n_surrogate=100, CI=95.0)
    ccm_sp2 = ccm_surrogate_ConfidenceInterval(df[['index', sp2, sp1]], E=int(best_embed_sp2), column=sp2, target=sp1, libSizes="10 210 20", method='phase_lock', obs_per_year=24, n_surrogate=100, CI=95.0)
    
    fig, ax = plt.subplots()
    ax.plot(range(0, 210, 20), ccm_sp1['rho'], color='red', label='actual{}:{}'.format(sp1, sp2))
    ax.fill_between(range(0, 210, 20), ccm_sp1['lower95'], ccm_sp1['upper95'], \
        facecolor='magenta', edgecolor='magenta', alpha=0.5, label='surrogate_{}:{}'.format(sp1, sp2))
    ax.plot(range(0,210,20), ccm_sp2['rho'], color='blue', label='actual_{}:{}'.format(sp2, sp1))
    ax.fill_between(range(0, 210, 20), ccm_sp2['lower95'], ccm_sp2['upper95'], \
        facecolor='cyan', edgecolor='cyan', alpha=0.5, label='surrogate_{}:{}'.format(sp2, sp1))
    ax.axhline(y=0.0, linestyle=':', color='black')
    plt.title('Relationship_{}_And_{}'.format(sp1, sp2))
    ax.set_xlabel('Library Size')
    ax.set_ylabel('Correlation ' + r'$\rho$')
    ax.legend()
    plt.savefig('{}_twinsurr_sp1_{}_sp2_{}.png'.format(areafito,sp1, sp2))
    sns.despine()
    plt.show(block=False)
    plt.clf()
    plt.close()

