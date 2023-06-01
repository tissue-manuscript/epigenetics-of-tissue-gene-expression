import csv
import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt


filename = "8tissues_all_male_female_samples_final.csv"
df = pd.read_csv(filename, sep=',')
df = df.fillna(0)
dft = df.transpose()
dft = dft.iloc[1:]

dft.insert(0, "target", [
    'Mus (M)', 'Mus (M)', 'Mus (M)', 'Mus (M)', 'Mus (M)', 'Mus (F)', 'Mus (F)', 'Mus (F)', 'Mus (F)', 'Mus (F)', 
    'Hea (M)', 'Hea (M)', 'Hea (M)', 'Hea (M)', 'Hea (M)', 'Hea (F)', 'Hea (F)', 'Hea (F)', 'Hea (F)', 'Hea (F)', 
    'Liv (M)', 'Liv (M)', 'Liv (M)', 'Liv (M)', 'Liv (M)', 'Liv (F)', 'Liv (F)', 'Liv (F)', 'Liv (F)', 'Liv (F)', 
    'Lun (M)', 'Lun (M)', 'Lun (M)', 'Lun (M)', 'Lun (M)', 'Lun (F)', 'Lun (F)', 'Lun (F)', 'Lun (F)', 'Lun (F)', 
    'Kid (M)', 'Kid (M)', 'Kid (M)', 'Kid (M)', 'Kid (M)', 'Kid (F)', 'Kid (F)', 'Kid (F)', 'Kid (F)', 'Kid (F)', 
    'Hip (M)', 'Hip (M)', 'Hip (M)', 'Hip (M)', 'Hip (M)', 'Hip (F)', 'Hip (F)', 'Hip (F)', 'Hip (F)', 'Hip (F)', 
    'BA (M)', 'BA (M)', 'BA (M)', 'BA (M)', 'BA (M)', 'BA (F)', 'BA (F)', 'BA (F)', 'BA (F)', 'BA (F)', 
    'WA (M)', 'WA (M)', 'WA (M)', 'WA (M)', 'WA (M)', 'WA (F)', 'WA (F)', 'WA (F)', 'WA (F)', 'WA (F)'
])

# dft.insert(0, "target", ["None"] * 30650)

print(dft)


features = list(range(0, 15533))

x = dft.loc[:, features].values# Separating out the target
y = dft.loc[:,['target']].values# Standardizing the features
x = StandardScaler().fit_transform(x)


pca = PCA(n_components=3)
principalComponents = pca.fit_transform(x)
principalDf = pd.DataFrame(data = principalComponents, columns = ['principal component 1', 'principal component 2', 'principal component 3'])

principalDf.insert(3, "target", [
    'Mus (M)', 'Mus (M)', 'Mus (M)', 'Mus (M)', 'Mus (M)', 'Mus (F)', 'Mus (F)', 'Mus (F)', 'Mus (F)', 'Mus (F)', 
    'Hea (M)', 'Hea (M)', 'Hea (M)', 'Hea (M)', 'Hea (M)', 'Hea (F)', 'Hea (F)', 'Hea (F)', 'Hea (F)', 'Hea (F)', 
    'Liv (M)', 'Liv (M)', 'Liv (M)', 'Liv (M)', 'Liv (M)', 'Liv (F)', 'Liv (F)', 'Liv (F)', 'Liv (F)', 'Liv (F)', 
    'Lun (M)', 'Lun (M)', 'Lun (M)', 'Lun (M)', 'Lun (M)', 'Lun (F)', 'Lun (F)', 'Lun (F)', 'Lun (F)', 'Lun (F)', 
    'Kid (M)', 'Kid (M)', 'Kid (M)', 'Kid (M)', 'Kid (M)', 'Kid (F)', 'Kid (F)', 'Kid (F)', 'Kid (F)', 'Kid (F)', 
    'Hip (M)', 'Hip (M)', 'Hip (M)', 'Hip (M)', 'Hip (M)', 'Hip (F)', 'Hip (F)', 'Hip (F)', 'Hip (F)', 'Hip (F)', 
    'BA (M)', 'BA (M)', 'BA (M)', 'BA (M)', 'BA (M)', 'BA (F)', 'BA (F)', 'BA (F)', 'BA (F)', 'BA (F)', 
    'WA (M)', 'WA (M)', 'WA (M)', 'WA (M)', 'WA (M)', 'WA (F)', 'WA (F)', 'WA (F)', 'WA (F)', 'WA (F)'
])
finalDf = principalDf

print(finalDf)
print(pca.explained_variance_ratio_)

maleDf = finalDf[finalDf["target"].str.contains("(M)", na=False)]
femaleDf = finalDf[finalDf["target"].str.contains("(F)", na=False)]

targets = [
    'Mus (M)', 'Mus (F)', 
    'Hea (M)', 'Hea (F)',
    'Liv (M)', 'Liv (F)',
    'Lun (M)', 'Lun (F)',
    'Kid (M)', 'Kid (F)',
    'Hip (M)', 'Hip (F)',
    'BA (M)', 'BA (F)', 
    'WA (M)', 'WA (F)',
]
targets_m = [t for t in targets if t.endswith("(M)")]
targets_f = [t for t in targets if t.endswith("(F)")]

colors = [
    'r', 'cornflowerblue', 'aqua', 'tab:brown', 
    'm', 'y', 'k', 'lime',
    'darkslategrey', 'darkgreen', 'indigo', 'firebrick',
    'darksalmon', 'olivedrab', 'silver', 'yellowgreen'
]


fig = plt.figure(figsize = (8,8))
ax = fig.add_subplot(1,1,1) 
ax.set_xlabel(f'PC1 [variance: {pca.explained_variance_ratio_[0]*100:.2f}%]', fontsize = 15)
ax.set_ylabel(f'PC2 [variance: {pca.explained_variance_ratio_[1]*100:.2f}%]', fontsize = 15)

for target, color in zip(targets_m, colors):
    indicesToKeep = finalDf['target'] == target
    ax.scatter(finalDf.loc[indicesToKeep, 'principal component 1']
               , finalDf.loc[indicesToKeep, 'principal component 2']
               , c = color
               , s = 150
               , marker = ".")

for target, color in zip(targets_f, colors):
    indicesToKeep = finalDf['target'] == target
    ax.scatter(finalDf.loc[indicesToKeep, 'principal component 1']
               , finalDf.loc[indicesToKeep, 'principal component 2']
               , c = color
               , s = 150
               , marker = "+")

ax.legend(targets, bbox_to_anchor=(1.2, 1.05))
ax.grid()

plt.savefig('RNA-8tissues-1-2-MF.png', bbox_inches='tight', dpi=300)


fig = plt.figure(figsize = (8,8))
ax = fig.add_subplot(1,1,1) 
ax.set_xlabel(f'PC1 [variance: {pca.explained_variance_ratio_[0]*100:.2f}%]', fontsize = 15)
ax.set_ylabel(f'PC3 [variance: {pca.explained_variance_ratio_[2]*100:.2f}%]', fontsize = 15)

for target, color in zip(targets_m, colors):
    indicesToKeep = finalDf['target'] == target
    ax.scatter(finalDf.loc[indicesToKeep, 'principal component 1']
               , finalDf.loc[indicesToKeep, 'principal component 3']
               , c = color
               , s = 150
               , marker = ".")

for target, color in zip(targets_f, colors):
    indicesToKeep = finalDf['target'] == target
    ax.scatter(finalDf.loc[indicesToKeep, 'principal component 1']
               , finalDf.loc[indicesToKeep, 'principal component 3']
               , c = color
               , s = 150
               , marker = "+")

ax.legend(targets, bbox_to_anchor=(1.2, 1.05))
ax.grid()

plt.savefig('RNA-8tissues-1-3-MF.png', bbox_inches='tight', dpi=300)



fig = plt.figure(figsize = (8,8))
ax = fig.add_subplot(1,1,1) 
ax.set_xlabel(f'PC2 [variance: {pca.explained_variance_ratio_[1]*100:.2f}%]', fontsize = 15)
ax.set_ylabel(f'PC3 [variance: {pca.explained_variance_ratio_[2]*100:.2f}%]', fontsize = 15)

for target, color in zip(targets_m, colors):
    indicesToKeep = finalDf['target'] == target
    ax.scatter(finalDf.loc[indicesToKeep, 'principal component 2']
               , finalDf.loc[indicesToKeep, 'principal component 3']
               , c = color
               , s = 150
               , marker = ".")

for target, color in zip(targets_f, colors):
    indicesToKeep = finalDf['target'] == target
    ax.scatter(finalDf.loc[indicesToKeep, 'principal component 2']
               , finalDf.loc[indicesToKeep, 'principal component 3']
               , c = color
               , s = 150
               , marker = "+")

ax.legend(targets, bbox_to_anchor=(1.0, 1.05))
ax.grid()

plt.savefig('RNA-8tissues-2-3-MF.png', bbox_inches='tight', dpi=300)
