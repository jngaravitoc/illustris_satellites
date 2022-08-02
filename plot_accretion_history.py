"""
Script to analyze MW analogues with three mergers > 1/10. 

"""


import numpy as np
import matplotlib.pyplot as plt
import illustris_python as il


def plot_accretion_histories(bp, groupids, figname):
    """

    """

    fields = ['SubhaloMass','SubfindID','SnapNum']
    for n in groupids:
      lb_arr = []
      tree = il.sublink.loadTree(bp, 135, n, fields=fields,onlyMPB=True)
      snap_t = tree["SnapNum"]
      for i in snap_t:
        lb_arr.append(lb_time[np.where(snap==i)[0]])
      plt.plot(lb_arr, tree['SubhaloMass'],'-', alpha=0.4, lw=0.5)
      plt.fill_betweenx(np.logspace(-1.5, 3, 10), np.ones(10)*9, np.ones(10)*12, color='C0', alpha=0.3)
      plt.fill_betweenx(np.logspace(-1.5, 3, 10), np.ones(10)*6, np.ones(10)*8, color='C1', alpha=0.3)
      plt.fill_betweenx(np.logspace(-1.5, 3, 10), np.ones(10)*1.5, np.ones(10)*3, color='C2', alpha=0.3)
      plt.yscale('log')
      plt.xlabel('$Lookback   t [Gyr]$')
      plt.ylabel('Total Subhalo Mass [code units]')
      plt.ylim(1, 3E2)
    plt.savefig(figname, bbox_inches='tight')
    plt.close()


def get_analogues(sat1_lbt, sat2_lbt, sat3_lbt, mrat1, mrat2, mrat3, id_halos, analogue, write_data=False):
    """
    By definiton sat1 is the most recent merger, sat2 the intermediate and sat3
    the oldest merger. 

    """
    GSE = 9
    SAG = 6.5
    LMC = 2.0
    Mlmc_min = 0.08
    Mlmc_max = 0.35
    Msag_min = 0.05
    Msag_max = 0.2
    Mgse_min = 0.08
    Mgse_max = 0.4

   
   
    if analogue == "LMC":
        mw_analogues = np.where((sat1_lbt < LMC+1) & (sat1_lbt > LMC-1) & (mrat1<Mlmc_max) & (mrat1>Mlmc_min))
        filename = "LMC_analogues.txt"

    elif analogue == "LMC+SAG":
        mw_analogues = np.where(((sat1_lbt < LMC+1) & (sat1_lbt > LMC-1)) &
            ((sat2_lbt < SAG+1.5) & (sat2_lbt>SAG-1.5)) & (mrat2<Msag_max) &
            (mrat2>Msag_min) & (mrat1<Mlmc_max) & (mrat1>Mlmc_min))
        filename = "LMC_SAG_analogues.txt"

    elif analogue == "LMC+GSE":
        mw_analogues = np.where(((sat1_lbt < LMC+1) & (sat1_lbt > LMC-1)) &
            ((sat3_lbt < GSE+2) & (sat3_lbt > GSE-2)) & (mrat1<Mlmc_max) &  (mrat1>Mlmc_min) &
            (mrat3>Mgse_min) & (mrat3<Mgse_max))
        filename = "LMC_GSE_analogues.txt"
    
    elif analogue == "LMC+SAG+GSE":
        mw_analogues = np.where(((sat1_lbt < LMC+1) & (sat1_lbt > LMC-1)) &
            ((sat2_lbt < SAG+1.5) & (sat2_lbt>SAG-1.5)) & ((sat3_lbt < GSE+2) &
              (sat3_lbt > GSE-2)) & ((mrat2<Msag_max) & (mrat2>Msag_min)) &
            ((mrat1<Mlmc_max) & (mrat1>Mlmc_min)) &((mrat3>Mgse_min) & (mrat3<Mgse_max)))
        filename = "LMC_SAG_GSE_analogues.txt"
    
    elif analogue == "SAG":
        mw_analogues = np.where(((sat2_lbt < SAG+1.5) & (sat2_lbt>SAG-1.5)))
        filename = "SAG_analogues.txt"
    
    elif analogue == "SAG+GSE":
        mw_analogues = np.where(((sat2_lbt < SAG+1.5) & (sat2_lbt>SAG-1.5)) &  ((sat3_lbt < GSE+2) & (sat3_lbt > GSE-2)))
        filename = "SAG_GSE_analogues.txt"
    
    elif analogue == "GSE":
        mw_analogues = np.where(((sat3_lbt < GSE+2) & (sat3_lbt > GSE-2)))
        filename = "GSE_analogues.txt"
       
    N_analogues = len(mw_analogues[0])
    IDhalos = id_halos[mw_analogues]
    lbtimes = [sat1_lbt[mw_analogues], sat2_lbt[mw_analogues], sat3_lbt[mw_analogues]]
    mratios = [mrat1[mw_analogues[0]], mrat2[mw_analogues[0]], mrat3[mw_analogues[0]]]

    if write_data==True :         
        head = "ID host halo, mass ratio sat1, mass ratio sat 2, mass ratio sat 3, looback t sat1, lookback t sat 2, lookback sat 3"
        np.savetxt("flex_def_"+"tophat_"+filename, np.array([IDhalos, mratios[0], mratios[1],
              mratios[2], lbtimes[0], lbtimes[1], lbtimes[2]]).T,
              fmt=['%d','%.4f', '%.4f', "%.4f", "%.4f", "%.4f", "%.4f"],
              header=head )


    return N_analogues, IDhalos, mratios, lbtimes


####### DATA ###########

# Illustris sims path
basePath = '/mnt/home/nico/ceph/illutris/Illustris-3/output'

# MW analogues with 3 mergers 
mw_analogues = np.loadtxt("halos_three_mergers_1to20.txt")
GroupFirstSub = mw_analogues[:,0]

nsnap_sat1 = mw_analogues[:,5]
nsnap_sat2 = mw_analogues[:,6]
nsnap_sat3 = mw_analogues[:,7]
mrat1 = mw_analogues[:,14]
mrat2 = mw_analogues[:,15]
mrat3 = mw_analogues[:,16]

# snap, redshift age lookback time
snap_z = np.loadtxt("snaps_z_age.txt")
snap = snap_z[:,0]
z = snap_z[:,1]
age = snap_z[:,2]
lb_time = snap_z[:,3]


# Defining lookback times of 3 mergers in the MW

#GSE = 10
#SAG = 7
#LMC = 2.2


# Getting lookback times for satellites and sorting them 

sat1lbt=np.zeros_like(nsnap_sat1)
sat2lbt=np.zeros_like(nsnap_sat1)
sat3lbt=np.zeros_like(nsnap_sat1)

sat1_mrat=np.zeros_like(nsnap_sat1)
sat2_mrat=np.zeros_like(nsnap_sat1)
sat3_mrat=np.zeros_like(nsnap_sat1)

## Sort mergers as a function of loolback time

# Sat1 -> most recent merger
# Sat3 -> oldest merger

for j in range(len(nsnap_sat1)):
  s1 = lb_time[np.where(snap==nsnap_sat1[j])[0]]
  s2 = lb_time[np.where(snap==nsnap_sat2[j])[0]]
  s3 = lb_time[np.where(snap==nsnap_sat3[j])[0]]
  
  all_s = np.array([s1[0], s2[0], s3[0]])
  sat1lbt[j] = np.sort(all_s)[0]
  sat2lbt[j] = np.sort(all_s)[1]
  sat3lbt[j] = np.sort(all_s)[2]
  
  assert sat1lbt[j]<=sat2lbt[j]<=sat3lbt[j], "! Error lookback times are not ordered "
  all_mrat = np.array([mrat1[j], mrat2[j], mrat3[j]])
  sat1_mrat[j] = all_mrat[np.argsort(all_s)[0]]
  sat2_mrat[j] = all_mrat[np.argsort(all_s)[1]]
  sat3_mrat[j] = all_mrat[np.argsort(all_s)[2]]
######## ANALOGUES  ########

LMC_analogues = get_analogues(sat1lbt, sat2lbt, sat3lbt, sat1_mrat, sat2_mrat,  sat3_mrat, GroupFirstSub, "LMC", write_data=True)
LMC_SAG_analogues = get_analogues(sat1lbt, sat2lbt, sat3lbt, sat1_mrat, sat2_mrat, sat3_mrat, GroupFirstSub, "LMC+SAG", write_data=True)
LMC_GSE_analogues = get_analogues(sat1lbt, sat2lbt, sat3lbt, sat1_mrat, sat2_mrat, sat3_mrat, GroupFirstSub, "LMC+GSE")
LMC_SAG_GSE_analogues = get_analogues(sat1lbt, sat2lbt, sat3lbt, sat1_mrat,  sat2_mrat, sat3_mrat, GroupFirstSub, "LMC+SAG+GSE", write_data=True)
#SAG_analogues = get_analogues(sat1lbt, sat2lbt, sat3lbt, sat1_mrat, sat2_mrat, sat3_mrat, GroupFirstSub, "SAG", write_data=True)
#SAG_GSE_analogues = get_analogues(sat1lbt, sat2lbt, sat3lbt, sat1_mrat, sat2_mrat, sat3_mrat, GroupFirstSub, "SAG+GSE", write_data=True)
#GSE_analogues = get_analogues(sat1lbt, sat2lbt, sat3lbt, sat1_mrat, sat2_mrat,  sat3_mrat, GroupFirstSub, "GSE", write_data=True)


######### FIGURES #############

# Histograms of when all the three mergers happen 

"""
fig, ax = plt.subplots(1, 3, figsize=(12, 4), sharey=True)
ax[0].hist(sat1lbt, align='left', rwidth=0.8, color='k', alpha=0.8)
ax[1].hist(sat2lbt, align='left', rwidth=0.8, color='k', alpha=0.8)
ax[2].hist(sat3lbt, align='left', rwidth=0.8, color='k', alpha=0.8)
ax[1].set_xlabel("Lookback time [Gyr]")
ax[0].set_ylabel("N MW Analogues")
ax[0].set_title("3rd merger")
ax[1].set_title("2nd merger")
ax[2].set_title("1st merger")

plt.savefig("histograms_lookback_times_satellites_tophat.png", bbox_inches='tight')
plt.close()
"""


#LMC_SAG_analogues[0]
#LMC_SAG_analogues[0]
#LMC_SAG_analogues[0]

"""
sat1_lbt = SAG_GSE_analogues[3][0]
sat2_lbt = SAG_GSE_analogues[3][1]
sat3_lbt = SAG_GSE_analogues[3][2]

fig, ax = plt.subplots(1, 3, figsize=(12, 4), sharey=True)
ax[0].hist(sat1_lbt, align='left', rwidth=0.8, color='k', alpha=0.8)
ax[1].hist(sat2_lbt, align='left', rwidth=0.8, color='k', alpha=0.8)
ax[2].hist(sat3_lbt, align='left', rwidth=0.8, color='k', alpha=0.8)
ax[1].set_xlabel("Lookback time [Gyr]")
ax[0].set_title("3rd merger")
ax[1].set_title("2nd merger")
ax[2].set_title("1st merger")
plt.savefig("histograms_lookback_times_sag_gse_analogues.png", bbox_inches='tight')
plt.close()
"""


