import numpy as np
import matplotlib.pyplot as plt
import illustris_python as il
from illustris_python.sublink import maxPastMass
from illustris_python.util import partTypeNum 
from illustris_python.groupcat import loadHeader

# 1. Get all halos with three mergers
# 2. Find mergers max mass and snap at which they have max mass, ids in the halo
# catalogue. 
# 3. Find times and snaps where the halos enter the halo and when they merge,
# ids in halo catalgue.


def get_redshifts(basePath, max_num_snapshots=136):
    # modified from: https://www.tng-project.org/data/forum/topic/369/snapshots-and-redshifts/
    redshifts = np.zeros((max_num_snapshots, 2), dtype='float32' )

    for i in range(0, max_num_snapshots):
        h = loadHeader(basePath,i)
        redshifts[0, i] = h['Redshift']
        redshifts[1, i] = i

    np.savetxt("redshifts_snaps_illustris-3.txt", redshifts)
    return 0

def maxPastMass(tree, index, partType='stars'):
      """ Get maximum past mass (of the given partType) along the main branch of
      a subhalo
      specified by index within this tree. """
      ptNum = partTypeNum(partType)

      branchSize = tree['MainLeafProgenitorID'][index] -  tree['SubhaloID'][index] + 1
      masses = tree['SubhaloMassType'][index: index + branchSize, ptNum]
      max_mass = np.argmax(masses)
      return np.max(masses), index+max_mass

def numMergersIDs(tree, minMassRatio=1e-10, massPartType='stars', index=0):
    """ 
    Calculate the number of mergers in this sub-tree (optionally above some mass ratio threshold). 
    
    Parameters:
    -----------

    index : 

    """
    # verify the input sub-tree has the required fields
    reqFields = ['SubhaloID', 'NextProgenitorID', 'MainLeafProgenitorID',
                 'FirstProgenitorID', 'SubhaloMassType']

    if not set(reqFields).issubset(tree.keys()):
        raise Exception('Error: Input tree needs to have loaded fields: '+', '.join(reqFields))

    numMergers   = 0
    invMassRatio = 1.0 / minMassRatio
    id_subhalos  = []
    snapnum_subhalos = []
    snapnum_fp_max = []
    ratios = []
    id_subhalos_fpmax = []
    # walk back main progenitor branch
    rootID = tree['SubhaloID'][index]
    fpID   = tree['FirstProgenitorID'][index]

    while fpID != -1:
        fpIndex = index + (fpID - rootID)
        fpMass  = maxPastMass(tree, fpIndex, massPartType)

        # explore breadth
        npID = tree['NextProgenitorID'][fpIndex]
        snapnum = tree['SnapNum'][fpIndex]
        pos = tree['SubhaloPos'][fpIndex]

        while npID != -1:
            npIndex = index + (npID - rootID)
            npMass, max_mass_index  = maxPastMass(tree, npIndex, massPartType)
            subID = tree['SubfindID'][max_mass_index]
            snapnum = tree['SnapNum'][max_mass_index]
            # mass of FP branch at max sat. mass
            branchSize = tree['MainLeafProgenitorID'][fpIndex] -  tree['SubhaloID'][fpIndex] + 1
            snapnum_FP = tree['SnapNum'][fpIndex: fpIndex + branchSize]
   
            # index in tree of the FP at the time of max sat mass.
            # First progenitor index
            if len(snapnum_FP[snapnum_FP>=snapnum])>1:         
                FP_max_index = np.argmin(snapnum_FP[snapnum_FP>=snapnum])+fpIndex 
                fpMass_satmax = tree['SubhaloMassType'][FP_max_index, 4] # 4 is parttype
                fpmax_snapnum = tree['SnapNum'][FP_max_index]
                fpmax_subID = tree['SubfindID'][FP_max_index]
        
                # count if both masses are non-zero, and ratio exceeds threshold
                if fpMass_satmax > 0.0 and npMass > 0.0 :
                    ratio = npMass / fpMass_satmax

                    if ratio >= minMassRatio and ratio <= invMassRatio:
                        numMergers += 1
                        id_subhalos.append(subID)
                        snapnum_subhalos.append(snapnum)
                        snapnum_fp_max.append(fpmax_snapnum)
                        id_subhalos_fpmax.append(fpmax_subID)
                        ratios.append(ratio)
                npID = tree['NextProgenitorID'][npIndex]
            else : 
                npID = -1       

        fpID = tree['FirstProgenitorID'][fpIndex]

    return numMergers, id_subhalos, snapnum_subhalos, id_subhalos_fpmax,  snapnum_fp_max, ratios


#basePath = '/mnt/home/nico/ceph/illutris/Illustris-3/output'
#basePath = '/mnt/home/nico/ceph/illutris/Illustris-1/output'
basePath = '/mnt/home/nico/ceph/illutris/TNG100-1/output'


# z=0 snap nmumber:
# 135 -> Illustris 
# 99 -> TNG

snap = 99

Grouphalos = il.groupcat.loadHalos(basePath, snap, fields=['GroupFirstSub'])
#GroupFirstSub_Mass = il.groupcat.loadHalos(basePath,snap, fields=['Group_M_Crit200'])

GroupFirstSub_Mass = il.groupcat.loadHalos(basePath, snap, fields=['Group_M_TopHat200'])
Group_mass = GroupFirstSub_Mass * 1E10/0.704


Groupsubhalos = il.groupcat.loadSubhalos(basePath, snap, fields=['SubhaloMass'])

# Select MW-like halos
mw_cuts = np.where((Group_mass < 2E12) & (Group_mass > 0.5E12))
print(mw_cuts[0])
print("N MW-like halos", len(mw_cuts[0]))

Nmw_analogues = len(mw_cuts[0])



ratio = 1/20.
fields = ['SubhaloMass','SubfindID','SnapNum','SubhaloID','NextProgenitorID','MainLeafProgenitorID','FirstProgenitorID','SubhaloMassType','SubhaloPos']

all_halos = np.zeros((Nmw_analogues, 19))
k=0

Nmergers = np.zeros(Nmw_analogues)


for i in range(100, Nmw_analogues):
    print(i, Grouphalos[i])
    tree = il.sublink.loadTree(basePath, snap, Grouphalos[i], fields=fields)
    #numMergers = il.sublink.numMergers(tree, minMassRatio=ratio)
    x = numMergersIDs(tree, minMassRatio=ratio)
    Nmergers[i] = x[0]
    if x[0]==3:
       all_halos[k] = np.array([Grouphalos[mw_cuts][i], x[0], x[1][0], x[1][1], x[1][2], x[2][0],
            x[2][1], x[2][2], x[3][0], x[3][1], x[3][2], x[4][0], x[4][1],
            x[4][2], x[5][0], x[5][1], x[5][2], tree["SubfindID"][0], tree["SnapNum"][0]])

       k+=1
np.savetxt("halos_three_mergers_1to20_Il-1.txt", all_halos[:k], fmt=["%d", "%d", "%d",
            "%d", "%d", "%d", "%d", "%d", "%d", "%d", "%d", "%d", "%d", "%d", "%.3f",
           "%.3f", "%.3f", "%d", "%d"], header="Group halo MW analogue, N mergers, id sub1, id sub2, idsub3, snap sub1, snap sub2, snap sub3, id fp1 max, id fp2 max, id fp3 max, merger rat 1, merger rat2, merger rat3, subfindID, SnapNum0")        


np.savetxt("TNG-100_nmergers_mw_analogues_tophat.txt", Nmergers, fmt="%d", header="N mergers")
