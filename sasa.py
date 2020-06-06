import matplotlib.pyplot as plt

from pyrosetta import *
init("-holes:dalphaball /net/shared/rosetta_builds/bcov_stable/20-04-27-broken/main/source/external/DAlpahBall/DAlphaBall.gcc")

import time
import pandas as pd
df = pd.DataFrame(columns=['name', 'type', 'sasa'])
daltimes = pd.DataFrame(columns=['name', 'Residue Count', 'Time (s)'])
voxtimes = pd.DataFrame(columns=['name', 'Residue Count', 'Time (s)'])

def get_per_atom_sasa_dalpha_ball(pose, probe_size=1.1):
    atoms = rosetta.core.id.AtomID_Map_bool_t()
    atoms.resize(pose.size())
    for i in range(1, pose.size()+1):
        atoms.resize( i, pose.residue(i).natoms(), True)
    surf_vol = rosetta.core.scoring.packing.get_surf_vol( pose, atoms, probe_size)
    return surf_vol

def get_per_atom_sasa_normal(pose, probe_size=1.1):
    sasas = rosetta.core.id.AtomID_Map_double_t()
    rsd_sasa = rosetta.utility.vector1_double()
    rosetta.core.scoring.calc_per_atom_sasa(pose, sasas, rsd_sasa, probe_size, False)
    return sasas

pdbs = []

with open('short_pdbs.csv', 'r') as file:
    pdbs = [pdb.rstrip() for pdb in file.readlines()]


start = time.time()
for pdb in pdbs:
    
    pose = pose_from_file(pdb)
    
    dalstart = time.time()
    surf_vol = get_per_atom_sasa_dalpha_ball(pose)

    dal = 0
    
    for i in range(1, pose.size()+1):
        res = pose.residue(i)
        for j in range(1, res.natoms()+1):
            dal += surf_vol.surf(i, j)

    dalend = time.time()
    
    df = df.append(pd.Series([pdb, 'Power Method', dal], index=df.columns), ignore_index=True)
    daltimes = daltimes.append(pd.Series([pdb, pose.size(), dalend - dalstart], index=daltimes.columns), ignore_index=True)

mid = time.time()
for pdb in pdbs:
    
    pose = pose_from_file(pdb)
    
    voxstart = time.time()
    normal_sasa = get_per_atom_sasa_normal(pose)

    sasa = 0
    for i in range(1, pose.size()+1):
        res = pose.residue(i)
        for j in range(1, res.natoms()+1):
            sasa += normal_sasa(i, j)
    
    voxend = time.time()
    
    df = df.append(pd.Series([pdb, 'Shrake-Rupley', sasa], index=df.columns), ignore_index=True)
    voxtimes = voxtimes.append(pd.Series([pdb, pose.size(), voxend - voxstart], index=voxtimes.columns), ignore_index=True)

end = time.time()


times = open('times.txt', 'w')
times.write("Power Method Time: {}\n".format(mid - start))
times.write("Shrake-Rupley Time: {}\n".format(end - mid))
times.close()

df.to_csv('values.csv', index=False, header=False)
daltimes.to_csv('daltimes.csv', index=False, header=False)
voxtimes.to_csv('voxtimes.csv', index=False, header=False)


print("Power Method Time: {}".format(mid - start))
print("Shrake-Rupley Time: {}".format(end - mid))

fig, axes = plt.subplots(nrows=1, ncols=2, sharex=True, sharey=True)

# drop sharex, sharey, layout & add ax=axes
df.hist(column='sasa',by='type', ax=axes)

# set title and axis labels
plt.suptitle('Computing SASA: Power Method vs. Shrake-Rupley', x=0.5, y=1.00, ha='center', fontsize='large')
fig.text(0.5, 0.02, 'Calculated SASA', ha='center')
fig.text(0.01, 0.5, 'Number of Structures', va='center', rotation='vertical')

plt.savefig('sasa.png')


fig2 =  daltimes.plot.scatter(x='Residue Count', y='Time (s)')
plt.title("Power Method Time vs. Residue Count")
plt.savefig('daltimes.png')

fig3 = voxtimes.plot.scatter(x='Residue Count', y='Time (s)')
plt.title("Shrake-Rupley Time vs. Residue Count")
plt.savefig('voxtimes.png')
