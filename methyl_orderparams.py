import mdtraj as mt
import numpy as np
import pandas as pd
import sys
import argparse

parser = argparse.ArgumentParser()

#-db DATABASE -u USERNAME -p PASSWORD -size 20000
parser.add_argument("-t", "--topology", dest = "topology", default = "../dinamica/reducida/redSqrR.prmtop", help="Topology for MD trajectory")
parser.add_argument("-x", "--trajectory", dest = "trajectory", default = ['../dinamica/reducida/redSqrR.MD2.nc'], help="Python list of trayectory names")
parser.add_argument("-o", "--output",dest ="outname",default = 'orderparams', help="Output name")
parser.add_argument("-r", "--residues",dest ="residues",default = 'ILVAM', help="Methyl aminoacids")

args = parser.parse_args()

#import matplotlib.pyplot as plt

#Cuáles son los metilos presentes en una proteina?
#ILVAM: Isoleucina, leucina, valina, alanina y metionina
#
#Vamos a tener esos vectores que unen Cprevio a Cmetilo y normalizarlo. Luego usamos la formula de https://www.cell.com/action/showPdf?pii=S0006-3495%2811%2900783-1
#
#Hay que ser cuidadoso con el alineamiento. Tratar de que sea lo mejor posible
#S^2 = 3/2 [(avg(X^2)^2+(avg(X^2)^2+(avg(Y^2)^2+(avg(Z^2)^2(avg(X*Y)^2+(avg(X*Z)^2+(avg(Y*Z)^2] - 1/2
#Recomendable hacerlo en trayectorias de aprox 1us

#Armamos dict con atomos correspondientes a cada AA
if args.residues=="ILVAM":
    methyl_vectors={"ILE":[["CB","CG2"],["CG1","CD1"]],"LEU":[["CG","CD1"],["CG","CD2"]],"VAL":[["CB","CG1"],["CB","CG2"]],"ALA":[["CA","CB"]],"MET":[["SD","CE"]]}
else:
    print("Only ILVAM scheme supported")
    exit()

trajectories=args.trajectory
topology=args.topology

protein_atoms=range(0,3200)
#Importamos trayectorias
trajs= mt.load(trajectories, atom_indices=protein_atoms,top=topology)
#Alineamo
trajs.superpose(reference=trajs)
df=pd.DataFrame()
#df.columns=["resid","resname","n_methyl","x","y","z"]
#identificamos ILVAM y levantamos nombres de atomos correspondientes
table, bonds = trajs.topology.to_dataframe()
correct_atoms=table.merge(pd.DataFrame([{'resName': k, 'name': o} 
                       for k, v in methyl_vectors.items() for i in v for o in i])).sort_values("resSeq")
gen_list=[]
#ahora que tenemos los átomos de ILVAM involucrados en metilos, buscamos su XYZ
for index, atom in correct_atoms.iterrows():
   atomselect="resid "+str(atom.resSeq)+" and name "+atom[1]+""
   atomindex=trajs.topology.select(atomselect)[0]
   xyz=trajs.xyz[:,atomindex,:]
   gen_list.append([atom.resSeq+1,atom.resName,atom[1],xyz[:,0],xyz[:,1],xyz[:,2]])
gen_list=pd.DataFrame(gen_list)
gen_list.columns=["resSeq","resName","atName","x","y","z"]

#Ahora hacemos un matcheo entre atomos que forman enlaces C-C (o S-C en metionina), calculamos el vector normalizado correspondiente a dicho enlace
# y directamente calculamos el parametro de orden
vector_list=[]
order_params_methyl=[]
for index, atom in gen_list.iterrows():
    for option,pair in methyl_vectors.items():
        if atom.resName==option:
            for bond in pair:
               if atom.atName==bond[0]:
                   appendable=[a for a in atom]+list(gen_list.loc[ (gen_list.resSeq==atom.resSeq) & (gen_list.atName==bond[1])].values[0])
                   X=(appendable[3]-appendable[9])/np.sqrt((appendable[3]-appendable[9])**2+(appendable[4]-appendable[10])**2+(appendable[5]-appendable[11])**2)
                   Y=(appendable[4]-appendable[10])/np.sqrt((appendable[3]-appendable[9])**2+(appendable[4]-appendable[10])**2+(appendable[5]-appendable[11])**2)
                   Z=(appendable[5]-appendable[11])/np.sqrt((appendable[3]-appendable[9])**2+(appendable[4]-appendable[10])**2+(appendable[5]-appendable[11])**2)

                   #pasamos XYZ a parametro de orden
                   #S^2 = 3/2 [(avg(X^2)^2+(avg(X^2)^2+(avg(Y^2)^2+(avg(Z^2)^2(avg(X*Y)^2+(avg(X*Z)^2+(avg(Y*Z)^2] - 1/2


                   S2=1.5*(np.mean(X**2)**2+np.mean(Y**2)**2+np.mean(Z**2)**2+np.mean(X*Z)**2+np.mean(X*Y)**2+np.mean(Y*Z)**2)-0.5

                   normalized_vector=[appendable[0],appendable[1],appendable[2],appendable[8],S2] 
                   order_params_methyl.append(normalized_vector)

#Grabamos todo en un csv
opm_df=pd.DataFrame(order_params_methyl)
opm_df.columns=["resid","resname","atom1","atom2","S**2"]
opm_df.to_csv(args.outname+'.csv',sep='\t')


#armamos grafico


