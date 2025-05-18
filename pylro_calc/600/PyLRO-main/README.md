# PyLRO
The goal of the amorpohous miller plot function is to calculate and display the relative periodicity of miller planes in an amorphous solid. 

Workflow:
* The function accepts an ASE readable crystal structure file of an amorphous crystal. The atomic number of the atom species to be used needs to be specified. Currently, the program only supports tracking atoms that appear only once in the unit cell.
* A sample perfect lattice of points is created using the shape (a,b,c vectors) of the structure. This sample lattice is normalized to unit lengths.
* An optimization scheme then is used to find a scaling factor for this lattice that minimizes the distance between the structure's atomic locations and the closest lattice point. The closest lattice points are the best estimate of where the atomic location would be if the structure were a perfect crystal.
* With atomic locations <x,y,z> and corresponding integer lattice locations <n_1,n_2,n_2>, the products of <x,y,z> and a (scaled) <n_1,n_2,n_3> with any miller plane will give the structure's projections and the ideal lattice's projections.
* The average difference between these projections is considered the order parameter. A smaller order parameter means smaller average difference in projections between the structure's atoms and a perfectly periodic crystal.
* To display this information, an order parameter is calculated for every miller plane of a set resolution
*  Then a unit sphere (<x,y,z> representing <h,k,l>) is plotted where each <x,y,z> on the sphere is inversely scaled to the order parameter. This sphere is also color coded to represent relative order.
* An optional 2 dimensional cross section of the sphere for the xy, xz, and yz planes may also be plotted.

The function will be useful for displaying the process of disorder created in exploration of free energy landscapes. It should also help display directions that may retain bulk periodic properties in an amorphous solids.

Contained in S812 are 812 poscar files for amorphous structures of AlPO4 within the swarm of low symmetry structures found through the stochastic surface walking method described in ["Sheng-Cai Zhu, 2022"](https://pubs.acs.org/doi/10.1021/jacs.2c01717#:~:text=This%20unique%20case%20of%20memory,with%20the%20same%20topological%20linkage.). This program was designed in collaboration with authors of this publication with the primary purpose of corroborating their findings regarding long range topical order in AlPO4.

```
>>> import PyLRO
>>> amor=PyLRO.pylro('S812/POSCAR.1',13)
>>> amor.lattice_fit()
>>> amor.plane_order([0,0,1])
0.16073258564239593
>>> amor.miller_sphere_plot()
```


![miller_sphere](https://github.com/user-attachments/assets/9389d07a-6d29-43d6-b60e-114eea3a4b32)


