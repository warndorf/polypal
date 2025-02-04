#!/bin/bash
#SBATCH -J test
#SBATCH -N 1 
#SBATCH -n 4
#SBATCH --mem=150G
#SBATCH --partition=xeon-p8

export PATH=$PATH:/home/gridsan/mwarndorf/software/bin
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/gridsan/mwarndorf/software/bin
source ~/gromacs/2023.2-cpu/bin/GMXRC


which gmx_mpi

gmx_mpi grompp -f minim.mdp -c ini.gro -p 8mer_final.top -o em.tpr -maxwarn 2

mpirun gmx_mpi mdrun -v -deffnm em

gmx_mpi grompp -f equil_nvt_300K_100ps.mdp -c em.gro -p 8mer_final.top -o nvt.tpr -maxwarn 2

mpirun gmx_mpi mdrun -v -deffnm nvt

gmx_mpi grompp -f equil_npt_200ps_100bar_300K.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p 8mer_final.top -o npt100bar.tpr -maxwarn 2

mpirun gmx_mpi mdrun -v -deffnm npt100bar

gmx_mpi grompp -f equil_npt_200ps_1bar_300K.mdp -c npt100bar.gro -p 8mer_final.top -o npt1bar.tpr -maxwarn 2

mpirun gmx_mpi mdrun -v -deffnm npt1bar

gmx_mpi grompp -f equil_npt_200ps_1000bar_300K.mdp -c npt1bar.gro -r npt1bar.gro -t npt1bar.cpt -p 8mer_final.top -o npt1000bar.tpr -maxwarn 2

mpirun gmx_mpi mdrun -v -deffnm npt1000bar

gmx_mpi grompp -f anneal_300_to_300.mdp -c npt1000bar.gro -r npt1000bar.gro -t npt1000bar.cpt -p 8mer_final.top -o anneal1.tpr -maxwarn 2

mpirun gmx_mpi mdrun -v -deffnm anneal1

gmx_mpi grompp -f anneal_300_to_300.mdp -c anneal1.gro -r anneal1.gro -t anneal1.cpt -p 8mer_final.top -o anneal2.tpr -maxwarn 2

mpirun gmx_mpi mdrun -v -deffnm anneal2

gmx_mpi grompp -f anneal_300_to_300.mdp -c anneal2.gro -r anneal2.gro -t anneal2.cpt -p 8mer_final.top -o anneal3.tpr -maxwarn 2

mpirun gmx_mpi mdrun -v -deffnm anneal3

gmx_mpi grompp -f anneal_300_to_300.mdp -c anneal3.gro -r anneal3.gro -t anneal3.cpt -p 8mer_final.top -o anneal4.tpr -maxwarn 2

mpirun gmx_mpi mdrun -v -deffnm anneal4

gmx_mpi grompp -f anneal_300_to_300.mdp -c anneal4.gro -r anneal4.gro -t anneal4.cpt -p 8mer_final.top -o anneal5.tpr -maxwarn 2

mpirun gmx_mpi mdrun -v -deffnm anneal5

gmx_mpi grompp -f anneal_300_to_300.mdp -c anneal5.gro -r anneal5.gro -t anneal5.cpt -p 8mer_final.top -o anneal6.tpr -maxwarn 2

mpirun gmx_mpi mdrun -v -deffnm anneal6

gmx_mpi grompp -f anneal_300_to_300.mdp -c anneal6.gro -r anneal6.gro -t anneal6.cpt -p 8mer_final.top -o anneal7.tpr -maxwarn 2

mpirun gmx_mpi mdrun -v -deffnm anneal7

gmx_mpi grompp -f anneal_300_to_300.mdp -c anneal7.gro -r anneal7.gro -t anneal7.cpt -p 8mer_final.top -o anneal8.tpr -maxwarn 2

mpirun gmx_mpi mdrun -v -deffnm anneal8

gmx_mpi grompp -f anneal_300_to_300.mdp -c anneal8.gro -r anneal8.gro -t anneal8.cpt -p 8mer_final.top -o anneal9.tpr -maxwarn 2

mpirun gmx_mpi mdrun -v -deffnm anneal9

gmx_mpi grompp -f anneal_300_to_300.mdp -c anneal9.gro -r anneal9.gro -t anneal9.cpt -p 8mer_final.top -o anneal10.tpr -maxwarn 2

mpirun gmx_mpi mdrun -v -deffnm anneal10

gmx_mpi grompp -f npt_5ns_1bar_300K.mdp -c anneal10.gro -r anneal10.gro -t anneal10.cpt -p 8mer_final.top -o production5ns.tpr -maxwarn 2

mpirun gmx_mpi mdrun -v -deffnm production5ns
