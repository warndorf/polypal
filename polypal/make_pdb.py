import sys
import os

def read_atomlist(itp_filename):
    atomlist = []
    with open(itp_filename, 'r') as file:
        for line in file:
            if line.startswith('[ atoms ]'):
                break
        next(file)  # Skip header line
        for line in file:
            if line.startswith('[ bonds ]'):
                break
            parts = line.split()
            if len(parts) >= 8:
                atom_name = parts[4]
                atom_index = parts[0]
                atomlist.append((atom_index, atom_name))
    return atomlist

def read_coordinates(xyz_filename):
    coordinates = []
    with open(xyz_filename, 'r') as file:
        lines = file.readlines()
        for line in lines[2:]:  # Skip the first two lines
            parts = line.strip().split()
            if len(parts) == 4:
                atom, x, y, z = parts
                coordinates.append((atom, float(x), float(y), float(z)))
    return coordinates

def write_pdb(atomlist, coordinates, output_filename):
    with open(output_filename, 'w') as file:
        for (atom_index, atom_name), (atom, x, y, z) in zip(atomlist, coordinates):
            file.write(f"ATOM  {atom_index:5}  {atom_name:<3} TRI X   1    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           {atom}\n")

def AA_to_ZZ_pdb(pdb_filename, output_pdb_filename, AA_input, ZZ_input):
    AA_name = []
    AA_final_atom = []
    ZZ_final_atom = []
    if len(AA_input) == 1:
        with open(pdb_filename, 'r') as file:
            AA_pdb_lines = file.readlines()
        with open(output_pdb_filename, 'w') as file:
            for line in AA_pdb_lines:
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    atom_name = line[12:16].strip()
                    if atom_name == "AA":
                        line = line[:12] + f" {AA_input}0 " + line[16:]
                        AA_final_atom.append(AA_input+'0')
                    AA_name.append(AA_input+'0')
                file.write(line)
    elif len(AA_input) ==2:
        with open(pdb_filename, 'r') as file:
            AA_pdb_lines = file.readlines()
        with open(output_pdb_filename, 'w') as file:
            for line in AA_pdb_lines:
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    atom_name = line[12:16].strip()
                    if atom_name == "AA":
                        line = line[:12] + f" {AA_input}0" + line[16:]
                        AA_final_atom.append(AA_input+'0')
                    AA_name.append(AA_input+'0')
                file.write(line)
    if len(ZZ_input) == 1:
        with open(output_pdb_filename, 'r') as file:
            ZZ_pdb_lines = file.readlines()
        with open(output_pdb_filename, 'w') as file:
            AA_cant_match = AA_name[0]
            for line in ZZ_pdb_lines:
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    atom_name = line[12:16].strip()
                    if atom_name == "ZZ":
                        if ZZ_input.isalpha():
                            ZZ_change = ZZ_input+'0'
                            if ZZ_change == AA_cant_match:
                                line = line[:12] + f" {ZZ_input}00" + line[16:]
                                ZZ_final_atom.append(ZZ_input+'00')
                            elif ZZ_change != AA_cant_match:
                                line = line[:12] + f" {ZZ_input}0 " + line[16:]
                                ZZ_final_atom.append(ZZ_input+'0')
                file.write(line)
    elif len(ZZ_input) == 2:
        with open(output_pdb_filename, 'r') as file:
            ZZ_pdb_lines = file.readlines()
        with open(output_pdb_filename, 'w') as file:
            AA_cant_match = AA_name[0]
            for line in ZZ_pdb_lines:
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    atom_name = line[12:16].strip()
                    if atom_name == "ZZ":
                        if ZZ_input.isalpha():
                            ZZ_change = ZZ_input+'0'
                            if ZZ_change == AA_cant_match:
                                line = line[:12] + f" {ZZ_input}x" + line[16:]
                                ZZ_final_atom.append(ZZ_input+'x') 
                            elif ZZ_change != AA_cant_match:
                                line = line[:12] + f" {ZZ_input}0" + line[16:]
                                ZZ_final_atom.append(ZZ_input+'0')
                file.write(line)
    print()
    print(f'AA changed to "{AA_final_atom[0]}" and ZZ changed to "{ZZ_final_atom[0]}" in {pdb_filename}, and saved to {output_pdb_filename}') 
    print("Please ensure GROMACS can read these atoms correctly, and add any atom names to vdwraddi.dat if it cannot")
    print()


def make_pdb(command=None):
    mode_flag = None
    pdb_filename = None
    qforce_itp = None
    xyz_filename = None
    output_pdb_filename = None

    if command:

        if len(command) < 2:
            print("ERROR:")
            print("Could not perform conversion. Please see polypal pdb -h for help.")
            print() 
            exit(1)

        if command[0].startswith('-') and command[0] != '-p' and command[0] != '-q' and command[0] != '-x':
            mode_flag = command[0]
            start_index = 1 # 
        else:
            start_index = 0 

        for i in range(start_index, len(command), 2):
            flag = command[i]
            file_or_dir = command[i + 1]
            if flag == '-p':
                if not os.path.isfile(file_or_dir):
                    print("ERROR:")
                    print(f"could not find: {file_or_dir}")
                    print()
                    exit(1)
                pdb_filename = file_or_dir
                output_pdb_filename = "".join(pdb_filename.split('.')[:-1])+'_final.pdb'
            elif flag == '-x':
                if not os.path.isfile(file_or_dir):
                    print("ERROR:")
                    print(f"could not find: {file_or_dir}")
                    print()
                    exit(1)
                xyz_filename = file_or_dir
            elif flag == '-q':
                if not os.path.isfile(file_or_dir):
                    print("ERROR:")
                    print(f"could not find: {file_or_dir}")
                    print()
                    exit(1)
                qforce_itp = file_or_dir 

    if mode_flag == '-xyz' and qforce_itp and xyz_filename:
        pdb_from_xyz_filename = os.path.splitext(os.path.basename(xyz_filename))[0] +'.pdb'
        atomlist = read_atomlist(qforce_itp)
        coordinates = read_coordinates(xyz_filename)
        write_pdb(atomlist, coordinates, pdb_from_xyz_filename)
        print(f"Output written to {pdb_from_xyz_filename}")
        print()

    if mode_flag == '-gromacs' and pdb_filename:
        while True:
            AA_input = input("Atom 'AA' found. Please enter single letter corresponding to atom element (ex: C for carbon): ").strip().upper()
            if AA_input.isalpha():
                break
            else:
                print("ERROR:")
                print("Invalid input. Please enter valid element name.")
                print()
        while True:
            ZZ_input = input("Atom 'ZZ' found. Please enter single letter corresponding to atom element: ").strip().upper()
            if ZZ_input.isalpha():
                break
            else:
                print("ERROR:")
                print("Invalid input. Please enter a valid element name.")
                print()
        AA_to_ZZ_pdb(pdb_filename, output_pdb_filename, AA_input, ZZ_input)
        print(f"Output written to {output_pdb_filename}")
        print()

    
    #else:
    #    print("ERROR:")
    #    print("Could not perform conversion. Please see polypal pdb -h for help.")
    #    print()   
        

#if __name__ == "__main__":
#    main()
