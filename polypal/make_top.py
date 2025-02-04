import sys
import re
import os

def read_itp_file(itp_filename):
    itp_atoms = []
    with open(itp_filename, 'r') as file:
        for line in file:
            if line.startswith('[ atoms ]'):
                break
        next(file)  # Skip header line
        for line in file:
            if line.startswith('['):
                break
            parts = line.split()
            if len(parts) > 4:
                index = int(parts[0])
                atom_name = parts[4]
                atom_type = parts[1]
                itp_atoms.append((index, atom_name, atom_type))
    return itp_atoms

def read_pdb_file(pdb_filename):
    pdb_atoms = []
    with open(pdb_filename, 'r') as file:
        for line in file:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                atom_name = line[12:16].strip()
                pdb_atoms.append(atom_name)
    return pdb_atoms

def find_matching_atoms(itp_atoms, pdb_atoms, output_file):
    full_list = []
    pdb_atoms_set = set(pdb_atoms)  # Convert list to set for faster lookup
    for index, atom_name, atom_type in itp_atoms:
        if atom_name in pdb_atoms_set:
            full_list.append((index, atom_name, atom_type))
    with open(output_file, 'w') as file:
        file.write(f"[ mapping ]\n" )
    with open(output_file, 'a') as file:
        for (index, atom_name, atom_type) in full_list:
            file.write(f"   {atom_name}   {atom_type} \n" )   
    return full_list

def process_bonds(assignment_file, start_bond, end_bond, full_list):
    section_bonds = []
    with open(assignment_file, 'r') as file:
        lines = file.readlines()
    recording = False
    for line in lines:
        if start_bond in line:
            recording = True
            continue
        elif end_bond in line:
            break
        if recording:
            section_bonds.append(line.strip())
    matching_bonds = []
    set_atom = [] 
    set_bond_column = {}
    for line in section_bonds:
        columns = line.split()
        if len(columns) >= 8: 
            bond_column = columns[0]  
            atom1 = columns[4]  
            atom2 = columns[5]  
            set_atom.append((atom1, atom2))
            set_bond_column[(atom1, atom2)] = bond_column
    set_name = {}
    for line in full_list:
        if isinstance(line, tuple):
            line = " ".join(map(str, line))  
        parts = line.split()
        if len(parts) == 3:
            atom_index = parts[0]  
            atom_name = parts[1]  
            set_name[atom_index] = atom_name  
    matching_bonds = []
    for atom1, atom2 in set_atom:
        if atom1 in set_name and atom2 in set_name:
            bond_column = set_bond_column.get((atom1, atom2), "")
            matching_bonds.append((set_name[atom1], set_name[atom2], bond_column))
    return matching_bonds

def process_angles(assignment_file, start_angle, end_angle, full_list):
    section_angles = []
    with open(assignment_file, 'r') as file:
        lines = file.readlines()
    recording = False
    for line in lines:
        if start_angle in line:
            recording = True
            continue
        elif end_angle in line:
            break
        if recording:
            section_angles.append(line.strip())
    matching_angles = []
    set_atom = {}
    for line in section_angles:
            columns = line.split()
            if len(columns) == 12:
                    angle_column = columns[0]
                    atom1 = int(columns[6])
                    atom2 = int(columns[7])
                    atom3 = int(columns[8])
                    set_atom[(atom1, atom2, atom3)] = angle_column
            if len(columns) == 10:
                    angle_column = columns[0]
                    atom1 = int(columns[4])
                    atom2 = int(columns[5])
                    atom3 = int(columns[6])
                    set_atom[(atom1, atom2, atom3)] = angle_column
    set_name = {}
    for line in full_list:
        if isinstance(line, tuple):
            line = " ".join(map(str, line))
        parts = line.split()
        if len(parts) == 3:
            atom_index = int(parts[0])
            atom_name = parts[1]
            set_name[atom_index] = atom_name
    for (atom1, atom2, atom3), angle_column in set_atom.items():
        if atom1 in set_name and atom2 in set_name and atom3 in set_name:
            matching_angles.append((set_name[atom1], set_name[atom2], set_name[atom3], angle_column))
    return matching_angles

def process_dihedrals(assignment_file, start_dihedral, end_dihedral, full_list):
    section_dihedrals = []
    with open(assignment_file, 'r') as file:
        lines = file.readlines()
    recording = False
    for line in lines:
        if start_dihedral in line:
            recording = True
            continue
        elif end_dihedral in line:
            break
        if recording:
            section_dihedrals.append(line.strip())
    matching_dihedrals = []
    set_atom = {}
    for line in section_dihedrals:
            columns = line.split()
            if len(columns) >= 8:
                    dihedral_column = columns[0]
                    atom1 = int(columns[4])
                    atom2 = int(columns[5])
                    atom3 = int(columns[6])
                    atom4 = int(columns[7])
                    set_atom[(atom1, atom2, atom3, atom4)] = dihedral_column
    set_name = {}
    for line in full_list:
        if isinstance(line, tuple):
            line = " ".join(map(str, line))
        parts = line.split()
        if len(parts) == 3:
            atom_index = int(parts[0])
            atom_name = parts[1]
            set_name[atom_index] = atom_name
    for (atom1, atom2, atom3, atom4), dihedral_column in set_atom.items():
        if atom1 in set_name and atom2 in set_name and atom3 in set_name and atom4 in set_name:
            matching_dihedrals.append((set_name[atom1], set_name[atom2], set_name[atom3], set_name[atom4], dihedral_column))
    return matching_dihedrals

def process_flex_diheds(assignment_file, start_flex_diheds, end_flex_diheds, full_list):
    section_flex_diheds = []
    with open(assignment_file, 'r') as file:
        lines = file.readlines()
    recording = False
    for line in lines:
        if start_flex_diheds in line:
            recording = True
            continue
        elif end_flex_diheds in line:
            break
        if recording:
            section_flex_diheds.append(line.strip())
    matching_flex_diheds = []
    set_atom = {}
    for line in section_flex_diheds:
            columns = line.split()
            if len(columns) >= 12:
                    dihedral_column = columns[0]
                    atom1 = int(columns[8])
                    atom2 = int(columns[9])
                    atom3 = int(columns[10])
                    atom4 = int(columns[11])
                    set_atom[(atom1, atom2, atom3, atom4)] = dihedral_column
    set_name = {}
    for line in full_list:
        if isinstance(line, tuple):
            line = " ".join(map(str, line))
        parts = line.split()
        if len(parts) == 3:
            atom_index = int(parts[0])
            atom_name = parts[1]
            set_name[atom_index] = atom_name
    for (atom1, atom2, atom3, atom4), dihedral_column in set_atom.items():
        if atom1 in set_name and atom2 in set_name and atom3 in set_name and atom4 in set_name:
            matching_flex_diheds.append((set_name[atom1], set_name[atom2], set_name[atom3], set_name[atom4], dihedral_column))
    return matching_flex_diheds

def write_formatted_bonds(matching_bonds, output_file):
    with open(output_file, 'a') as file:
        # Determine the maximum width for each column
        col_widths = [max(len(str(item)) for item in col) for col in zip(*matching_bonds)]
        # Write header (optional, for clarity)
        header = ['[ bonds ]']
        file.write(" ".join(f"{h:<{w}}" for h, w in zip(header, col_widths)) + '\n')
        # Write each line with formatting
        for bonds in matching_bonds:
            atom1, atom2, bond_type = bonds
            file.write("   {:<5} {:<10} {:<5}\n".format(atom1, atom2, bond_type))

def write_formatted_angles(matching_angles, output_file):
    with open(output_file, 'a') as file:
        # Determine the maximum width for each column
        col_widths = [max(len(str(item)) for item in col) for col in zip(*matching_angles)]
        # Write header (optional, for clarity)
        header = ['[ angles ]']
        file.write(" ".join(f"{h:<{w}}" for h, w in zip(header, col_widths)) + '\n')
        # Write each line with formatting
        for angles in matching_angles:
            atom1, atom2, atom3, angle_name = angles
            file.write("   {:<5} {:<5} {:<5} {:<5}\n".format(atom1, atom2, atom3, angle_name))

def write_formatted_dihedrals(matching_dihedrals, matching_flex_diheds, output_file):
    with open(output_file, 'a') as file:
        # Determine the maximum width for each column
        col_widths = [max(len(str(item)) for item in col) for col in zip(*matching_dihedrals)]
        # Write header (optional, for clarity)
        header = ['[ dihedrals ]']
        file.write(" ".join(f"{h:<{w}}" for h, w in zip(header, col_widths)) + '\n')
        # Write each line with formatting
        for angles in matching_dihedrals:
            atom1, atom2, atom3, atom4, dihedral_name = angles
            file.write("   {:<5} {:<5} {:<5} {:<5} {:<5}\n".format(atom1, atom2, atom3, atom4, dihedral_name))
    with open(output_file, 'a') as file:
        for angles in matching_flex_diheds:
            atom1, atom2, atom3, atom4, dihedral_name = angles
            file.write("   {:<5} {:<5} {:<5} {:<5} {:<5}\n".format(atom1, atom2, atom3, atom4, dihedral_name))  

def write_terminal_types(output_file, final_output_file, AA_name, ZZ_name):
    with open(output_file, 'r') as file:
        recording = False
        for line in file:
            parts = line.strip().split()
            if len(parts) == 2:
                recording = True
                name_from_list = parts[0]
                atom_type = parts[1]
                if AA_name.strip() == name_from_list.strip():
                    with open(final_output_file, 'a') as file:
                        file.write("[ nterminal ]" + "\n")
                        file.write("   {:<5} {:<5}\n".format('AA', atom_type))
                if ZZ_name.strip() == name_from_list.strip():
                    with open(final_output_file, 'a') as file:
                        file.write("[ cterminal ]" + "\n")
                        file.write("   {:<5} {:<5}\n".format('ZZ', atom_type))
            elif recording:
                break


def replace_AA_and_ZZ(output_file, final_output_file, pdb_file, renamed_pdb):
    with open(output_file, 'r') as file:
        content_AA = file.read()
        print(f'Processing {pdb_file}')
        print()
        while True:
            AA_name = input("Enter AA name:")
            pattern_AA = r'\b' + re.escape(AA_name) + r'\b'
            if re.search(pattern_AA, content_AA):
                updated_AA_content = re.sub(pattern_AA, "AA", content_AA)
                with open(final_output_file, 'w') as file:
                    file.write(updated_AA_content)
                break
            else:
                print()                
                print("ERROR:")
                print(f'Could not find atom: {AA_name}')
                print()
    with open(final_output_file, 'r') as file:
        content_ZZ = file.read()
        while True:
            ZZ_name = input("Enter ZZ name:")
            pattern_ZZ = r'\b' + re.escape(ZZ_name) + r'\b'
            if re.search(pattern_ZZ, content_ZZ):
                updated_ZZ_content = re.sub(pattern_ZZ, "ZZ", content_ZZ)
                with open(final_output_file, 'w') as file:
                    file.write(updated_ZZ_content + " ")
                break
            else:
                print()                
                print("ERROR:")               
                print(f'Could not find atom: {ZZ_name}')
                print()
    with open(pdb_file, 'r') as file:
        content = file.readlines()
        with open(renamed_pdb, 'w') as output_file:
                for line in content:
                    if line.startswith("ATOM"):
                        current_atom_name = line[12:16].strip()
                        if current_atom_name == AA_name:
                            updated_line = line[:12] + " AA " + line[16:]
                            output_file.write(updated_line)
                        else:
                            output_file.write(line)
                    else:
                        output_file.write(line)
    with open(renamed_pdb, 'r') as file:
        content = file.readlines()
        with open(renamed_pdb, 'w') as output_file:
                for line in content:
                    if line.startswith("ATOM"):
                        current_atom_name = line[12:16].strip()
                        if current_atom_name == ZZ_name:
                            updated_line = line[:12] + " ZZ " + line[16:]
                            output_file.write(updated_line)
                        else:
                            output_file.write(line)
                    else:
                        output_file.write(line)   
    return AA_name, ZZ_name

def make_top(command=None):

    pdb_file = None
    itp_file = None
    assignment_file = None

    if command:
        if len(sys.argv) < 1:
            print("ERROR:")
            print("Could not perform conversion. Please see polypal pdb -h for help.")
            print() 
            exit(1)

        if command[0].startswith('-') and command[0] != '-p' and command[0] != '-q' and command[0] != '-a':
            start_index = 1
        else:
            start_index = 0 

        for i in range(start_index, len(command), 2):
            flag = command[i]
            file_or_dir = command[i + 1]
            if flag == '-p':
                if not os.path.isfile(file_or_dir):
                    print("ERROR:")                
                    print(f"Could not find PDB file input: {file_or_dir}")
                    print()
                    exit(1)
                pdb_file = file_or_dir
            elif flag == '-q':
                if not os.path.isfile(file_or_dir):
                    print("ERROR:")
                    print(f"Could not find Q-Force .itp file: {file_or_dir}")
                    print()
                    exit(1)
                itp_file = file_or_dir
            elif flag == '-a':
                if not os.path.isfile(file_or_dir):
                    print("ERROR:")
                    print(f"Could not find assignments file: {file_or_dir}")
                    print(f"polypal_ff can be used to generate the necessary assignment file")
                    print()
                    exit(1)
                assignment_file = file_or_dir

        if assignment_file == None:
            files = os.listdir()
            for file in files:
                if file.endswith('_assignments.txt'):
                    assignment_file = os.path.join(file)
                    break
            if assignment_file == None:
                print("ERROR:")
                print('Could not find assignments file. Please generate with polypal ff or specify location with -a flag if it not in current directory.')
                print()
                exit(1)

        if pdb_file == None:
            print("ERROR:")
            print('PDB file required for Assemble! topolgy generator')
            print('Use -h for help')
            print()
        if itp_file == None:
            print("ERROR:")
            print('polypal_top requires the itp file generated by Q-Force')
            print('Use -h for help')
            print()

        renamed_pdb = "".join(pdb_file.split('.')[:-1])+'_assemble.pdb'
        output_file = '.'+"".join(pdb_file.split('.')[:-1])+'_temp_top'
        final_output_file = "".join(pdb_file.split('.')[:-1])+'_top'
        start_bond = ';bonds'
        end_bond = ';angles'
        start_angle = ';angles'
        end_angle = ';ridgid dihedrals'
        start_dihedral = ';ridgid dihedrals'
        end_dihedral = ';flexible dihedrals'
        start_flex_diheds = ';flexible dihedrals'
        end_flex_diheds = '[ atomtypes ]'

        # Extract atom numbers from PDB file
        itp_atoms = read_itp_file(itp_file)
        pdb_atoms = read_pdb_file(pdb_file)
        full_list = find_matching_atoms(itp_atoms, pdb_atoms, output_file)

        #Match all parameters from force field to corresponding parameters
        matching_bonds = process_bonds(assignment_file, start_bond, end_bond, full_list)
        matching_angles = process_angles(assignment_file, start_angle, end_angle, full_list)
        matching_dihedrals = process_dihedrals(assignment_file, start_dihedral, end_dihedral, full_list)
        matching_flex_diheds = process_flex_diheds(assignment_file, start_flex_diheds, end_flex_diheds, full_list)

        # Write formatted lines to the output file
        write_formatted_bonds(matching_bonds, output_file)
        write_formatted_angles(matching_angles, output_file)
        write_formatted_dihedrals(matching_dihedrals, matching_flex_diheds, output_file)

        AA_name, ZZ_name = replace_AA_and_ZZ(output_file, final_output_file, pdb_file, renamed_pdb)
        write_terminal_types(output_file, final_output_file, AA_name, ZZ_name)

        print()
        print(f'Topolgy for Assemble! written to {final_output_file}')
        print(f'PDB file for Assemble! written to {renamed_pdb}')
        print()

        if os.path.isfile(output_file):
            os.remove(output_file)

 #   else:
 #       print("ERROR:")
 #       print("Could not perform conversion. Please see polypal pdb -h for help.")
 #       print()  

#if __name__ == "__main__":
#    main()
