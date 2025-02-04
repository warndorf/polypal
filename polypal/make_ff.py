import re
import os 
import pkg_resources

def extract_atom_names(itp_filename):
    atom_names = []
    with open(itp_filename, 'r') as file:
        for line in file:
            if line.startswith('[ atoms ]'):
                break
        next(file)
        for line in file:
            if line.startswith('['):
                break
            parts = line.split()
            if len(parts) > 4:
                index = int(parts[0])
                atom_name = parts[4]
                atom_type = parts[1]
                atom_names.append((index, atom_name, atom_type))
    return atom_names

def extract_atomtypes(ext_lj_file, nonbonded_itp):
    with open(ext_lj_file, 'r') as ext_lj:
        names = set(line.strip() for line in ext_lj)
    atom_types = set()
    in_atomtypes_section = False
    with open(nonbonded_itp, 'r') as nonbond_ff:
        for line in nonbond_ff:
            stripped_line = line.strip()
            if stripped_line == '[ atomtypes ]':
                in_atomtypes_section = True
                continue
            if in_atomtypes_section and stripped_line == '':
                break
            if in_atomtypes_section:
                tokens = stripped_line.split()
                if tokens and tokens[0] in names:
                    atom_types.add(line)
    return list(atom_types)

def bonds_in_itp(itp_filename):
    bonds = []
    bond_data = []
    with open(itp_filename, 'r') as file:
        for line in file:
            if line.startswith('[ bonds ]'):
                break
        next(file)  # Skip header line
        for line in file:
            if line.startswith('['):
                break
            parts = line.split()
            if len(parts) >= 5:
                atom1 = parts[0]
                atom2 = parts[1]
                bond_type = parts[2]
                bond_length = parts[3]
                bond_force_const = parts[4]
                bonds.append((atom1, atom2, bond_type, bond_length, bond_force_const))
    for i, bond in enumerate(bonds, start=1):
        bond_name = f"b{i}"
        bond_data.append((bond_name, *bond))
    return bond_data

def angles_in_itp(itp_filename):
    angles = []
    angle_data =[]
    with open(itp_filename, 'r') as file:
        for line in file:
            if line.startswith('[ angles ]'):
                break
        next(file)  # Skip header line
        for line in file:
            if line.startswith('['):
                break
            parts = line.split()
            if len(parts) == 8:
                atom1 = parts[0]
                atom2 = parts[1]
                atom3 = parts[2]
                angle_type = parts[3]
                angle_theta0 = parts[4]
                angle_k_theta = parts[5]
                r0 = parts[6]
                k_bond = parts[7]
                angles.append((atom1, atom2, atom3, angle_type, angle_theta0, angle_k_theta, r0, k_bond))
            if len(parts) == 6:
                atom1 = parts[0]
                atom2 = parts[1]
                atom3 = parts[2]
                angle_type = parts[3]
                r0 = ('   ')
                k_bond = ('   ')
                angles.append((atom1, atom2, atom3, angle_type, angle_theta0, angle_k_theta, r0, k_bond))                
    for i, angle in enumerate(angles, start=1):
        angle_name = f"a{i}"
        angle_data.append((angle_name, *angle))
    return angle_data

def dihedrals_in_itp(itp_filename):
    ridgid_dihedrals = []
    ridgid_dihed_data =[]
    with open(itp_filename, 'r') as file:
        for line in file:
            if line.startswith('; rigid dihedrals '):
                break
        for line in file:
            if line.startswith('[ pairs ]') or line.startswith('; flexible dihedrals') or line.startswith('; improper dihedrals') or line.startswith('; inversion dihedrals'):
                break
            parts = line.split()
            if len(parts) >= 5 and parts[0].isdigit():
                atom1 = parts[0]
                atom2 = parts[1]
                atom3 = parts[2]
                atom4 = parts[3]
                dihed_type = parts[4]
                dihed_theta0 = parts[5]
                dihed_k_theta = parts[6]
                ridgid_dihedrals.append((atom1, atom2, atom3, atom4, dihed_type, dihed_theta0, dihed_k_theta))
    for i, dihed in enumerate(ridgid_dihedrals, start=1):
        dihed_name = f"d{i}"
        ridgid_dihed_data.append((dihed_name, *dihed))
    return ridgid_dihed_data

def flex_diheds_in_itp(itp_filename):
    flex_dihedrals = []
    flex_dihed_data =[]
    with open(itp_filename, 'r') as file:
        for line in file:
            if line.startswith('; flexible dihedrals '):
                break
        for line in file:
            if line.startswith('[ pairs ]') or line.startswith('; improper dihedrals') or line.startswith('; rigid dihedrals') or line.startswith('; inversion dihedrals'):
                break
            parts = line.split()
            if len(parts) >= 5 and parts[0].isdigit():
                atom1 = parts[0]
                atom2 = parts[1]
                atom3 = parts[2]
                atom4 = parts[3]
                dihed_type = parts[4]
                c0 = parts[5]
                c1 = parts[6]
                c2 = parts[7]
                c3 = parts[8]
                c4 = parts[9]
                c5 = parts[10]
                flex_dihedrals.append((atom1, atom2, atom3, atom4, dihed_type, c0, c1, c2, c3, c4, c5))
    for i, dihed in enumerate(flex_dihedrals, start=1):
        dihed_name = f"f{i}"
        flex_dihed_data.append((dihed_name, *dihed))
    return flex_dihed_data

def impropers_in_itp(itp_filename):
    impropers = []
    impropers_data =[]
    with open(itp_filename, 'r') as file:
        for line in file:
            if line.startswith('; improper dihedrals'):
                break
        for line in file:
            if line.startswith('[ pairs ]') or line.startswith('; flexible dihedrals') or line.startswith('; rigid dihedrals') or line.startswith('; inversion dihedrals'):
                break
            parts = line.split()
            if len(parts) == 7 and parts[0].isdigit():
                atom1 = parts[0]
                atom2 = parts[1]
                atom3 = parts[2]
                atom4 = parts[3]
                dihed_type = parts[4]
                dihed_theta0 = parts[5]
                dihed_k_theta = parts[6]
                impropers.append((atom1, atom2, atom3, atom4, dihed_type, dihed_theta0, dihed_k_theta))
    for i, improper in enumerate(impropers, start=1):
        improper_name = f"i{i}"
        impropers_data.append((improper_name, *improper))
    return impropers_data

def inversions_in_itp(itp_filename):
    inversion_dihedrals = []
    inversions_data =[]
    with open(itp_filename, 'r') as file:
        for line in file:
            if line.startswith('; inversion dihedrals'):
                break
        for line in file:
            if line.startswith('[ pairs ]') or line.startswith('; flexible dihedrals') or line.startswith('; rigid dihedrals') or line.startswith('; improper dihedrals'):
                break
            parts = line.split()
            if len(parts) >= 5 and parts[0].isdigit():
                atom1 = parts[0]
                atom2 = parts[1]
                atom3 = parts[2]
                atom4 = parts[3]
                dihed_type = parts[4]
                c0 = parts[5]
                c1 = parts[6]
                c2 = parts[7]
                c3 = parts[8]
                c4 = parts[9]
                c5 = parts[10]
                inversion_dihedrals.append((atom1, atom2, atom3, atom4, dihed_type, c0, c1, c2, c3, c4, c5))
    for i, inversion in enumerate(inversion_dihedrals, start=1):
        inversion_name = f"v{i}"
        inversions_data.append((inversion_name, *inversion))
    return inversions_data


def write_data_to_file(bond_data, angle_data, ridgid_dihed_data, flex_dihed_data, impropers_data, inversions_data, atom_names, atom_types, temp_output_filename):
    with open(temp_output_filename, 'w') as file:
        bonded_types_header = "{:<6} {:<6} {:<7}  {:<7}\n".format(";bonds","angles", "dihedrals", "impropers")
        bonded_types_info =  "  {:<6} {:<6} {:<7}  {:<7}\n".format("1","5", "2", "2")
        file.write(" [ bondedtypes ]\n")
        file.write(bonded_types_header)
        file.write(bonded_types_info)
    with open(temp_output_filename, 'a') as file:
        atom_name1 = []
        atom_name2 = []
        bond_header = "{:<6}{:<10} {:<7}\n".format(";Name","r0", "k_bond")
        file.write(";"+"=" * len(bond_header) + "\n")
        file.write(bond_header)
        file.write(";bonds" + "\n")
        for bond in bond_data:
            bond_name, atom1, atom2, bond_type, bond_length, bond_force_constant = bond
            for name_info in atom_names:
                index, atom_name, atom_type = name_info
                if str(atom1) == str(index):
                    atom_name1 = atom_name
                    break
            for name_info in atom_names:
                index, atom_name, atom_type = name_info
                if str(atom2) == str(index):
                    atom_name2 = atom_name
                    break
            file.write("{:<5} {:<10} {:<10} {:<0}{:<3} {:<5} {:<5} {:<5}\n".format(bond_name, bond_length, bond_force_constant, ";  ", atom1, atom2, str(atom_name1), str(atom_name2)))
        file.write("\n")
    with open(temp_output_filename, 'a') as file:
        atom_name1 = []
        atom_name2 = []
        atom_name3 = []
        angle_header = "{:<5} {:<10}{:<10} {:<10}  {:<10}\n".format(";Name", "theta0", "k_theta", "r0", "k_bond")
        file.write(";"+"=" * len(angle_header) + "\n")
        file.write(angle_header)
        file.write(";angles" + "\n")
        for angle in angle_data:
            angle_name, atom1, atom2, atom3, angle_type, angle_theta0, angle_k_theta, r0, k_bond = angle
            for name_info in atom_names:
                index, atom_name, atom_type = name_info
                if str(atom1) == str(index):
                    atom_name1 = atom_name
                    break
            for name_info in atom_names:
                index, atom_name, atom_type = name_info
                if str(atom2) == str(index):
                    atom_name2 = atom_name
                    break
            for name_info in atom_names:
                index, atom_name, atom_type = name_info
                if str(atom3) == str(index):
                    atom_name3 = atom_name
                    break             
            file.write("{:<5} {:<10} {:<10} {:<10} {:<10} {:<0}{:<5}{:<5}{:<5}{:<5}{:<5}{:<5}\n".format(angle_name, angle_theta0, angle_k_theta, r0, k_bond,";  ", atom1, atom2, atom3, str(atom_name1), str(atom_name2), str(atom_name3)))
        file.write("\n")
    with open(temp_output_filename, 'a') as file:
        atom_name1 = []
        atom_name2 = []
        atom_name3 = []
        atom_name4 = []
        ridgid_dihed_header = "{:<5} {:<7}   {:<10}\n".format(";Name", "theta0", "k_theta")
        file.write(";"+"=" * len(ridgid_dihed_header) + "\n")
        file.write(ridgid_dihed_header)  # Underline the headers
        file.write(";ridgid dihedrals" + "\n")
        for diheds in ridgid_dihed_data:
            dihed_name, atom1, atom2, atom3, atom4, dihed_type, dihed_theta0, dihed_k_theta = diheds
            for name_info in atom_names:
                index, atom_name, atom_type = name_info
                if str(atom1) == str(index):
                    atom_name1 = atom_name
                    break
            for name_info in atom_names:
                index, atom_name, atom_type = name_info
                if str(atom2) == str(index):
                    atom_name2 = atom_name
                    break
            for name_info in atom_names:
                index, atom_name, atom_type = name_info
                if str(atom3) == str(index):
                    atom_name3 = atom_name
                    break   
            for name_info in atom_names:
                index, atom_name, atom_type = name_info
                if str(atom4) == str(index):
                    atom_name4 = atom_name
                    break                          
            file.write("{:<5} {:<10} {:<10} {:<0}{:<5} {:<5}{:<5}{:<5} {:<5} {:<5}{:<5}{:<5}\n".format(dihed_name, dihed_theta0, dihed_k_theta,";  ", atom1, atom2, atom3, atom4, str(atom_name1), str(atom_name2), str(atom_name3), str(atom_name4)))
        file.write("\n")
    with open(temp_output_filename, 'a') as file:
        impropers_header = "{:<5} {:<7}   {:<10}\n".format(";Name", "theta0", "k_theta")
        file.write(";"+"=" * len(impropers_header) + "\n")
        file.write(impropers_header)  # Underline the headers
        file.write(";improper dihedrals" + "\n")
        for impropers in impropers_data:
            improper_name, atom1, atom2, atom3, atom4, dihed_type, dihed_theta0, dihed_k_theta = impropers
            for name_info in atom_names:
                index, atom_name, atom_type = name_info
                if str(atom1) == str(index):
                    atom_name1 = atom_name
                    break
            for name_info in atom_names:
                index, atom_name, atom_type = name_info
                if str(atom2) == str(index):
                    atom_name2 = atom_name
                    break
            for name_info in atom_names:
                index, atom_name, atom_type = name_info
                if str(atom3) == str(index):
                    atom_name3 = atom_name
                    break   
            for name_info in atom_names:
                index, atom_name, atom_type = name_info
                if str(atom4) == str(index):
                    atom_name4 = atom_name
                    break  
            file.write("{:<5} {:<10} {:<10} {:<0}{:<5} {:<5}{:<5}{:<5} {:<5} {:<5}{:<5}{:<5}\n".format(improper_name, dihed_theta0, dihed_k_theta,";  ", atom1, atom2, atom3, atom4, str(atom_name1), str(atom_name2), str(atom_name3), str(atom_name4)))
        file.write("\n")  
    with open(temp_output_filename, 'a') as file:
        flex_dihed_header = "{:<0}{:<5}{:<7}{:<10}{:<10}{:<10}{:<10}{:<6}\n".format(";","name","c0","c1","c2","c3","c4","c5")
        file.write(";"+"=" * len(flex_dihed_header) + "\n")
        file.write(flex_dihed_header)
        file.write(";flexible dihedrals" + "\n")
        for flex_diheds in flex_dihed_data:
            dihed_name, atom1, atom2, atom3, atom4, dihed_type, c0, c1, c2, c3, c4, c5 = flex_diheds
            for name_info in atom_names:
                index, atom_name, atom_type = name_info
                if str(atom1) == str(index):
                    atom_name1 = atom_name
                    break
            for name_info in atom_names:
                index, atom_name, atom_type = name_info
                if str(atom2) == str(index):
                    atom_name2 = atom_name
                    break
            for name_info in atom_names:
                index, atom_name, atom_type = name_info
                if str(atom3) == str(index):
                    atom_name3 = atom_name
                    break   
            for name_info in atom_names:
                index, atom_name, atom_type = name_info
                if str(atom4) == str(index):
                    atom_name4 = atom_name
                    break 
            file.write("{:<5} {:<7} {:<10} {:<10} {:<10}{:<10}{:<10}{:<0}{:<5}{:<5}{:<5}{:<5}{:<5}{:<5}{:<5}{:<5}\n".format(dihed_name, c0, c1, c2, c3, c4, c5,";  ", atom1, atom2, atom3, atom4, str(atom_name1), str(atom_name2), str(atom_name3), str(atom_name4)))
        file.write("\n")
    with open(temp_output_filename, 'a') as file:
        inversion_dihed_header = "{:<0}{:<5}{:<7}{:<10}{:<10}{:<10}{:<10}{:<6}\n".format(";","name","c0","c1","c2","c3","c4","c5")
        file.write(";"+"=" * len(inversion_dihed_header) + "\n")
        file.write(inversion_dihed_header)
        file.write(";inversion dihedrals" + "\n")
        for inversion_diheds in inversions_data:
            inversion_name, atom1, atom2, atom3, atom4, dihed_type, c0, c1, c2, c3, c4, c5 = inversion_diheds
            for name_info in atom_names:
                index, atom_name, atom_type = name_info
                if str(atom1) == str(index):
                    atom_name1 = atom_name
                    break
            for name_info in atom_names:
                index, atom_name, atom_type = name_info
                if str(atom2) == str(index):
                    atom_name2 = atom_name
                    break
            for name_info in atom_names:
                index, atom_name, atom_type = name_info
                if str(atom3) == str(index):
                    atom_name3 = atom_name
                    break   
            for name_info in atom_names:
                index, atom_name, atom_type = name_info
                if str(atom4) == str(index):
                    atom_name4 = atom_name
                    break 
            file.write("{:<5} {:<7} {:<10} {:<10} {:<10}{:<10}{:<10}{:<0}{:<5}{:<5}{:<5}{:<5}{:<5}{:<5}{:<5}{:<5}\n".format(inversion_name, c0, c1, c2, c3, c4, c5,";  ", atom1, atom2, atom3, atom4, str(atom_name1), str(atom_name2), str(atom_name3), str(atom_name4)))
        file.write("\n")
    with open(temp_output_filename, 'a') as file:
        lj_header = "{:<0}{:<10}{:<10}   {:<10}  {:<10}{:<10}    {:<10}{:<6}\n".format(";","name","at_num","mass","charge","type","sigma","epsilon")
        file.write(";"+"=" * 60 + "\n")
        file.write("[ atomtypes ]"+ "\n")
        file.write(lj_header)
        column_width=12
        for line in atom_types:
            columns = line.split()
            if len(columns) > 1:
                columns.pop(1)  # Remove the second column
            final_lj_line = [word.ljust(column_width) for word in columns]
            file.write(''.join(final_lj_line)+ '\n')
    with open(temp_output_filename, 'a') as file:
        file.write("\n" + ";"+"=" * 60 + "\n")
        file.write("[ defaults ]"+ "\n")
        combo_header = "{:<0}{:<10}{:<10}{:<10}{:<10}{:<10}\n".format(";", "nbfunc", "comb-rule", "gen-pairs", "fudgeLJ", "fudgeQQ")
        combo_rule = "{:<10}{:<10}{:<10}{:<10}{:<10}\n".format("1","3","yes","0.5","0.5")
        file.write(combo_header)
        file.write(combo_rule)

def final_file_processing(temp_output_filename, output_filename):
    with open(temp_output_filename, 'r') as file:
        lines = file.readlines()
    with open(output_filename, 'w') as output_file:
        for line in lines:
            if re.search(r';\s+(\d+\s+)+([A-Za-z0-9]+\s+)+$', line):
                new_line = re.sub(r';\s+(\d+\s+)+([A-Za-z0-9]+\s+)+$', '', line).rstrip() + '\n'
                output_file.write(new_line)
            else:
                output_file.write(line)

def make_ff(command=None):
    itp_filename = None
    ext_lj_file = None
    nonbonded_itp = None    


    if command:
        if len(command) < 1:
            print("ERROR:")
            print("Could not perform conversion. Please see polypal ff -h for help.")
            print() 
            exit(1)

        if command[0].startswith('-') and command[0] != '-q' and command[0] != '-lj' and command[0] != '-n':
            start_index = 1
        else:
            start_index = 0 

        for i in range(start_index, len(command), 2):
            flag = command[i]
            file_or_dir = command[i + 1]
            if flag == '-q':
                if not os.path.isfile(file_or_dir):
                    print("ERROR:")
                    print(f"Could not find Q-Force .itp file: {file_or_dir}")
                    print()
                    exit(1)
                itp_filename = file_or_dir
                temp_output_filename = itp_filename.replace('.itp', '_assignments.txt')
                output_filename = itp_filename.replace('.itp', '_ff.txt')
            elif flag == '-lj':
                if not os.path.isfile(file_or_dir):
                    print("ERROR:")
                    print(f"Could not find ext_lj file: {file_or_dir}")
                    print()
                    exit(1)
                ext_lj_file = file_or_dir
            elif flag == '-n':
                if not os.path.isfile(file_or_dir):
                    print("ERROR:")
                    print(f"Could not find OPLS non-bonded .itp file: {file_or_dir}")
                    print()
                    exit(1)
                nonbonded_itp = file_or_dir        
    

        if ext_lj_file == None:
            ext_lj_file = 'ext_lj'
            if not os.path.isfile(ext_lj_file):
                print("ERROR:")
                print(f"Could not find ext_lj, please specify with -lj")
                print()
                exit(1)

        if nonbonded_itp == None:
            docs_folder = os.path.join(os.path.dirname(pkg_resources.resource_filename('polypal', '__init__.py')), 'data')
            nonbonded_itp = os.path.join(docs_folder, 'opls_ffnonbonded.itp')
                
        if itp_filename and ext_lj_file:
            bond_data = bonds_in_itp(itp_filename)
            angle_data = angles_in_itp(itp_filename)
            ridgid_dihed_data = dihedrals_in_itp(itp_filename)
            flex_dihed_data = flex_diheds_in_itp(itp_filename)
            impropers_data=impropers_in_itp(itp_filename)
            inversions_data=inversions_in_itp(itp_filename)
            atom_types = extract_atomtypes(ext_lj_file, nonbonded_itp)
            atom_names = extract_atom_names(itp_filename)
            write_data_to_file(bond_data, angle_data,ridgid_dihed_data, flex_dihed_data, impropers_data, inversions_data, atom_names, atom_types, temp_output_filename)
            final_file_processing(temp_output_filename, output_filename)

            print(f"FF file for use in Assemble! written to {os.path.basename(output_filename)}")
            print(f"Atom assignments written to {os.path.basename(temp_output_filename)}")
            print()


    #    else:
    #        print("ERROR:")
    #        print("Could not perform conversion. Please see polypal pdb -h for help.")
    #        print()  

#if __name__ == "__main__":
#    main()
