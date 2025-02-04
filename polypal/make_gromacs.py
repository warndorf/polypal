import shutil
import sys
import os
import re

def process_charge_file(charge_dir, charge_flag):
    expected_charge_folders = {
        '-txt': 'txt',
        '-cif': 'cif',
        '-mol2': 'mol2',
        '-pqr': 'pqr',
    }
    for subfolder in expected_charge_folders.values():
        subfolder_path = os.path.join(charge_dir, subfolder)
        if not os.path.isdir(subfolder_path):
            program_header()
            print()
            print("ERROR:")
            print(f"Required subfolder not found: {subfolder}")
            print()
    target_folder = expected_charge_folders[charge_flag]
    target_folder_path = os.path.join(charge_dir, target_folder)
    target_file = None

    for filename in os.listdir(target_folder_path):
        file_path = os.path.join(target_folder_path, filename)
        if os.path.isfile(file_path):
            target_file = file_path
            break
        if not target_file:
            program_header()
            print()
            print("ERROR:")
            print(f"No file found in {target_folder_path}")
            print()
    return target_file

def read_charges_txt(txt_charges_file):
    with open(txt_charges_file, 'r') as file:
        file.readline() 
        charges_line = file.readline().strip()
        charges = charges_line.split() 
    return charges

def read_charges_cif(cif_charges_file):
    with open(cif_charges_file, 'r') as file:
        lines = file.readlines() 
    charges= []
    charges_section = False
    for line in lines:
        if line.strip().startswith('_partial_atomic_charges.charge'):
            charges_section = True
            continue
        elif line.strip().startswith('loop_') and charges_section:
            charges_section = False
            break
        if charges_section:
            parts = line.split()
            if len(parts) == 3:
                charge_value = line.split()[-1]
            charges.append(charge_value)
    return charges

def read_charges_pqr(pqr_charges_file):
    with open(pqr_charges_file, 'r') as file:
        lines = file.readlines() 
    charges= []
    for line in lines:
        parts = line.split()
        if len(parts) == 11:
            charge_value = line.split()[-2]
        charges.append(charge_value)
    return charges

def update_atoms_with_charges(charges, itp_output_file):
    with open(itp_output_file, 'r') as file:
        lines = file.readlines()
    updated_lines = []
    in_atoms_section = False
    charge_index = 0
    for line in lines:
        if line.strip().startswith('[ atoms ]'):
            in_atoms_section = True
            updated_lines.append(line)
            continue
        elif line.strip().startswith('[ bonds ]'):
            in_atoms_section = False
        if in_atoms_section and not line.strip().startswith(';'):
            parts = line.split()
            if len(parts) > 5 and charge_index < len(charges):
                parts[6] = charges[charge_index] 
                charge_index += 1
                updated_line = f"{parts[0]:>6} {parts[1]:<10} {parts[2]:>4} {parts[3]:<4} {parts[4]:<4} {parts[5]:>4} {parts[6]:>10} {parts[7]:>10}\n"
                updated_lines.append(updated_line)
            else:
                updated_lines.append(line)
        else:
            updated_lines.append(line)
    with open(itp_output_file, 'w') as file:
        file.writelines(updated_lines)

def fix_gro_file(gro_filename, gro_output, AA_input, ZZ_input):
    AA_name = []
    with open(gro_filename, 'r') as file:
        AA_gro_lines = file.readlines()
    with open(gro_output, 'w') as file:
        for line in AA_gro_lines:
            parts = line.split()
            if len(parts) >= 4:
                atom_name = line[12:16].strip()
                if atom_name == "AA":
                    line = line[:12] + f" {AA_input}0 " + line[16:]
                AA_name.append(AA_input+'0')
            file.write(line)            
    with open(gro_output, 'r') as file:
        ZZ_gro_lines = file.readlines()
    with open(gro_output, 'w') as file:
        AA_cant_match = AA_name[0]
        for line in ZZ_gro_lines:
            parts = line.split()
            if len(parts) >= 4:               
                atom_name = line[12:16].strip()
                if atom_name == "ZZ":
                    ZZ_change = ZZ_input+'0'
                    if ZZ_change == AA_cant_match:
                        line = line[:12] + f"{ZZ_input}00 " + line[16:]
                    elif ZZ_change != AA_cant_match:
                        line = line[:12] + f" {ZZ_input}0 " + line[16:]
                file.write(line)
    print()
    print(f'GRO file "{gro_filename}" updated to "{gro_output}"')
    print()

def fix_itp_file(itp_filename, itp_output_file, AA_input, ZZ_input):
    with open(itp_filename, 'r') as file:
        lines = file.readlines()
    updated_lines = []
    in_dihedral_section = False
    for line in lines:
        if line.startswith(' [ dihedrals ]'):
            in_dihedral_section = True
            updated_lines.append(line)
            continue
        elif line.startswith('#ifdef POSRES') or line.startswith(' [ impropers ]'):
            in_dihedral_section = False
        if in_dihedral_section and not line.strip().startswith(';'):
            parts = line.split()
            if len(parts) == 11:
                parts[4] = '3'
                updated_line = f"{parts[0]:>5} {parts[1]:>6} {parts[2]:>6} {parts[3]:>6} {parts[4]:>6} {parts[5]:>6} {parts[6]:>8} {parts[7]:>8} {parts[8]:>8} {parts[9]:>8} {parts[10]:>8}\n"
                updated_lines.append(updated_line)
            else:
                updated_lines.append(line)
        else:
            updated_lines.append(line)
    with open(itp_output_file, 'w') as file:
        file.writelines(updated_lines) 
    with open(itp_output_file, 'r') as file:
        lines = file.readlines()

    updated_lines2 = []
    in_angle_section = False
    for line in lines:
        if line.startswith(' [ angles ]'):
            in_angle_section = True
            updated_lines2.append(line)
            continue
        elif line.startswith(' [ dihedrals ]') or line.startswith(' [ impropers ]') or line.startswith('#ifdef POSRES'):
            in_angle_section = False
        if in_angle_section and not line.strip().startswith(';'):
            parts = line.split()
            if len(parts) == 6:
                parts[3] = '1'
                updated_line = f"{parts[0]:>5} {parts[1]:>6} {parts[2]:>6} {parts[3]:>6} {parts[4]:>6} {parts[5]:>6}\n"
                updated_lines2.append(updated_line)
            else:
                updated_lines2.append(line)
        else:
            updated_lines2.append(line)
    with open(itp_output_file, 'w') as file:
        file.writelines(updated_lines2)     

    AA_name = []
    with open(itp_output_file, 'r') as file:
        AA_itp_lines = file.readlines()
    with open(itp_output_file, 'w') as file:
        for line in AA_itp_lines:
            parts = line.split()
            if len(parts) == 8:
                atom_name = line[35:39].strip()
                if atom_name == "AA":
                    line = line[:35] + f" {AA_input}0 " + line[39:]
                AA_name.append(AA_input+'0')
            file.writelines(line)
                    
    with open(itp_output_file, 'r') as file:
        ZZ_itp_lines = file.readlines()
    with open(itp_output_file, 'w') as file:
        AA_cant_match = AA_name[0]
        for line in ZZ_itp_lines:
            parts = line.split()
            if len(parts) == 8:               
                atom_name = line[35:39].strip()
                if atom_name == "ZZ":
                    ZZ_change = ZZ_input+'0'
                    if ZZ_change == AA_cant_match:
                        line = line[:35] + f"{ZZ_input}00 " + line[39:] 
                    elif ZZ_change != AA_cant_match:
                        line = line[:35] + f" {ZZ_input}0 " + line[39:]
            file.writelines(line)
    print()

    
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


def fix_top_file(top_filename, qforce_itp, top_output):
    nonbond_params = []
    atom_type_params = []
    nonbond_section = False
    atom_type_section = False
    with open(qforce_itp, 'r') as file:
        for line in file:
            if line.strip().startswith('[ nonbond_params ]'):
                nonbond_section = True
                nonbond_params.append(line)
                continue
            elif line.strip().startswith('[ moleculetype ]'):
                nonbond_section = False
            if nonbond_section:
                nonbond_params.append(line)
    with open(top_filename, 'r') as file:
        for line in file:
            if line.strip().startswith('[ atomtypes ]'):
                atom_type_section = True
                atom_type_params.append(line)
                continue
            elif line.strip().startswith('#include'):
                atom_type_section = False
            if atom_type_section:
                incorrect_pattern = r'(\sA)(\d+\.\d+e[+-]\d+)'
                atom_type_params.append(re.sub(incorrect_pattern, r'\1 \2', line))
    with open(top_filename, 'r') as in_file, open(top_output, 'w') as out_file:
        for line in in_file:
            if line.strip().startswith('#include'):
                out_file.writelines(nonbond_params)
                new_include_line = line.replace('.itp"', '_final.itp"')
                out_file.write(new_include_line)
            elif line.strip().startswith('; generated with'):
                new_title_line = line.replace('2018', '2018, and updated by PolyPal')
                out_file.write(new_title_line)
            else:
                out_file.write(line)
    temp_output = top_output + '.temp'
    with open(top_output, 'r') as out_file, open(temp_output, 'w') as file:
        atom_section = False
        for line in out_file:
            if line.strip().startswith('[ atomtypes ]'):
                atom_section = True
                #file.write('[ atomtypes ]\n')  # Keep the section header
                continue
            if atom_section and line.strip().startswith('[ nonbond_params ]'):
                atom_section = False
                file.writelines(atom_type_params)  # Add the new section content
                file.write('[ nonbond_params ]\n')  # Keep the section footer
                continue
            if not atom_section:
                file.write(line)
    shutil.move(temp_output, top_output)

def make_gromacs(command=None):
    charge_flag = None
    charge_dir = None
    itp_filename = None
    top_filename = None
    gro_filename = None
    pdb_filename = None
    qforce_itp = None

    if command:

        if len(sys.argv) < 1:
            print("ERROR:")
            print("Could not perform conversion. Please see polypal gmx -h for help.")
            print() 
            exit(1)

        if command[0].startswith('-') and command[0] != '-c' and command[0] != '-f' and command[0] != '-g' and command[0] != '-t' and command[0] != '-p' and command[0] != '-q':
            charge_flag = command[0]
            start_index = 1 
        else:
            start_index = 0  

        for i in range(start_index, len(command), 2):
            flag = command[i]
            file_or_dir = command[i + 1]
            if flag == '-c':
                if not os.path.isdir(file_or_dir):
                    print()
                    print("ERROR:")
                    print(f"could not find: {file_or_dir}")
                    print()
                    exit(1)
                charge_dir = file_or_dir
            elif flag == '-f':
                if not os.path.isfile(file_or_dir):
                    print()
                    print("ERROR:")
                    print(f"could not find: {file_or_dir}")
                    print()
                    exit(1)
                itp_filename = file_or_dir
            elif flag == '-g':
                if not os.path.isfile(file_or_dir):
                    print()
                    print("ERROR:")
                    print(f"could not find: {file_or_dir}")
                    print()
                    exit(1)
                gro_filename = file_or_dir
            elif flag == '-p':
                if not os.path.isfile(file_or_dir):
                    print()
                    print("ERROR:")
                    print(f"could not find: {file_or_dir}")
                    print()
                    exit(1)
                pdb_filename = file_or_dir
                output_pdb_filename = "".join(pdb_filename.split('.')[:-1])+'_final.pdb'
            elif flag == '-t':
                if not os.path.isfile(file_or_dir):
                    print()
                    print("ERROR:")
                    print(f"could not find: {file_or_dir}")
                    print()
                    exit(1)
                top_filename = file_or_dir
            elif flag == '-q':
                if not os.path.isfile(file_or_dir):
                    print()
                    print("ERROR:")
                    print(f"could not find: {file_or_dir}")
                    print()
                    exit(1)
                qforce_itp = file_or_dir           

        if itp_filename or pdb_filename or gro_filename:
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
            
        if pdb_filename:
            AA_to_ZZ_pdb(pdb_filename, output_pdb_filename, AA_input, ZZ_input)

        if gro_filename:
            gro_output = "".join(gro_filename.split('.')[:-1])+'_final.gro'
            fix_gro_file(gro_filename, gro_output, AA_input, ZZ_input)

        if top_filename:
            if qforce_itp == None:
                print("ERROR:")
                print('To fix Assemble! top file for Gromacs, please provide the itp forcefield file created by Q-Force')
                print()
            top_output = "".join(top_filename.split('.')[:-1])+'_final.top'
            fix_top_file(top_filename, qforce_itp, top_output)
            print(f'Updated topology for Gromacs written to: {top_output}')

        if charge_dir and itp_filename:
            if charge_flag:
                target_file = process_charge_file(charge_dir, charge_flag)
                if charge_flag == '-txt':
                    charges = read_charges_txt(target_file)
                elif charge_flag == '-pqr':
                    charges = read_charges_pqr(target_file)
                elif charge_flag == '-mol2':
                    print("ERROR:")
                    print("ACC II did not provide file for mol2, please use other option")
                    print()
                elif charge_flag == '-cif':
                    charges = read_charges_cif(target_file)
            elif charge_flag == None:
                charge_flag = '-txt'
                target_file = process_charge_file(charge_dir, charge_flag)
                charges = read_charges_txt(target_file)

            itp_output_file = "".join(itp_filename.split('.')[:-1])+'_final.itp'
            fix_itp_file(itp_filename, itp_output_file, AA_input, ZZ_input)
            update_atoms_with_charges(charges, itp_output_file)
            print(f'ITP file "{itp_filename}" updated with charges from directory {charge_dir}, and saved to "{itp_output_file}"')
            print()

        elif charge_dir == None and itp_filename:
            itp_output_file = "".join(itp_filename.split('.')[:-1])+'_final.itp'
            fix_itp_file(itp_filename, itp_output_file, AA_input, ZZ_input)
            print(f'No charges were provides for "{itp_filename}". The ITP file with no charges was saved to "{itp_output_file}"')
            print()
        
        elif charge_dir and itp_filename == None:
            print(f'The charge directory {charge_dir} was provided, but no ITP file was provided. Please provide one with -f if you would like to update charges')
            print()

