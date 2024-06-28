import os
import numpy as np
from decimal import Decimal

def format_e(n):
    a = '%E' % n
    return a.split('E')[0].rstrip('0').rstrip('.') + 'E' + a.split('E')[1]



def adjust_variables(base_file, freq_value, dt_0_value, output_name):
    # Read the base file
    with open(base_file, 'r') as f:
        lines = f.readlines()

    # Adjust the frequency and time step
    freq = format_e(Decimal(str(freq_value)))
    dt =  format_e(Decimal(dt_0_value))
    freq_str = format_e(Decimal(str(freq_value))).replace('+', '')

    # Modify the lines containing the variables
    for i, line in enumerate(lines):
        if "vals = '0.02 w'" in line:
            lines[i] = f"    vals = '0.02 {freq_str}'\n"
        elif 'dt = dti' in line:
            lines[i] = f"    dt = {dt}\n"
        elif 'file_base = outputs' in line:
            lines[i] = f"   file_base = out_aBTO_8-200nm_0.5_T298K_VzEx_E2e-2_w{freq_str}\n"

    # Write the modified content to a new file
    output_file = f"{output_name}.i"
    with open(output_file, 'w') as f:
        f.writelines(lines)

    print(f"File '{output_file}' has been created.")

# Example usage
base_file = 'induce_aBTO_8-120nm_0.5_T298K_VzEx_E2e-2_INPUT.i'  # Update with your base file name
dt_0 = 1.0e-10
#frequencies = [1 ,1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 6.5, 7, 7.5, 8, 8.5, 9, 9.5]
frequencies = np.arange(1, 10, 0.25)
scales = [1e+08, 1e+09, 1e+10, 1e+11, 1e+12, 1e+13]

for j in scales:
    for i in frequencies:
        freq_value = i*j  # Update with the desired frequency value
        dt_value = dt_0 / (freq_value/1e9)  # Update with the desired time step value
        output_name = 'induce_aBTO_8-120nm_0.5_T298K_VzEx_E2e-2_w'+format_e(Decimal(str(freq_value))).replace('+', '').replace('.', '_')  # Update with the desired output file name
        adjust_variables(base_file, freq_value, dt_value, output_name)
