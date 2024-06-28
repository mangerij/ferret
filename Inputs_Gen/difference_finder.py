def compare_files(file1_path, file2_path):
    with open(file1_path, 'r') as file1:
        lines1 = file1.readlines()
    with open(file2_path, 'r') as file2:
        lines2 = file2.readlines()

    num_lines1 = len(lines1)
    num_lines2 = len(lines2)

    # if num_lines1 != num_lines2:
    #     print("Files have different number of lines.")
    #     return

    for i in range(num_lines1):
        if lines1[i] != lines2[i]:
            print(f"Difference found at line {i + 1}:")
            print(f"File 1: {lines1[i]}")
            print(f"File 2: {lines2[i]}")
            return
# Example usage:
file1_path = 'induce_aBTO_8-120nm_0.5_T298K_VzEx_E2e-2_INPUT.i'
file2_path = 'induce_aBTO_8-120nm_0.5_T298K_VzEx_E2e-2_w3E08.i'
compare_files(file1_path, file2_path)
