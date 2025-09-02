#!/usr/bin/python3
def remove_lines(input_file, output_file, special_char):  
    with open(input_file, 'r', encoding='utf-8') as f_in:  
        lines = f_in.readlines()  
    with open(output_file, 'w', encoding='utf-8') as f_out: 
        for line in lines:
            line_content = str(line)
            column1 = str(input_file[:8])
            if special_char not in line:
                f_out.write(column1 +','+ line_content)

         
