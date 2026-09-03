import os
import re

def remove_slash_labels(line):
    """
    Removes internal node labels containing '/', by scanning for a right parenthesis ')'
    and then checking the substring up to the next boundary character (one of ':),;').

    Examples:
      - )0/47:0.123  -> ):0.123
      - )12.5/55);   -> ));      (if slash is found, remove label, keep boundary )
      - )12.5/55;    -> );       (if slash is found, remove label, keep boundary ;)
      - )label       -> )label   (if no slash, do nothing)
    """

    result = []
    i = 0
    length = len(line)

    # The set of boundary chars after a parenthesis-based label
    boundary_chars = {':', ')', ',', ';'}

    while i < length:
        char = line[i]

        if char == ')':
            # Look ahead to find the earliest occurrence of one of the boundary chars
            boundary_index = -1
            earliest_pos = None
            for bc in boundary_chars:
                pos = line.find(bc, i+1)
                if pos != -1:
                    if earliest_pos is None or pos < earliest_pos:
                        earliest_pos = pos

            if earliest_pos is not None:
                # Substring between ')' and that boundary
                between = line[i+1:earliest_pos]
                if '/' in between:
                    # Remove the entire label portion, but keep the boundary char
                    # We'll add ")" if boundary is ) or add ":" if boundary is :, etc.
                    # Actually, we keep just the initial ")" and then the boundary char
                    # Example:
                    #   )12.5/55) -> ))   (we keep both closing parentheses)
                    #   )12.5/55: -> ):   (keep the colon and what's after)
                    #   )12.5/55; -> );   (keep semicolon)
                    #   )12.5/55, -> ),   (keep comma)
                    boundary_char = line[earliest_pos]

                    # Insert a single ')'
                    result.append(')')

                    # Insert the boundary char (if it's not a second ')',
                    # we also need to handle the case )12.5/55) -> we want '))'
                    # So let's handle that logic:
                    if boundary_char in [':', ',', ';']:
                        # e.g. )12.5/55: -> ):   or )12.5/55; -> );
                        result.append(boundary_char)
                    else:
                        # boundary_char == ')'
                        # So )12.5/55) => we want '))'
                        result.append(')')
                    
                    # Now skip past the boundary char
                    i = earliest_pos + 1
                else:
                    # No slash => keep the entire segment as is
                    # We'll just add the character ')' and move on
                    result.append(')')
                    i += 1
            else:
                # No boundary found => just add ')' and continue
                result.append(')')
                i += 1
        else:
            result.append(char)
            i += 1

    return "".join(result)

def process_nexus_file(input_path, output_path):
    """
    1) Remove the second line if it contains '[R-package PHYTOOLS'.
    2) Remove everything from 'BEGIN TAXA;' up to (and including) the corresponding 'END;'.
    3) Remove node labels containing '/', handling both ')label/num:' and ')label/num)' (etc.).
    4) Write the cleaned lines to 'output_path'.
    """

    with open(input_path, 'r', encoding='utf-8') as infile:
        lines = infile.readlines()

    cleaned_lines = []
    skip_block = False
    line_count = 0
    taxa_block_removed = False

    for line in lines:
        line_count += 1

        # (1) Skip the second line if it contains '[R-package PHYTOOLS'
        if line_count == 2 and '[R-package PHYTOOLS' in line:
            continue

        # We'll convert to uppercase for block detection
        check_upper = line.strip().upper()

        # (2) Remove everything from 'BEGIN TAXA;' up to 'END;'
        if not taxa_block_removed:
            if check_upper.startswith("BEGIN TAXA;"):
                skip_block = True
                continue
            if skip_block:
                if check_upper.startswith("END;"):
                    skip_block = False
                    taxa_block_removed = True
                continue  # skip lines in the TAXA block

        # (3) For lines we keep, remove slash-based labels, including
        #     cases that end in a colon, a parenthesis, semicolon, etc.
        line_modified = remove_slash_labels(line)

        cleaned_lines.append(line_modified)

    # (4) Write out the final lines
    with open(output_path, 'w', encoding='utf-8') as outfile:
        outfile.writelines(cleaned_lines)

    print(f"Finished. Output written to: {output_path}")
