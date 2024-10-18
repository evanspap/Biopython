import pandas as pd
import re
import sys

file_contents = open(sys.argv[1], "r")

# Initialize lists to store the extracted data
protein_names = []
pocket_counts = []
pocket_details = []

# Regular expressions to extract protein name, score, coordinates, and number of pockets
protein_name_pattern = re.compile(r'==> (.*?) <==')
score_coord_pattern = re.compile(r'score: ([\d\.]+) predicted center: \[([\d\.\s\-]+)\]')
pocket_count_pattern = re.compile(r'number of all predicted pockets: (\d+)')

current_protein = None
current_pockets = []

# Iterate over each line to extract the data
for line in file_contents:
    # Check for protein name
    protein_match = protein_name_pattern.search(line)
    if protein_match:
        # Save the previous protein data (if any) before moving to the next
        if current_protein and current_pockets:
            protein_names.append(current_protein)
            pocket_counts.append(len(current_pockets))
            pocket_details.append(current_pockets)
        
        # Start new protein entry
        current_protein = protein_match.group(1)
        current_pockets = []
    
    # Check for pocket score and coordinates
    score_coord_match = score_coord_pattern.search(line)
    if score_coord_match:
        score = float(score_coord_match.group(1))
        coords = [float(x) for x in score_coord_match.group(2).split()]
        current_pockets.append((score, coords))
    
    # Check for number of pockets (not strictly necessary as we're counting the pockets)
    pocket_count_match = pocket_count_pattern.search(line)

# Append the last protein data
if current_protein and current_pockets:
    protein_names.append(current_protein)
    pocket_counts.append(len(current_pockets))
    pocket_details.append(current_pockets)

# Create a dataframe from the extracted data
data = {
    'Protein Name': protein_names,
    'Number of Pockets': pocket_counts,
    'Pocket Details': pocket_details
}

df = pd.DataFrame(data)

#df.to_csv('druggability.csv', index=False)

#for pname in data['Protein Name']
#     print(data['Pocket Details (Score and Coordinates)'][data['Protein Name'].index(pname)][0][0])
for i in range(len(data['Protein Name'])):
    print(data['Protein Name'][i], end=" > ")
    print(data['Number of Pockets'][i], end=" > ")
    for j in range((data['Number of Pockets'][i])):
        print(data['Pocket Details'][i][j][0], end=" > ")
        #print(data['Pocket Details'][i][j][1], end=" > ")
    print()
        


#print(df)

# Display the first few rows of the dataframe
#df.head()
