import os

# Directory containing the images
image_dir = "images"
output_html = "pdb_gallery.html"

# Get all PNG files in the directory
pdb_files = sorted([f for f in os.listdir(image_dir) if f.endswith(".png")])

# Generate HTML file
with open(output_html, "w") as html_file:
    html_file.write("<!DOCTYPE html>\n<html lang='en'>\n<head>\n")
    html_file.write("    <meta charset='UTF-8'>\n")
    html_file.write("    <meta name='viewport' content='width=device-width, initial-scale=1.0'>\n")
    html_file.write("    <title>PDB Image Gallery</title>\n")
    html_file.write("    <style>\n")
    html_file.write("        body { font-family: Arial, sans-serif; text-align: center; }\n")
    html_file.write("        table { width: 100%; border-collapse: collapse; }\n")
    html_file.write("        td { padding: 10px; border: 1px solid #ccc; text-align: center; }\n")
    html_file.write("        img { width: 400px; height: auto; display: block; margin: auto; }\n")
    html_file.write("    </style>\n</head>\n<body>\n")
    html_file.write("    <h1>PDB Image Gallery</h1>\n")
    html_file.write("    <table>\n")

    for index, pdb_file in enumerate(pdb_files):
        pdb_basename = os.path.splitext(pdb_file)[0]  # Get filename without extension
        if index % 4 == 0:
            html_file.write("        <tr>\n")
        html_file.write(f"            <td><img src='{image_dir}/{pdb_file}' alt='{pdb_file}'><br>{pdb_basename}</td>\n")
        if index % 4 == 3:
            html_file.write("        </tr>\n")
    
    html_file.write("    </table>\n")
    html_file.write("</body>\n</html>")

print(f"Gallery generated: {output_html}")
