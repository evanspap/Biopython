#!/usr/bin/env python3
"""
Generate_Pdb_Image_List_clickable.py

This script generates an HTML file displaying a gallery of PNG images from a specified directory.
Each image can be clicked to open a full-size version in the same browser tab with the filename
displayed below, and it returns to the gallery at the same scroll position when going back.

Usage:
    python Generate_Pdb_Image_List_clickable.py [image_dir] [output_html]

Arguments:
    image_dir   - Directory containing PNG images (default: "images")
    output_html - Output HTML file (default: "pdb_gallery.html")

Example:
    python Generate_Pdb_Image_List_clickable.py pdb_images gallery.html
"""

import os
import sys

# Default values
image_dir = "images"
output_html = "pdb_gallery.html"

# Get arguments if provided
if len(sys.argv) > 1:
    image_dir = sys.argv[1]
if len(sys.argv) > 2:
    output_html = sys.argv[2]

# Check if the directory exists
if not os.path.isdir(image_dir):
    print(f"Error: Directory '{image_dir}' does not exist.")
    sys.exit(1)

# Get all PNG files in the directory
pdb_files = sorted([f for f in os.listdir(image_dir) if f.endswith(".png")])

# Generate the HTML
with open(output_html, "w") as html_file:
    html_file.write("<!DOCTYPE html>\n")
    html_file.write("<html lang='en'>\n<head>\n")
    html_file.write("    <meta charset='UTF-8'>\n")
    html_file.write("    <meta name='viewport' content='width=device-width, initial-scale=1.0'>\n")
    html_file.write("    <title>PDB Image Gallery</title>\n")

    # Some basic CSS
    html_file.write("    <style>\n")
    html_file.write("        body { font-family: Arial, sans-serif; margin: 0; padding: 0; text-align: center; }\n")
    html_file.write("        .header { position: fixed; top: 0; width: 100%; background: white; padding: 10px 0; "
                    "box-shadow: 0px 4px 2px -2px gray; z-index: 1000; }\n")
    html_file.write("        .content { margin-top: 60px; }\n")
    html_file.write("        table { width: 100%; border-collapse: collapse; }\n")
    html_file.write("        td { padding: 10px; border: 1px solid #ccc; text-align: center; vertical-align: middle; }\n")
    html_file.write("        img { width: 400px; height: auto; display: block; margin: auto; cursor: pointer; }\n")
    html_file.write("        /* Overlay for full image view */\n")
    html_file.write("        #overlay {\n")
    html_file.write("            position: fixed;\n")
    html_file.write("            top: 0; left: 0; right: 0; bottom: 0;\n")
    html_file.write("            background: white;\n")
    html_file.write("            display: none;\n")
    html_file.write("            z-index: 2000;\n")
    html_file.write("            overflow: auto;\n")
    html_file.write("        }\n")
    html_file.write("        #overlay-content {\n")
    html_file.write("            max-width: 1200px;\n")
    html_file.write("            margin: 40px auto;\n")
    html_file.write("        }\n")
    html_file.write("        #overlay img {\n")
    html_file.write("            width: 100%;\n")
    html_file.write("            max-width: 1200px;\n")
    html_file.write("            height: auto;\n")
    html_file.write("        }\n")
    html_file.write("        #backButton {\n")
    html_file.write("            padding: 10px 20px;\n")
    html_file.write("            font-size: 16px;\n")
    html_file.write("            margin-top: 20px;\n")
    html_file.write("        }\n")
    html_file.write("    </style>\n</head>\n")

    html_file.write("<body>\n")
    html_file.write("    <div class='header'><h1>PDB Image Gallery</h1></div>\n")

    # The gallery container
    html_file.write("    <div class='content' id='gallery'>\n")
    html_file.write("    <table>\n")

    for i, pdb_file in enumerate(pdb_files):
        pdb_basename = os.path.splitext(pdb_file)[0]
        if i % 4 == 0:
            html_file.write("        <tr>\n")

        html_file.write("            <td>\n")
        # Each image calls showOverlay() on click
        html_file.write(f"                <img src='{image_dir}/{pdb_file}' alt='{pdb_file}' "
                        f"onclick=\"showOverlay('{image_dir}/{pdb_file}','{pdb_basename}')\" >\n")
        html_file.write(f"                <br>{pdb_basename}\n")
        html_file.write("            </td>\n")

        if i % 4 == 3:
            html_file.write("        </tr>\n")

    # Close the last row if needed
    if len(pdb_files) % 4 != 0:
        html_file.write("        </tr>\n")

    html_file.write("    </table>\n")
    html_file.write("    </div>\n")

    # Overlay for showing the big image
    html_file.write("    <div id='overlay'>\n")
    html_file.write("        <div id='overlay-content'>\n")
    html_file.write("            <img id='overlay-img' src='' alt='Full Image'>\n")
    html_file.write("            <p id='overlay-caption' style='text-align: center; font-size: 20px;'></p>\n")
    html_file.write("            <button id='backButton' onclick='hideOverlay()'>Back</button>\n")
    html_file.write("        </div>\n")
    html_file.write("    </div>\n")

    # Script for the overlay & history
    html_file.write("    <script>\n")
    html_file.write("        let overlay = document.getElementById('overlay');\n")
    html_file.write("        let overlayImg = document.getElementById('overlay-img');\n")
    html_file.write("        let overlayCaption = document.getElementById('overlay-caption');\n")
    html_file.write("        let lastScroll = 0;\n\n")

    html_file.write("        // Show the overlay with a big image\n")
    html_file.write("        function showOverlay(imgSrc, imgName) {\n")
    html_file.write("            // Save current scroll\n")
    html_file.write("            lastScroll = window.scrollY;\n")
    html_file.write("            // Update the browser history state\n")
    html_file.write("            history.pushState({img: imgSrc, name: imgName, scroll: lastScroll}, '', '');\n")
    html_file.write("            // Display the overlay\n")
    html_file.write("            overlayImg.src = imgSrc;\n")
    html_file.write("            overlayCaption.textContent = imgName;\n")
    html_file.write("            overlay.style.display = 'block';\n")
    html_file.write("            document.body.style.overflow = 'hidden';\n")
    html_file.write("        }\n\n")

    html_file.write("        // Hide the overlay (triggered by the Back button)\n")
    html_file.write("        function hideOverlay() {\n")
    html_file.write("            // Use the browser's back functionality\n")
    html_file.write("            history.back();\n")
    html_file.write("        }\n\n")

    html_file.write("        // Called whenever the user clicks Back/Forward or we programmatically call history.back()\n")
    html_file.write("        window.onpopstate = function(event) {\n")
    html_file.write("            // If there's no state or no image in the state, we assume we want to show the gallery\n")
    html_file.write("            if (!event.state || !event.state.img) {\n")
    html_file.write("                // Hide the overlay\n")
    html_file.write("                overlay.style.display = 'none';\n")
    html_file.write("                document.body.style.overflow = '';\n")
    html_file.write("                // Restore scroll if available\n")
    html_file.write("                if (event.state && event.state.scroll) {\n")
    html_file.write("                    window.scrollTo(0, event.state.scroll);\n")
    html_file.write("                }\n")
    html_file.write("            } else {\n")
    html_file.write("                // The state has an image, show the overlay again\n")
    html_file.write("                overlayImg.src = event.state.img;\n")
    html_file.write("                overlayCaption.textContent = event.state.name;\n")
    html_file.write("                overlay.style.display = 'block';\n")
    html_file.write("                document.body.style.overflow = 'hidden';\n")
    html_file.write("            }\n")
    html_file.write("        }\n")

    html_file.write("    </script>\n")
    html_file.write("</body>\n</html>\n")

print(f"Gallery generated: {output_html}")
