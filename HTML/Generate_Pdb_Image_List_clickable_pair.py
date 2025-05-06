#!/usr/bin/env python3
"""
Generate_Pdb_Image_List_clickable.py

This script generates an HTML file displaying a gallery of PNG images from a specified directory.
For each "base" image (e.g., ADSS2.png), if there's a matching "*_relative_simple_heatmap.png"
(e.g., ADSS2_relative_simple_heatmap.png), then upon clicking the base image, both images will be
shown side by side in an overlay. The overlay also displays the filename below the images.

We maintain:
- Four columns at ~24% width in the main gallery.
- Images keep their original aspect ratio.
- On click, an overlay opens with up to two images side by side (each ~49% max width), also
  preserving their aspect ratios.
- The scroll position is restored upon returning to the gallery.
- Browser back button and a custom Back button in the overlay are supported.

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

# Check if directory exists
if not os.path.isdir(image_dir):
    print(f"Error: Directory '{image_dir}' does not exist.")
    sys.exit(1)

# Gather all PNG files
all_pngs = sorted([f for f in os.listdir(image_dir) if f.endswith('.png')])
all_pngs_set = set(all_pngs)

# We'll treat each file that does NOT end with '_relative_simple_heatmap.png' as a 'main' image.
# For each main image, check if there's a partner named <basename>_relative_simple_heatmap.png.

pairs = []  # Each element: (main_image, heatmap_image_or_none, base_name)

for f in all_pngs:
    if f.endswith('_relative_simple_heatmap.png'):
        # Skip heatmap images here because
        # we only drive the gallery from the main image.
        continue
    base_name = os.path.splitext(f)[0]
    # Potential partner
    heatmap_candidate = base_name + '_relative_simple_heatmap.png'
    if heatmap_candidate in all_pngs_set:
        pairs.append((f, heatmap_candidate, base_name))
    else:
        pairs.append((f, None, base_name))

# Generate HTML
with open(output_html, 'w') as html_file:
    html_file.write("<!DOCTYPE html>\n")
    html_file.write("<html lang='en'>\n<head>\n")
    html_file.write("    <meta charset='UTF-8'>\n")
    html_file.write("    <meta name='viewport' content='width=device-width, initial-scale=1.0'>\n")
    html_file.write("    <title>PDB Image Gallery</title>\n")

    # CSS Styles
    html_file.write("    <style>\n")
    html_file.write("        body { font-family: Arial, sans-serif; margin: 0; padding: 0; text-align: center; }\n")
    html_file.write("        .header { position: fixed; top: 0; width: 100%; background: white; padding: 10px 0; box-shadow: 0px 4px 2px -2px gray; z-index: 1000; }\n")
    html_file.write("        .content { margin-top: 60px; }\n")

    # Table with four columns at ~24% each
    html_file.write("        table { width: 100%; border-collapse: collapse; table-layout: fixed; }\n")
    html_file.write("        td { width: 24%; padding: 10px; border: 1px solid #ccc; text-align: center; vertical-align: middle; }\n")

    # Images in gallery: preserve aspect ratio, up to the cell's width
    html_file.write("        .gallery-img {\n")
    html_file.write("            display: block;\n")
    html_file.write("            margin: auto;\n")
    html_file.write("            max-width: 100%;\n")
    html_file.write("            height: auto;\n")
    html_file.write("            cursor: pointer;\n")
    html_file.write("        }\n")

    # Overlay area
    html_file.write("        #overlay {\n")
    html_file.write("            position: fixed;\n")
    html_file.write("            top: 0; left: 0; right: 0; bottom: 0;\n")
    html_file.write("            background: white;\n")
    html_file.write("            display: none;\n")
    html_file.write("            z-index: 2000;\n")
    html_file.write("            overflow: auto;\n")
    html_file.write("            padding: 20px;\n")
    html_file.write("        }\n")
    html_file.write("        #overlay-content {\n")
    html_file.write("            max-width: 1600px;\n")
    html_file.write("            margin: 40px auto;\n")
    html_file.write("        }\n")
    html_file.write("        .overlay-img-container {\n")
    html_file.write("            display: flex;\n")
    html_file.write("            flex-wrap: wrap;\n")
    html_file.write("            gap: 20px;\n")
    html_file.write("            justify-content: center;\n")
    html_file.write("        }\n")

    # Each image is up to 49% wide, maintaining aspect ratio
    html_file.write("        .overlay-img {\n")
    html_file.write("            display: block;\n")
    html_file.write("            margin: auto;\n")
    html_file.write("            max-width: 49%;\n")
    html_file.write("            height: auto;\n")
    html_file.write("            box-sizing: border-box;\n")
    html_file.write("        }\n")

    html_file.write("        #backButton {\n")
    html_file.write("            padding: 10px 20px;\n")
    html_file.write("            font-size: 16px;\n")
    html_file.write("            margin-top: 20px;\n")
    html_file.write("        }\n")
    html_file.write("        #overlay-caption {\n")
    html_file.write("            text-align: center;\n")
    html_file.write("            font-size: 20px;\n")
    html_file.write("            margin-top: 20px;\n")
    html_file.write("        }\n")

    html_file.write("    </style>\n</head>\n")

    html_file.write("<body>\n")
    html_file.write("    <div class='header'><h1>PDB Image Gallery</h1></div>\n")

    # Gallery region
    html_file.write("    <div class='content' id='gallery'>\n")
    html_file.write("    <table>\n")

    for i, (main_img, heatmap_img, base_name) in enumerate(pairs):
        if i % 4 == 0:
            html_file.write("        <tr>\n")
        second_img_str = heatmap_img if heatmap_img else ''
        html_file.write("            <td>\n")
        html_file.write(
            f"                <img class='gallery-img' src='{image_dir}/{main_img}' alt='{main_img}' "
            f"onclick=\"showOverlay('{image_dir}/{main_img}','{image_dir}/{second_img_str}','{base_name}')\" >\n"
        )
        html_file.write(f"                <br>{base_name}\n")
        html_file.write("            </td>\n")
        if i % 4 == 3:
            html_file.write("        </tr>\n")

    if len(pairs) % 4 != 0:
        html_file.write("        </tr>\n")

    html_file.write("    </table>\n")
    html_file.write("    </div>\n")

    # The overlay for side-by-side images
    html_file.write("    <div id='overlay'>\n")
    html_file.write("        <div id='overlay-content'>\n")
    html_file.write("            <div class='overlay-img-container'>\n")
    html_file.write("                <img id='overlay-img1' class='overlay-img' src='' alt='Image'>\n")
    html_file.write("                <img id='overlay-img2' class='overlay-img' src='' alt='' style='display: none;'>\n")
    html_file.write("            </div>\n")
    html_file.write("            <p id='overlay-caption'></p>\n")
    html_file.write("            <button id='backButton' onclick='hideOverlay()'>Back</button>\n")
    html_file.write("        </div>\n")
    html_file.write("    </div>\n")

    # Script for single-page transitions
    html_file.write("    <script>\n")
    html_file.write("        let overlay = document.getElementById('overlay');\n")
    html_file.write("        let overlayImg1 = document.getElementById('overlay-img1');\n")
    html_file.write("        let overlayImg2 = document.getElementById('overlay-img2');\n")
    html_file.write("        let overlayCaption = document.getElementById('overlay-caption');\n")
    html_file.write("        let lastScroll = 0;\n\n")

    # showOverlay function
    html_file.write("        function showOverlay(imgSrc1, imgSrc2, baseName) {\n")
    html_file.write("            lastScroll = window.scrollY;\n")
    html_file.write("            // pushState so the browser Back/Forward buttons are recognized\n")
    html_file.write("            history.pushState({main: imgSrc1, heatmap: imgSrc2, name: baseName, scroll: lastScroll}, '', '');\n")
    html_file.write("            overlayImg1.src = imgSrc1;\n")
    html_file.write("            if (imgSrc2) {\n")
    html_file.write("                overlayImg2.src = imgSrc2;\n")
    html_file.write("                overlayImg2.style.display = 'block';\n")
    html_file.write("            } else {\n")
    html_file.write("                overlayImg2.src = '';\n")
    html_file.write("                overlayImg2.style.display = 'none';\n")
    html_file.write("            }\n")
    html_file.write("            overlayCaption.textContent = baseName;\n")
    html_file.write("            overlay.style.display = 'block';\n")
    html_file.write("            document.body.style.overflow = 'hidden';\n")
    html_file.write("        }\n")

    # hideOverlay function
    html_file.write("        function hideOverlay() {\n")
    html_file.write("            history.back();\n")
    html_file.write("        }\n")

    # onpopstate handling
    html_file.write("        window.onpopstate = function(event) {\n")
    html_file.write("            if (!event.state || !event.state.main) {\n")
    html_file.write("                // We are showing the gallery\n")
    html_file.write("                overlay.style.display = 'none';\n")
    html_file.write("                document.body.style.overflow = '';\n")
    html_file.write("                // Restore scroll if in state\n")
    html_file.write("                if (event.state && event.state.scroll) {\n")
    html_file.write("                    window.scrollTo(0, event.state.scroll);\n")
    html_file.write("                }\n")
    html_file.write("            } else {\n")
    html_file.write("                // Show overlay\n")
    html_file.write("                overlay.style.display = 'block';\n")
    html_file.write("                document.body.style.overflow = 'hidden';\n")
    html_file.write("                overlayImg1.src = event.state.main;\n")
    html_file.write("                if (event.state.heatmap) {\n")
    html_file.write("                    overlayImg2.src = event.state.heatmap;\n")
    html_file.write("                    overlayImg2.style.display = 'block';\n")
    html_file.write("                } else {\n")
    html_file.write("                    overlayImg2.src = '';\n")
    html_file.write("                    overlayImg2.style.display = 'none';\n")
    html_file.write("                }\n")
    html_file.write("                overlayCaption.textContent = event.state.name;\n")
    html_file.write("            }\n")
    html_file.write("        }\n")

    html_file.write("    </script>\n")

    html_file.write("</body>\n</html>\n")

print(f"Gallery generated: {output_html}")
