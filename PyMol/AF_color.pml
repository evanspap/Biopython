load $pdb_file
hide everything
show cartoon
select b_50, b < 50
color orange, b_50
select 50_b_70, b <70  & b >50
color yellow, 50_b_70
select 70_b_90, b <90 & b >70
color cyan, 70_b_90
select 90_b, b >90
color blue, 90_b
orient
bg_color white
zoom visible
png $output_png, width=1920, height=1080, dpi=300
quit
