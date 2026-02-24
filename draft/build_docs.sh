rm -rf ./public
cd draft/main
rm -rf _freeze
quarto render 
cd ../supp

# for svg_file in figures/supp/*.svg; 
# do
#     echo "Converting $svg_file to pdf"
#     no_ext=${svg_file%.svg}
#     rsvg-convert  --format pdf --output "$no_ext".pdf $svg_file
#     rsvg-convert  --format png --keep-aspect-ratio --dpi-x 396 --dpi-y 396 --output "$no_ext".png $svg_file
# done

# merge the header only file with the rawmd.not_md into the supplementary.qmd
quarto render 
cd ../
# python build_supplementary_pdf.py

# cp ./draft_main_manuscript.docx ../public/index.docx

exit