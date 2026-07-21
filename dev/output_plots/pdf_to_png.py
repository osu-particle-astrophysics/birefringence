import argparse
import os
from pdf2image import convert_from_path

def convert_pdfs_to_pngs(input_files, output_dir, dpi=300):
    """Convert multiple PDF files to PNG images."""
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    for pdf_file in input_files:
        if not os.path.isfile(pdf_file):
            print(f"File not found: {pdf_file}")
            continue

        try:
            print(f"Converting: {pdf_file}")
            # Convert PDF to images
            images = convert_from_path(pdf_file, dpi=dpi)
            
            # Save each page as a PNG file
            pdf_basename = os.path.splitext(os.path.basename(pdf_file))[0]
            for i, image in enumerate(images):
                output_path = os.path.join(output_dir, f"{pdf_basename}.png")
                image.save(output_path, "PNG")
                print(f"Saved: {output_path}")
        except Exception as e:
            print(f"Error processing {pdf_file}: {e}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert PDF files to PNG images.")
    parser.add_argument("--input_files", nargs='+', required=True, help="List of PDF files to convert.")
    parser.add_argument("--output_dir", default="converted_pngs", help="Directory to save the PNG files.")
    parser.add_argument("--dpi", type=int, default=300, help="Resolution of the output PNGs (default: 300).")
    
    args = parser.parse_args()
    
    convert_pdfs_to_pngs(args.input_files, args.output_dir, args.dpi)

