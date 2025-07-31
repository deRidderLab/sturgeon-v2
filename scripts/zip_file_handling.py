import os
import zipfile

def unpack_zip_file(zip_path):
    # Get the absolute path of the directory where this script is located
    script_dir = os.path.dirname(os.path.abspath(__file__))

    include_dir = os.path.join(script_dir, '..', 'include')

    # Get the base name of the zip file without extension
    zip_base = os.path.splitext(os.path.basename(zip_path))[0]

    # Construct the expected extraction path
    extract_dir = os.path.join(include_dir, zip_base)

    # Check if the directory already exists
    if os.path.exists(extract_dir) and os.path.isdir(extract_dir):
        return extract_dir

    # If not, unzip the contents
    with zipfile.ZipFile(zip_path, 'r') as zip_ref:
        zip_ref.extractall(extract_dir)

    return extract_dir

