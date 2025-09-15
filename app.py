import os
from flask import Flask, render_template, request, jsonify
from werkzeug.utils import secure_filename
import geopandas as gpd
from pathlib import Path
import json
import sys
import zipfile
import tempfile
import shutil

# This section is for the new pyngrok integration
from pyngrok import ngrok

# Add the current directory to the path to ensure modules can be imported
sys.path.insert(0, str(Path(__file__).parent.resolve()))

# Import the main analysis functions from your refactored scripts
try:
    from spatial_analysis import run_spatial_analysis
    from flowpath_analysis import run_flowpath_analysis
    print("Successfully imported analysis modules.")
except ImportError as e:
    print(f"Error importing analysis modules: {e}")
    run_spatial_analysis = None
    run_flowpath_analysis = None

app = Flask(__name__)

# NOTE: We have removed the flask_ngrok lines from here

# Create a dedicated, temporary folder for uploads
# Note: Using relative paths is more robust in Colab
app.config['UPLOAD_FOLDER'] = Path('./uploads')
app.config['UPLOAD_FOLDER'].mkdir(exist_ok=True)


def handle_zip_upload(file_storage, upload_folder):
    """
    Saves and unzips a shapefile archive, returning the path to the .shp file.
    """
    filename = secure_filename(file_storage.filename)
    temp_dir = Path(tempfile.mkdtemp(dir=upload_folder))
    zip_path = temp_dir / filename
    file_storage.save(zip_path)

    with zipfile.ZipFile(zip_path, 'r') as zip_ref:
        zip_ref.extractall(temp_dir)

    shp_files = list(temp_dir.glob('*.shp'))
    if not shp_files:
        shutil.rmtree(temp_dir)
        raise ValueError("The uploaded .zip file does not contain a .shp file.")

    return shp_files[0], temp_dir

# --- All your @app.route functions remain the same ---

@app.route('/')
def home():
    return render_template('home.html')

@app.route('/spatial_analysis_page')
def spatial_analysis_page():
    return render_template('spatial_analysis.html')

@app.route('/flowpath_analysis_page')
def flowpath_analysis_page():
    return render_template('flowpath_analysis.html')

@app.route('/get_shapefile_columns', methods=['POST'])
def get_shapefile_columns():
    if 'shapefile' not in request.files:
        return jsonify({"error": "No shapefile part in the request"}), 400
    file = request.files['shapefile']
    if file.filename == '':
        return jsonify({"error": "No selected file"}), 400
    if not file.filename.lower().endswith('.zip'):
        return jsonify({"error": "Please upload a .zip archive containing all shapefile components."}), 400
    temp_dir = None
    try:
        shp_filepath, temp_dir = handle_zip_upload(file, app.config['UPLOAD_FOLDER'])
        gdf = gpd.read_file(shp_filepath)
        columns = [col for col in gdf.columns if col.lower() != 'geometry']
        return jsonify({"columns": columns})
    except Exception as e:
        return jsonify({"error": f"Failed to read shapefile: {str(e)}"}), 500
    finally:
        if temp_dir and temp_dir.exists():
            shutil.rmtree(temp_dir)

@app.route('/run_spatial_analysis', methods=['POST'])
def handle_spatial_analysis():
    # This function remains the same
    if not run_spatial_analysis:
        return jsonify({"error": "Spatial analysis module not loaded."}), 500
    # ... (rest of the function code is the same)
    return jsonify({"message": "Placeholder for spatial analysis results"})


@app.route('/run_flowpath_analysis', methods=['POST'])
def handle_flowpath_analysis():
    # This function remains the same
    # ... (rest of the function code is the same)
    return jsonify({"message": "Placeholder for flowpath analysis results"})

# --- The main execution block is now updated to use pyngrok ---

if __name__ == '__main__':
    # Get the authtoken from Colab secrets
    # IMPORTANT: Add your token as a secret named NGROK_AUTHTOKEN in your Colab notebook
    authtoken = os.environ.get("NGROK_AUTHTOKEN")
    if authtoken:
        ngrok.set_auth_token(authtoken)
    else:
        print("!!! Ngrok authtoken not set. Get a token from https://dashboard.ngrok.com/get-started/your-authtoken and add it as a secret in Colab named NGROK_AUTHTOKEN !!!")

    # Open a http tunnel on the default port 5000
    public_url = ngrok.connect(5000)
    print(f" * ngrok tunnel available at: {public_url}")

    # Start the Flask app
    app.run()
