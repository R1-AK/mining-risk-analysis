import os
from flask import Flask, render_template, request, jsonify, send_from_directory, url_for
from werkzeug.utils import secure_filename
import geopandas as gpd
from pathlib import Path
import json
import sys
import zipfile
import tempfile
import shutil

# This section is for the pyngrok integration
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

# Use relative paths for robustness in Colab/other environments
app.config['UPLOAD_FOLDER'] = Path('./uploads')
app.config['UPLOAD_FOLDER'].mkdir(exist_ok=True)
Path('./outputs').mkdir(exist_ok=True) # Ensure outputs folder exists

def handle_zip_upload(file_storage, upload_folder):
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

@app.route('/')
def home():
    return render_template('home.html')

@app.route('/spatial_analysis_page')
def spatial_analysis_page():
    return render_template('spatial_analysis.html')

@app.route('/flowpath_analysis_page')
def flowpath_analysis_page():
    return render_template('flowpath_analysis.html')

@app.route('/download/<path:filename>')
def download(filename):
    """Serves a file for download from the 'outputs' directory."""
    return send_from_directory('outputs', filename, as_attachment=True)

@app.route('/get_shapefile_columns', methods=['POST'])
def get_shapefile_columns():
    if 'shapefile' not in request.files: return jsonify({"error": "No shapefile part"}), 400
    file = request.files['shapefile']
    if file.filename == '': return jsonify({"error": "No selected file"}), 400
    if not file.filename.lower().endswith('.zip'): return jsonify({"error": "Must be a .zip archive"}), 400
    temp_dir = None
    try:
        shp_filepath, temp_dir = handle_zip_upload(file, app.config['UPLOAD_FOLDER'])
        gdf = gpd.read_file(shp_filepath)
        columns = [col for col in gdf.columns if col.lower() != 'geometry']
        return jsonify({"columns": columns})
    except Exception as e:
        return jsonify({"error": f"Failed to read shapefile: {str(e)}"}), 500
    finally:
        if temp_dir and temp_dir.exists(): shutil.rmtree(temp_dir)

@app.route('/run_spatial_analysis', methods=['POST'])
def handle_spatial_analysis():
    if not run_spatial_analysis: return jsonify({"error": "Spatial module not loaded."}), 500
    temp_dirs = []
    try:
        if 'mining_data' not in request.files: return jsonify({"error": "Mining data required."}), 400
        mining_file = request.files['mining_data']
        mining_filepath, mining_temp_dir = handle_zip_upload(mining_file, app.config['UPLOAD_FOLDER'])
        temp_dirs.append(mining_temp_dir)
        buffer_kms = [int(km.strip()) for km in request.form.get('buffer_distances').split(',')]
        risk_data_config = json.loads(request.form.get('risk_data_config'))
        for config in risk_data_config:
            file_key = config['file_input_id']
            if file_key in request.files:
                risk_file = request.files[file_key]
                if risk_file.filename.lower().endswith('.zip'):
                    risk_filepath, risk_temp_dir = handle_zip_upload(risk_file, app.config['UPLOAD_FOLDER'])
                    temp_dirs.append(risk_temp_dir)
                else:
                    risk_filename = secure_filename(risk_file.filename)
                    risk_filepath = Path(app.config['UPLOAD_FOLDER']) / risk_filename
                    risk_file.save(risk_filepath)
                config['path'] = str(risk_filepath)
            else:
                return jsonify({"error": f"Missing file for {config['name']}"}), 400
        
        result_filename = run_spatial_analysis(str(mining_filepath), buffer_kms, risk_data_config)
        
        if "No results" in result_filename:
            return jsonify({"status": "success", "message": result_filename, "results_file": None})
        
        download_url = url_for('download', filename=result_filename, _external=True)
        return jsonify({"status": "success", "message": "Analysis complete!", "results_file": download_url})
    except Exception as e:
        import traceback
        traceback.print_exc()
        return jsonify({"error": f"An error occurred: {str(e)}"}), 500
    finally:
        for temp_dir in temp_dirs:
            if temp_dir and temp_dir.exists(): shutil.rmtree(temp_dir)
        for item in (app.config['UPLOAD_FOLDER']).glob('*'):
            if item.is_file(): item.unlink()

@app.route('/run_flowpath_analysis', methods=['POST'])
def handle_flowpath_analysis():
    if not run_flowpath_analysis: return jsonify({"error": "Flowpath module not loaded."}), 500
    gee_key_filepath, mining_filepath_flow = None, None
    mining_temp_dir = None
    try:
        if 'mining_data' not in request.files or 'gee_key_file' not in request.files:
            return jsonify({"error": "Mining data and GEE key required."}), 400
        mining_file = request.files['mining_data']
        mining_filepath_flow, mining_temp_dir = handle_zip_upload(mining_file, app.config['UPLOAD_FOLDER'])
        gee_key_file = request.files['gee_key_file']
        gee_key_filename = secure_filename(gee_key_file.filename)
        gee_key_filepath = Path(app.config['UPLOAD_FOLDER']) / gee_key_filename
        gee_key_file.save(gee_key_filepath)
        gee_credentials = {"SERVICE_ACCOUNT": request.form.get('gee_service_account'), "SERVICE_ACCOUNT_KEY_FILE": str(gee_key_filepath), "PROJECT_ID": request.form.get('gee_project_id')}
        parameters = {"SCALE": int(request.form.get('scale')), "RADIUS_M": int(request.form.get('radius_m'))}
        risk_indicators = request.form.getlist('indicators')
        
        result_filename = run_flowpath_analysis(str(mining_filepath_flow), gee_credentials, parameters, risk_indicators)
        
        if "No results" in result_filename:
            return jsonify({"status": "success", "message": result_filename, "results_file": None})
            
        download_url = url_for('download', filename=result_filename, _external=True)
        return jsonify({"status": "success", "message": "Flowpath analysis completed!", "results_file": download_url})
    except Exception as e:
        return jsonify({"error": f"An error occurred: {str(e)}"}), 500
    finally:
        if gee_key_filepath and gee_key_filepath.exists(): gee_key_filepath.unlink()
        if mining_temp_dir and mining_temp_dir.exists(): shutil.rmtree(mining_temp_dir)

if __name__ == '__main__':
    authtoken = os.environ.get("NGROK_AUTHTOKEN")
    if authtoken:
        ngrok.set_auth_token(authtoken)
    else:
        print("!!! Ngrok authtoken not set. Get a token from https://dashboard.ngrok.com/get-started/your-authtoken and add it as a secret in Colab named NGROK_AUTHTOKEN !!!")
    
    public_url = ngrok.connect(5000)
    print(f" * ngrok tunnel available at: {public_url}")
    app.run()
