import os
import shutil
import subprocess
import zipfile

from core.logger import get_logger
from core.error_capture import capture_error

BASE = os.path.dirname(os.path.abspath(__file__))

LOG_PATH = os.path.join(BASE, "logs", "arcgis_run.log")
ERROR_PATH = os.path.join(BASE, "logs", "error_traceback.log")

UPLOAD_DIR = os.path.join(BASE, "codex_sync", "upload_here")
ZIP_PATH = os.path.join(BASE, "codex_sync", "package.zip")

SCRIPT_TO_RUN = os.path.join(BASE, "SCRIPTS", "zone_crime_fire.py")

logger = get_logger(LOG_PATH)


def run_script():
    logger.info("==== Starting arc GIS RUN =====")

    try:
        result = subprocess.run(
            ["python", SCRIPT_TO_RUN],
            capture_output=True,
            text=True,
            check=False,
        )
        logger.info(result.stdout)

        if result.stderr:
            logger.error(result.stderr)

    except Exception:
        capture_error(ERROR_PATH)
        logger.error("Fatal error captured")


def export_to_codex():
    os.makedirs(UPLOAD_DIR, exist_ok=True)

    for file_path in [LOG_PATH, ERROR_PATH]:
        if os.path.exists(file_path):
            shutil.copy(file_path, UPLOAD_DIR)


def zip_for_mobile():
    os.makedirs(os.path.dirname(ZIP_PATH), exist_ok=True)
    with zipfile.ZipFile(ZIP_PATH, "w") as zipf:
        for folder, _, files in os.walk(UPLOAD_DIR):
            for file_name in files:
                full_path = os.path.join(folder, file_name)
                zipf.write(full_path, arcname=file_name)


def main():
    run_script()
    export_to_codex()
    zip_for_mobile()
    print("AMDS Complete -> Ready for mobile upload")


if __name__ == "__main__":
    main()