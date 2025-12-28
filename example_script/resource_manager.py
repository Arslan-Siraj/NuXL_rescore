import io
import zipfile
import requests
from pathlib import Path

def ensure_resources(nuxl_rescore_dir: Path, zip_url: str ) -> None:
    """
    Download and extract nuxl_rescore_resource.zip if directory is empty.
    Pure Python version (no Streamlit).
    """

    # If not empty → nothing to do
    if any(nuxl_rescore_dir.iterdir()):
        return

    print("Resources missing for rescoring. Downloading NuXL rescore resources...")

    # Stream download to avoid memory issues
    with requests.get(zip_url, stream=True) as r:
        r.raise_for_status()
        total_size = int(r.headers.get("Content-Length", 0))
        downloaded = 0
        zip_buffer = io.BytesIO()

        for chunk in r.iter_content(chunk_size=8192):
            if chunk:
                zip_buffer.write(chunk)
                downloaded += len(chunk)

                # Optional lightweight progress display
                if total_size > 0:
                    percent = downloaded * 100 // total_size
                    print(f"\rDownloading... {percent}%", end="")

    print("\nExtracting files...")

    # Extract ZIP
    zip_buffer.seek(0)
    with zipfile.ZipFile(zip_buffer) as z:
        z.extractall(nuxl_rescore_dir)

    print("Resources downloaded and extracted successfully.")