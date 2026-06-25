"""
Download one video from a URL.

Purpose:
    BoatFinder.ipynb will create a list of video URLs using getvideourls.py.
    This file handles the next step: take one URL from that list, download the
    video to a local folder, and return the local file path.

Typical use:
    from videodownload import download_video

    video_url = video_urls[num]
    video_path = download_video(video_url, output_folder="/content")

    # Run pythontester/process_video on video_path.
    # Delete video_path after processing so storage does not fill up.

Input:
    video_url:
        Full URL to one MP4 video.
    output_folder:
        Folder where the video should be downloaded. In Colab, "/content" is
        good temporary storage.

Output:
    The local path to the downloaded video file as a string.
"""

from pathlib import Path

import requests


def download_video(video_url, output_folder="/content"):
    """
    Download one video URL and return the local file path.
    """
    output_folder = Path(output_folder)
    output_folder.mkdir(parents=True, exist_ok=True)

    # Preserve the original filename from the URL.
    filename = video_url.rstrip("/").split("/")[-1]
    local_path = output_folder / filename

    # Stream in chunks so the whole video is never loaded into memory at once.
    with requests.get(video_url, stream=True, timeout=30) as response:
        response.raise_for_status()

        with open(local_path, "wb") as file:
            for chunk in response.iter_content(chunk_size=1024 * 1024):
                if chunk:
                    file.write(chunk)

    return str(local_path)
