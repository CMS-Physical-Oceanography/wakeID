from datetime import date, datetime, timedelta
from urllib.parse import urljoin
import re

import requests


DEFAULT_BASE_URL = "https://stage-ams.srv.axds.co/archive/mp4/uncw/cms_dock_south/"


def _parse_date(value):
    if isinstance(value, date):
        return value

    return datetime.strptime(value, "%Y-%m-%d").date()


def _date_range(start_date, end_date):
    current_date = _parse_date(start_date)
    final_date = _parse_date(end_date)

    if current_date > final_date:
        raise ValueError("start_date must be before or equal to end_date")

    while current_date <= final_date:
        yield current_date
        current_date += timedelta(days=1)


def get_video_urls(start_date, end_date, base_url=DEFAULT_BASE_URL):
    """
    Return a list of CMS dock MP4 URLs between start_date and end_date.

    start_date and end_date can be strings in "YYYY-MM-DD" format or Python
    date objects. The date range is inclusive, so both boundary dates are used.
    """
    base_url = base_url.rstrip("/") + "/"
    video_urls = []

    for current_date in _date_range(start_date, end_date):
        day_url = urljoin(
            base_url,
            f"{current_date:%Y}/{current_date:%m}/{current_date:%d}/",
        )

        response = requests.get(day_url, timeout=30)
        response.raise_for_status()

        filenames = list(dict.fromkeys(re.findall(r"cms.*?\.mp4", response.text)))

        for filename in filenames:
            video_urls.append(urljoin(day_url, filename))

    return video_urls
