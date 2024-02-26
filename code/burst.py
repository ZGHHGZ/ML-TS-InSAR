import copy
import logging
import re
import shutil
import time
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator, List, Optional, Tuple, Union

import asf_search
import isce  # noqa: F401
import requests
from isceobj.Sensor.TOPS.Sentinel1 import Sentinel1
from lxml import etree
from shapely import geometry


log = logging.getLogger(__name__)


URL = 'https://sentinel1-burst.asf.alaska.edu'


@dataclass
class BurstParams:
    """Parameters necessary to request a burst from the API."""

    granule: str
    swath: str
    polarization: str
    burst_number: int


class BurstMetadata:
    """Metadata for a burst."""

    def __init__(self, metadata: etree.Element, burst_params: BurstParams):
        self.safe_name = burst_params.granule
        self.swath = burst_params.swath
        self.polarization = burst_params.polarization
        self.burst_number = burst_params.burst_number
        self.manifest = metadata[0]
        self.manifest_name = 'manifest.safe'
        metadata = metadata[1]

        names = [file.attrib['source_filename'] for file in metadata]
        lengths = [len(name.split('-')) for name in names]
        swaths = [name.split('-')[length - 8] for name, length in zip(names, lengths)]
        products = [x.tag for x in metadata]
        swaths_and_products = list(zip(swaths, products))

        files = {'product': 'annotation', 'calibration': 'calibration', 'noise': 'noise'}
        for name in files:
            elem = metadata[swaths_and_products.index((self.swath.lower(), name))]
            content = copy.deepcopy(elem.find('content'))
            content.tag = 'product'
            setattr(self, files[name], content)
            setattr(self, f'{files[name]}_name', elem.attrib['source_filename'])

        file_paths = [elements.attrib['href'] for elements in self.manifest.findall('.//fileLocation')]
        pattern = f'^./measurement/s1.*{self.swath.lower()}.*{self.polarization.lower()}.*.tiff$'
        self.measurement_name = [Path(path).name for path in file_paths if re.search(pattern, path)][0]

        self.orbit_direction = self.manifest.findtext('.//{*}pass').lower()


def create_burst_request_url(params: BurstParams, content_type: str) -> str:
    """Create a URL to request a burst from the API.

    Args:
        params: The burst search parameters.
        content_type: The content type of the burst to request.

    Returns:
        A URL to request a burst from the API.
    """
    filetypes = {'metadata': 'xml', 'geotiff': 'tiff'}
    extension = filetypes[content_type]
    url = f'{URL}/{params.granule}/{params.swath}/{params.polarization}/{params.burst_number}.{extension}'
    return url


def wait_for_extractor(response: requests.Response, sleep_time: int = 15) -> bool:
    """Wait for the burst extractor to finish processing a burst.

    Args:
        response: The response from the burst extractor.
        sleep_time: The number of seconds to wait between checking the status of the burst.

    Returns:
        True if the burst was successfully downloaded, False otherwise.
    """
    if response.status_code == 202:
        time.sleep(15)
        return False

    response.raise_for_status()
    return True


def download_from_extractor(asf_session: requests.Session, burst_params: BurstParams, content_type: str) -> bytes:
    """Download burst data from the extractor.

    Args:
        asf_session: A requests session with an ASF URS cookie.
        burst_params: The burst search parameters.
        content_type: The type of content to download (metadata or geotiff).

    Returns:
        The downloaded content.
    """
    burst_request = {
        'url': create_burst_request_url(burst_params, content_type=content_type),
        'cookies': {'asf-urs': asf_session.cookies['asf-urs']},
    }

    for i in range(1, 11):
        log.info(f'Download attempt #{i} for {burst_request["url"]}')
        response = asf_session.get(**burst_request)
        downloaded = wait_for_extractor(response)
        if downloaded:
            break

    if not downloaded:
        raise RuntimeError('Download failed too many times')

    return response.content


def download_metadata(
    asf_session: requests.Session, burst_params: BurstParams, out_file: Union[Path, str] = None
) -> Union[etree._Element, str]:
    """Download burst metadata.

    Args:
        asf_session: A requests session with an ASF URS cookie.
        burst_params: The burst search parameters.
        out_file: The path to save the metadata to (if desired).

    Returns:
        The metadata as an lxml.etree._Element object or the path to the saved metadata file.
    """
    content = download_from_extractor(asf_session, burst_params, 'metadata')
    metadata = etree.fromstring(content)

    if not out_file:
        return metadata

    with open(out_file, 'wb') as f:
        f.write(content)

    return str(out_file)


def download_burst(asf_session: requests.Session, burst_params: BurstParams, out_file: Union[Path, str] = None) -> Path:
    """Download a burst geotiff.

    Args:
        asf_session: A requests session with an ASF URS cookie.
        burst_params: The burst search parameters.
        out_file: The path to save the geotiff to (if desired).

    Returns:
        The path to the saved geotiff file.
    """
    content = download_from_extractor(asf_session, burst_params, 'geotiff')

    if not out_file:
        out_file = (
            f'{burst_params.granule}_{burst_params.swath}_{burst_params.polarization}_{burst_params.burst_number}.tiff'
        ).lower()

    with open(out_file, 'wb') as f:
        f.write(content)

    return Path(out_file)


def spoof_safe(burst: BurstMetadata, burst_tiff_path: Path, base_path: Path = Path('.')) -> Path:
    """Spoof a Sentinel-1 SAFE file for a burst.

    The created SAFE file will be saved to the base_path directory. The SAFE will have the following structure:
    SLC.SAFE/
    ├── manifest.safe
    ├── measurement/
    │   └── burst.tif
    └── annotation/
        ├── annotation.xml
        └── calibration/
            ├── calibration.xml
            └── noise.xml

    Args:
        burst: The burst metadata.
        burst_tiff_path: The path to the burst geotiff.
        base_path: The path to save the SAFE file to.

    Returns:
        The path to the saved SAFE file.
    """
    safe_path = base_path / f'{burst.safe_name}.SAFE'
    annotation_path = safe_path / 'annotation'
    calibration_path = safe_path / 'annotation' / 'calibration'
    measurement_path = safe_path / 'measurement'
    paths = [annotation_path, calibration_path, measurement_path]
    for path in paths:
        path.mkdir(parents=True, exist_ok=True)

    et_args = {'encoding': 'UTF-8', 'xml_declaration': True}

    etree.ElementTree(burst.annotation).write(annotation_path / burst.annotation_name, **et_args)
    etree.ElementTree(burst.calibration).write(calibration_path / burst.calibration_name, **et_args)
    etree.ElementTree(burst.noise).write(calibration_path / burst.noise_name, **et_args)
    etree.ElementTree(burst.manifest).write(safe_path / 'manifest.safe', **et_args)

    shutil.move(str(burst_tiff_path), str(measurement_path / burst.measurement_name))

    return safe_path


def get_isce2_burst_bbox(params: BurstParams, base_dir: Optional[Path] = None) -> geometry.Polygon:
    """Get the bounding box of a Sentinel-1 burst using ISCE2.
    Using ISCE2 directly ensures that the bounding box is the same as the one used by ISCE2 for processing.

    args:
        params: The burst parameters.
        base_dir: The directory containing the SAFE file.
            If base_dir is not set, it will default to the current working directory.

    returns:
        The bounding box of the burst as a shapely.geometry.Polygon object.
    """
    if base_dir is None:
        base_dir = Path.cwd()

    s1_obj = Sentinel1()
    s1_obj.configure()
    s1_obj.safe = [str(base_dir / f'{params.granule}.SAFE')]
    s1_obj.swathNumber = int(params.swath[-1])
    s1_obj.parse()
    snwe = s1_obj.product.bursts[params.burst_number].getBbox()

    # convert from south, north, west, east -> minx, miny, maxx, maxy
    bbox = geometry.box(snwe[2], snwe[0], snwe[3], snwe[1])
    return bbox


def get_region_of_interest(
    ref_bbox: geometry.Polygon, sec_bbox: geometry.Polygon, is_ascending: bool = True
) -> Tuple[float]:
    """Get the region of interest for two bursts that will lead to single burst ISCE2 processing.

    For a descending orbit, the roi is in the lower left corner of the two bursts, and for an ascending orbit the roi is
    in the upper right corner.

    Args:
        ref_bbox: The reference burst's bounding box.
        sec_bbox: The secondary burst's bounding box.
        is_ascending: Whether the orbit is ascending or descending.

    Returns:
        The region of interest as a tuple of (minx, miny, maxx, maxy).
    """
    intersection = ref_bbox.intersection(sec_bbox)
    bounds = intersection.bounds

    x, y = (0, 1) if is_ascending else (2, 1)
    roi = geometry.Point(bounds[x], bounds[y]).buffer(0.02)
    return roi.bounds


def get_asf_session() -> requests.Session:
    """Get a requests session with an ASF URS cookie.

    requests will automatically use the netrc file:
    https://requests.readthedocs.io/en/latest/user/authentication/#netrc-authentication

    Returns:
        A requests session with an ASF URS cookie.
    """
    session = requests.Session()
    payload = {
        'response_type': 'code',
        'client_id': 'BO_n7nTIlMljdvU6kRRB3g',
        'redirect_uri': 'https://auth.asf.alaska.edu/login',
    }
    response = session.get('https://urs.earthdata.nasa.gov/oauth/authorize', params=payload)
    response.raise_for_status()
    return session


def download_bursts(param_list: Iterator[BurstParams]) -> List[BurstMetadata]:
    """Download bursts in parallel and creates SAFE files.

    For each burst:
        1. Download metadata
        2. Download geotiff
        3. Create BurstMetadata object
        4. Create directory structure
        5. Write metadata
        6. Move geotiff to correct directory

    Args:
        param_list: An iterator of burst search parameters.

    Returns:
        A list of BurstMetadata objects.
    """
    with get_asf_session() as asf_session:
        with ThreadPoolExecutor(max_workers=10) as executor:
            xml_futures = [executor.submit(download_metadata, asf_session, params) for params in param_list]
            tiff_futures = [executor.submit(download_burst, asf_session, params) for params in param_list]
            metadata_xmls = [future.result() for future in xml_futures]
            burst_paths = [future.result() for future in tiff_futures]

    bursts = []
    for params, metadata_xml, burst_path in zip(param_list, metadata_xmls, burst_paths):
        burst = BurstMetadata(metadata_xml, params)
        spoof_safe(burst, burst_path)
        bursts.append(burst)
    log.info('SAFEs created!')

    return burst




def get_burst_params(scene_name,swath_number,burst_number,) -> BurstParams:
    return BurstParams(
        granule=scene_name,
        swath=swath_number,
        polarization='VV',
        burst_number=burst_number,
    )
