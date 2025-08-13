import io
import sys
import json
import grass.script as gs
from grass.gunittest.case import TestCase
from grass.gunittest.main import test
from grass.pygrass.utils import get_lib_path
from grass.pygrass.vector.geometry import Point
from unittest.mock import patch, MagicMock
import base64
import tempfile
from pathlib import Path

path = get_lib_path(modname="i.overlap", libname="metalib")
if path is None:
    gs.fatal("Not able to find the metalib library directory.")
sys.path.append(path)

import metalib as ml


class TestIOverlap(TestCase):
    input = "test_input"
    output = "test_output"
    data_dir = "data"
    photos = sorted(Path(data_dir).glob("*.jpg"))
    create_transformer = ml.create_transformer()
    exif = {
        "GPSInfo": {
            "GPSVersionID": b"\x02\x03\x00\x00",
            "GPSLatitudeRef": "N",
            "GPSLatitude": (35.0, 37.0, 10.85),
            "GPSLongitudeRef": "W",
            "GPSLongitude": (82.0, 24.0, 1.12),
            "GPSAltitudeRef": b"\x00",
            "GPSAltitude": 1750.0,
            "GPSMapDatum": "WGS84",
        },
        "ResolutionUnit": 2,
        "ExifOffset": 216,
        "Make": "Canon",
        "Model": "Canon EOS 5DS R",
        "YResolution": 72.0,
        "Orientation": 1,
        "DateTime": "2024:10:03 17:41:16",
        "YCbCrPositioning": 2,
        "Copyright": "",
        "XResolution": 72.0,
        "Artist": "",
        "ExifVersion": b"0231",
        "ComponentsConfiguration": b"\x01\x02\x03\x00",
        "ShutterSpeedValue": 11.0,
        "DateTimeOriginal": "2024:10:03 17:41:16",
        "DateTimeDigitized": "1999:12:31 21:18:46",
        "ApertureValue": 5.0,
        "ExposureBiasValue": 0.0,
        "MeteringMode": 5,
        "Flash": 16,
        "FocalLength": 50.0,
        "ColorSpace": 1,
        "ExifImageWidth": 8688,
        "ExifInteroperabilityOffset": 10152,
        "FocalPlaneXResolution": 5991.724137931034,
        "FocalPlaneYResolution": 6002.072538860104,
        "OffsetTime": "-05:00",
        "OffsetTimeOriginal": "-05:00",
        "OffsetTimeDigitized": "-05:00",
        "SubsecTime": "17",
        "SubsecTimeOriginal": "254",
        "SubsecTimeDigitized": "17",
        "ExifImageHeight": 5792,
        "FocalPlaneResolutionUnit": 2,
        "ExposureTime": 0.0005,
        "FNumber": 5.6,
        "ExposureProgram": 1,
        "CustomRendered": 0,
        "ISOSpeedRatings": 1000,
        "ExposureMode": 1,
        "FlashPixVersion": b"0100",
        "SensitivityType": 2,
        "WhiteBalance": 1,
        "RecommendedExposureIndex": 1000,
        "CameraOwnerName": "",
        "BodySerialNumber": "384055000156",
        "LensSpecification": (50.0, 50.0, 0.0, 0.0),
        "LensModel": "EF50mm f/1.4 USM",
        "LensSerialNumber": "0000000000",
        "SceneCaptureType": 0,
    }

    @classmethod
    def setUpClass(cls):
        cls.use_temp_region()
        cls.runModule("g.region", n=10, s=0, e=10, w=0, res=1)

    @classmethod
    def tearDownClass(cls):
        # cls.runModule("g.remove", flags="f", type="raster", name=cls.input)
        cls.del_temp_region()

    def tearDown(cls):
        """Remove output map after each test method"""
        pass
        # cls.runModule("g.remove", flags="f", type="raster", name=cls.output)

    def test_get_exif_data(self):
        """Test get_exif_data function"""
        photo = self.photos[0]
        exif_data = ml.get_exif_data(photo)

        for key, value in exif_data.items():
            print(f"{key}: {value}")
            self.assertIn(key, exif_data, f"{key} should be in EXIF data")
            # self.assertEqual(
            #     exif_data[key], value, f"{key} value mismatch in EXIF data"
            # )
        self.assertIsNotNone(exif_data, "EXIF data should not be None")
        self.assertIn("GPSInfo", exif_data, "GPSInfo should be in EXIF data")

    # def test_get_coords(self):
    #     """Test get_coords function"""
    #     lon, lat, alt = get_coords(self.exif)
    #     self.assertAlmostEqual(lon, -82.4031, places=4)
    #     self.assertAlmostEqual(lat, 35.6197, places=4)
    #     self.assertAlmostEqual(alt, 1750.0, places=1)

    # def test_get_coords_invalid(self):
    #     """Test get_coords with invalid EXIF data"""
    #     invalid_exif = {"GPSInfo": {}}
    #     with self.assertRaises(ValueError) as cm:
    #         get_coords(invalid_exif)
    #     self.assertEqual(str(cm.exception), "Invalid GPS coordinates in EXIF data")

    # def test_output_exists(self):
    #     """Test output map exists"""
    #     self.assertModule("i.overlap", input=self.input, output=self.output)
    #     self.assertRasterExists(name=self.output, msg="Output was not created")

    # def test_get_focal_length(self):
    #     """Test focal length extraction from EXIF data"""
    #     focal_length = get_focal_length(self.exif)
    #     self.assertEqual(focal_length, 50.0, "Focal length should be 50.0 mm")

    # def test_get_focal_length_invalid(self):
    #     """Test focal length extraction with invalid EXIF data"""
    #     invalid_exif = {"FocalLength": None}
    #     focal_length = get_focal_length(invalid_exif)
    #     self.assertIsNone(focal_length, "Focal length should be None for invalid data")

    # def test_compute_sensor_size(self):
    #     """Test sensor size computation"""
    #     sensor_size = compute_sensor_size(self.exif)
    #     self.assertEqual(
    #         sensor_size, (36.0, 24.0), "Sensor size should be 36.0 x 24.0 mm"
    #     )

    # def test_compute_gsd(self):
    #     """Test GSD computation"""
    #     focal_length = get_focal_length(self.exif)
    #     sensor_size = compute_sensor_size(self.exif)
    #     gsd = self.create_transformer.compute_gsd(focal_length, sensor_size, self.input)
    #     self.assertGreater(gsd, 0, "GSD should be a positive value")

    # def test_transform(self):
    #     """Test create_transformer function"""
    #     gps = get_coords(self.exif)
    #     lon, lat, alt = gps
    #     transformer = self.create_transformer
    #     self.assertIsNotNone(transformer, "Transformer should not be None")
    #     rlon, rlat = transformer.transform(lon, lat)
    #     self.assertAlmostEqual(rlon, 0, place=4)
    #     self.assertAlmostEqual(rlat, 0, place=4)

    # def test_get_orientation(self):
    #     """Test orientation extraction from EXIF data"""
    #     orientation = get_orientation(self.exif)
    #     self.assertEqual(orientation, 1, "Orientation should be 1 (normal)")

    # def test_make_footprint(self):
    #     """Test footprint creation from EXIF data"""
    #     gps = get_coords(self.exif)
    #     lon, lat, alt = gps
    #     focal_length = get_focal_length(self.exif)
    #     sensor_size = compute_sensor_size(self.exif)
    #     footprint = make_footprint(lon, lat, alt, focal_length, sensor_size)
    #     self.assertIsNotNone(footprint, "Footprint should not be None")
    #     self.assertEqual(len(footprint), 4, "Footprint should have 4 corners")

    # def test_calculate_overlaps(self):
    #     """Test overlap calculation"""
    #     gps = get_coords(self.exif)
    #     lon, lat, alt = gps
    #     focal_length = get_focal_length(self.exif)
    #     sensor_size = compute_sensor_size(self.exif)
    #     footprint = make_footprint(lon, lat, alt, focal_length, sensor_size)
    #     overlaps = calculate_overlaps(footprint, self.input)
    #     self.assertIsNotNone(overlaps, "Overlaps should not be None")
    #     self.assertGreater(len(overlaps), 0, "There should be some overlaps calculated")

    # def test_write_vector(self):
    #     """Test vector writing"""
    #     gps = get_coords(self.exif)
    #     lon, lat, alt = gps
    #     focal_length = get_focal_length(self.exif)
    #     sensor_size = compute_sensor_size(self.exif)
    #     footprint = make_footprint(lon, lat, alt, focal_length, sensor_size)
    #     vector_name = "test_vector"
    #     write_vector(vector_name, footprint)
    #     self.assertVectorExists(name=vector_name, msg="Vector was not created")
    #     self.runModule("g.remove", flags="f", type="vector", name=vector_name)
    #     self.assertVectorExists(name=vector_name, msg="Vector was not created")


if __name__ == "__main__":
    test()
