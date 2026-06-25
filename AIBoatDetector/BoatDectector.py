# This script uses two modes to reduce YOLO computation:
# 1. SEARCH MODE checks the ROI once every few seconds with YOLOv8n.
# 2. TRACK MODE rewinds a few seconds, checks the ROI at a chosen stride with YOLOv8n,
#    and records when a tracked boat crosses the vertical line.

import os
import time
import pytesseract
import cv2
import numpy as np
import scipy.io
from ultralytics import YOLO


pytesseract.pytesseract.tesseract_cmd = r"C:\Program Files\Tesseract-OCR\tesseract.exe"

start_time = time.perf_counter()


def line_intersections_with_rect(point1, point2, rect):
    # Extend the line through two points until it intersects the ROI rectangle.
    x_min, y_min, x_max, y_max = rect
    x_a, y_a = point1
    x_b, y_b = point2
    dx = x_b - x_a
    dy = y_b - y_a
    intersections = []

    # Check intersections with the left and right edges of the rectangle.
    if dx != 0:
        for x_edge in (x_min, x_max):
            t = (x_edge - x_a) / dx
            y_edge = y_a + t * dy
            if y_min <= y_edge <= y_max:
                intersections.append((int(round(x_edge)), int(round(y_edge))))

    # Check intersections with the top and bottom edges of the rectangle.
    if dy != 0:
        for y_edge in (y_min, y_max):
            t = (y_edge - y_a) / dy
            x_edge = x_a + t * dx
            if x_min <= x_edge <= x_max:
                intersections.append((int(round(x_edge)), int(round(y_edge))))

    # Corner hits can be found twice, so remove duplicate intersection points.
    unique_intersections = []
    for point in intersections:
        if point not in unique_intersections:
            unique_intersections.append(point)

    if len(unique_intersections) >= 2:
        return unique_intersections[0], unique_intersections[1]

    return point1, point2


# Load the nano model for both sparse search and active tracking.
search_model = YOLO("yolov8n.pt")
track_model = YOLO("yolov8n.pt")

# Open the video file.
video_path = r"D:\Boatvideos\October\06\cms_dock_south-2022-10-06-164803Z.mp4"
cap = cv2.VideoCapture(video_path)

# Pressure sensor location in NC State Plane feet.
pressure_sensor_east_ft = 2344242.3
pressure_sensor_north_ft = 143745.445

# Load the MATLAB calibration file that contains matched image and UTM coordinates.
mat_path = r"C:\Users\colin\OneDrive - UNC-Wilmington\WakeProjData\Summer2022\ICW_survey_20220922\mats\data.mat"
mat_data = scipy.io.loadmat(mat_path, squeeze_me=True, struct_as_record=False)

# bothCoord contains one structure per GCP with image pixel and UTM coordinate fields.
both_coord = mat_data["bothCoord"]

# Store matching calibration pairs in separate lists:
# pixel_points[i] corresponds to utm_points[i].
pixel_points = []
utm_points = []

# Read GCP 1 through GCP 22 from bothCoord.
for gcp_number in range(1, 23):
    # Open the MATLAB structure named gcp1, gcp2, ..., gcp22.
    gcp = getattr(both_coord, f"gcp{gcp_number}")

    # Ximg/Yimg are the image pixel coordinates clicked in MATLAB.
    x_img = float(gcp.Ximg)
    y_img = float(gcp.Yimg)

    # East/North are the matching UTM coordinates for that same point.
    east = float(gcp.East)
    north = float(gcp.North)

    # Add the matched pair to the calibration lists.
    pixel_points.append([x_img, y_img])
    utm_points.append([east, north])

# Convert the calibration point lists into NumPy arrays for OpenCV.
pixel_points = np.array(pixel_points, dtype=np.float32)
utm_points = np.array(utm_points, dtype=np.float32)

# Build a homography that projects full-frame image pixels into UTM coordinates.
# RANSAC helps reduce the impact of any GCP point that is slightly off.
pixel_to_utm_matrix, homography_inlier_mask = cv2.findHomography(
    pixel_points,
    utm_points,
    cv2.RANSAC,
)

# Stop early if OpenCV cannot estimate the pixel-to-UTM projection.
if pixel_to_utm_matrix is None:
    raise RuntimeError("Could not compute pixel-to-UTM homography from GCP points.")


def project_pixel_to_utm(pixel_point):
    # OpenCV expects perspective-transform input as an array of points.
    pixel_array = np.array([[pixel_point]], dtype=np.float32)

    # Project the full-frame pixel point into UTM East/North coordinates.
    utm_array = cv2.perspectiveTransform(pixel_array, pixel_to_utm_matrix)

    # Return plain Python floats so the values store cleanly in crossing_events.
    east = float(utm_array[0][0][0])
    north = float(utm_array[0][0][1])
    return east, north


# Save one full-frame JPG for each detected crossing.
image_output_folder = r"C:\Users\colin\OneDrive - UNC-Wilmington\Documents\Python test\Boatimages"
os.makedirs(image_output_folder, exist_ok=True)

# Crop this full-frame region to read the camera timestamp.
timestamp_x1, timestamp_y1, timestamp_x2, timestamp_y2 = 24, 11, 276, 59

# The camera records at 5 frames per second.
video_fps = 5

# Search mode runs YOLO once every search_seconds seconds.
search_seconds = 2
search_interval_frames = video_fps * search_seconds

# When search mode finds a boat, track mode rewinds rewind_seconds seconds.
rewind_seconds = 3
rewind_frames = video_fps * rewind_seconds

# Track mode runs the ROI line-crossing logic every track_stride_frames frames.
track_stride_frames = 5

# Track mode switches back to search mode after this many seconds with no boats detected.
no_boat_timeout_seconds = 5
no_boat_timeout_frames = video_fps * no_boat_timeout_seconds

# YOLO's COCO class ID for boat is 8.
boat_class_id = 8

# ROI coordinates in the original full video frame.
x1, y1, x2, y2 = 315, 340, 750, 675

# These settings are for a simple custom tracker based on boat center positions.
active_tracks = {}
next_track_id = 1
max_match_distance = 250
max_missing_frames = 5
direction_penalty = 100
duplicate_x_center_threshold = 2.5
minimum_track_box_area = 160

# This stores each detected crossing event and its OCR camera timestamp.
crossing_events = []
crossing_count = 0

# The script starts by sparsely searching the ROI for any boat.
mode = "search"

# These values are set when track mode begins.
track_start_frame = None
last_boat_frame = None


while cap.isOpened():
    # Read the next frame from the video.
    ret, frame = cap.read()
    if not ret:
        break

    # OpenCV's current position is after the read, so subtract 1 for this frame's index.
    frame_index = int(cap.get(cv2.CAP_PROP_POS_FRAMES)) - 1
    seconds = frame_index / video_fps

    if mode == "search":
        # In search mode, skip frames until the next scheduled search check.
        if frame_index % search_interval_frames != 0:
            continue

        # Crop the ROI so search mode only looks inside the selected window.
        roi = frame[y1:y2, x1:x2]

        # Run YOLO on the ROI only to see whether any boat is present there.
        search_results = search_model.predict(
            roi,
            classes=[boat_class_id],
            verbose=False,
        )

        # Check whether YOLO found at least one boat in the ROI.
        boat_detected = any(len(result.boxes) > 0 for result in search_results)

        # If no boat was found, stay in search mode and wait for the next interval.
        if not boat_detected:
            continue

        # If a boat was found, rewind so track mode can inspect the lead-up.
        track_start_frame = max(0, frame_index - rewind_frames)

        # Search mode just saw a boat at frame_index, so use that as the
        # latest known boat frame while track mode rewinds and catches up.
        last_boat_frame = frame_index
        active_tracks.clear()
        next_track_id = 1

        # Move the video reader back to the start of the tracking window.
        cap.set(cv2.CAP_PROP_POS_FRAMES, track_start_frame)

        # Switch into active tracking mode.
        mode = "track"
        continue

    if mode == "track":
        # If no boat has been seen for the timeout period, return to search mode.
        if frame_index - last_boat_frame >= no_boat_timeout_frames:
            print(
                f"leaving track mode - frame {frame_index}, seconds {seconds:.2f}, "
                f"last_boat_frame {last_boat_frame}"
            )
            mode = "search"
            track_start_frame = None
            last_boat_frame = None
            active_tracks.clear()
            next_track_id = 1
            continue

        # In track mode, only run YOLO on stride frames to reduce computation.
        if (frame_index - track_start_frame) % track_stride_frames != 0:
            continue

        print(f"track yolo - frame {frame_index}, seconds {seconds:.2f}")

        # Crop the ROI so line-crossing work only happens inside the selected window.
        roi = frame[y1:y2, x1:x2]

        # Place the vertical crossing line halfway across the ROI from the left.
        roi_height, roi_width = roi.shape[:2]
        line_x = roi_width / 2

        # Run YOLO detection on the ROI. Custom center matching handles tracking.
        track_results = track_model.predict(
            roi,
            classes=[boat_class_id],
            verbose=False,
        )

        total_boxes = sum(len(result.boxes) for result in track_results)
        print(f"track detections - frame {frame_index}, boxes {total_boxes}")

        # Build a list of current boat centers from YOLO detections.
        current_detections = []
        for result in track_results:
            for box in result.boxes:
                # Get the boat's bounding box inside the ROI.
                x1b, y1b, x2b, y2b = box.xyxy[0]

                # Ignore tiny track-mode detections that are likely waves or noise.
                box_width = float(x2b - x1b)
                box_height = float(y2b - y1b)
                box_area = box_width * box_height
                if box_area < minimum_track_box_area:
                    continue

                # Store the center of this detected boat for custom matching.
                x_center = float((x1b + x2b) / 2)
                y_center = float((y1b + y2b) / 2)
                current_detections.append({
                    "x": x_center,
                    "y": y_center,
                    "box": (
                        int(x1b),
                        int(y1b),
                        int(x2b),
                        int(y2b),
                    ),
                })

                print(
                    f"track box - frame {frame_index}, "
                    f"x_center {x_center:.1f}, line_x {line_x:.1f}"
                )

        # Remove near-identical detections so one boat is less likely to become two tracks.
        filtered_detections = []
        for detection in current_detections:
            duplicate_detection = any(
                abs(detection["x"] - kept_detection["x"]) <= duplicate_x_center_threshold
                for kept_detection in filtered_detections
            )

            if not duplicate_detection:
                filtered_detections.append(detection)

        current_detections = filtered_detections

        unmatched_detection_indices = set(range(len(current_detections)))
        unmatched_track_ids = set(active_tracks.keys())
        matches = []

        # Match each current detection to the closest predicted active track.
        for detection_index, detection in enumerate(current_detections):
            best_track_id = None
            best_score = None

            for track_id in unmatched_track_ids:
                track = active_tracks[track_id]
                # Predict where the track should be based on its last velocity.
                predicted_x = track["x"] + track["velocity_x"]
                predicted_y = track["y"] + track["velocity_y"]
                dx = detection["x"] - predicted_x
                dy = detection["y"] - predicted_y
                match_distance = (dx ** 2 + dy ** 2) ** 0.5
                score = match_distance

                # Penalize matches that reverse the track's current x-direction.
                movement_x = detection["x"] - track["x"]
                if track["velocity_x"] > 0 and movement_x < 0:
                    score += direction_penalty
                elif track["velocity_x"] < 0 and movement_x > 0:
                    score += direction_penalty

                if best_score is None or score < best_score:
                    best_score = score
                    best_track_id = track_id

            if best_track_id is not None and best_score <= max_match_distance:
                matches.append((best_track_id, detection_index))
                unmatched_track_ids.remove(best_track_id)
                unmatched_detection_indices.remove(detection_index)

        # Check line crossings for detections matched to existing custom tracks.
        for track_id, detection_index in matches:
            track = active_tracks[track_id]
            detection = current_detections[detection_index]
            prev_x = track["x"]
            current_x = detection["x"]
            crossed_lr = prev_x < line_x and current_x >= line_x
            crossed_rl = prev_x > line_x and current_x <= line_x

            print(
                f"crossing check - frame {frame_index}, track_id {track_id}, "
                f"prev_x {prev_x:.1f}, current_x {current_x:.1f}, "
                f"velocity_x {track['velocity_x']:.1f}, "
                f"crossed_lr {crossed_lr}, crossed_rl {crossed_rl}"
            )

            if crossed_lr or crossed_rl:
                crossing_count += 1
                direction = "left_to_right" if crossed_lr else "right_to_left"
                image_path = os.path.join(image_output_folder, f"crossing_{crossing_count}.jpg")

                # Draw debugging overlays on a full-frame copy before saving the JPG.
                annotated_frame = frame.copy()

                # Convert previous/current boxes from ROI coordinates to full-frame coordinates.
                prev_box = track["box"]
                current_box = detection["box"]
                prev_box_full = (
                    x1 + prev_box[0],
                    y1 + prev_box[1],
                    x1 + prev_box[2],
                    y1 + prev_box[3],
                )
                current_box_full = (
                    x1 + current_box[0],
                    y1 + current_box[1],
                    x1 + current_box[2],
                    y1 + current_box[3],
                )
                # Average the previous and current YOLO box sizes for this crossing event.
                prev_box_width = prev_box_full[2] - prev_box_full[0]
                prev_box_height = prev_box_full[3] - prev_box_full[1]
                current_box_width = current_box_full[2] - current_box_full[0]
                current_box_height = current_box_full[3] - current_box_full[1]
                average_box_width_pixels = (prev_box_width + current_box_width) / 2
                average_box_height_pixels = (prev_box_height + current_box_height) / 2

                # Get the bottom-left and bottom-right pixel coordinates for the previous box.
                prev_bottom_left = (prev_box_full[0], prev_box_full[3])
                prev_bottom_right = (prev_box_full[2], prev_box_full[3])

                # Get the bottom-left and bottom-right pixel coordinates for the current box.
                current_bottom_left = (current_box_full[0], current_box_full[3])
                current_bottom_right = (current_box_full[2], current_box_full[3])

                # Project the previous box bottom edge into UTM coordinates.
                prev_left_east, prev_left_north = project_pixel_to_utm(prev_bottom_left)
                prev_right_east, prev_right_north = project_pixel_to_utm(prev_bottom_right)

                # Project the current box bottom edge into UTM coordinates.
                current_left_east, current_left_north = project_pixel_to_utm(current_bottom_left)
                current_right_east, current_right_north = project_pixel_to_utm(current_bottom_right)

                # Use Pythagorean theorem to estimate previous-box bottom-edge length in feet.
                previous_box_boat_length_ft = (
                    (prev_right_east - prev_left_east) ** 2
                    + (prev_right_north - prev_left_north) ** 2
                ) ** 0.5

                # Use Pythagorean theorem to estimate current-box bottom-edge length in feet.
                current_box_boat_length_ft = (
                    (current_right_east - current_left_east) ** 2
                    + (current_right_north - current_left_north) ** 2
                ) ** 0.5

                # Average the two bottom-edge length estimates for the crossing event.
                average_boat_length_ft = (
                    previous_box_boat_length_ft + current_box_boat_length_ft
                ) / 2

                # Use bottom-center points because they are closest to the boat/water contact point.
                prev_bottom_center = (
                    int((prev_box_full[0] + prev_box_full[2]) / 2),
                    prev_box_full[3],
                )
                current_bottom_center = (
                    int((current_box_full[0] + current_box_full[2]) / 2),
                    current_box_full[3],
                )

                # Convert the ROI-relative vertical crossing line to full-frame x coordinates.
                full_line_x = x1 + line_x

                # Find where the boat path line intersects the vertical crossing line.
                # This uses the two bottom-center points from the previous and current boxes.
                path_dx = current_bottom_center[0] - prev_bottom_center[0]
                path_dy = current_bottom_center[1] - prev_bottom_center[1]

                # Most crossings have x movement; if not, fall back to the current point's y.
                if path_dx != 0:
                    intersection_t = (full_line_x - prev_bottom_center[0]) / path_dx
                    crossing_intersection_y = prev_bottom_center[1] + intersection_t * path_dy
                else:
                    crossing_intersection_y = current_bottom_center[1]

                # Store the crossing-line intersection as a full-frame pixel point.
                crossing_intersection_pixel = (
                    int(round(full_line_x)),
                    int(round(crossing_intersection_y)),
                )

                # Project the crossing-line intersection into UTM coordinates.
                crossing_intersection_east, crossing_intersection_north = project_pixel_to_utm(
                    crossing_intersection_pixel
                )

                # Calculate distance from the crossing point to the pressure sensor in feet.
                distance_to_pressure_sensor_ft = (
                    (crossing_intersection_east - pressure_sensor_east_ft) ** 2
                    + (crossing_intersection_north - pressure_sensor_north_ft) ** 2
                ) ** 0.5

                # Project both bottom-center points from image pixels into UTM coordinates.
                prev_east, prev_north = project_pixel_to_utm(prev_bottom_center)
                current_east, current_north = project_pixel_to_utm(current_bottom_center)

                # Calculate the real-world movement between the two projected boat points.
                crossing_step_distance_ft = (
                    (current_east - prev_east) ** 2
                    + (current_north - prev_north) ** 2
                ) ** 0.5

                # Convert the frame gap between previous/current detections into seconds.
                crossing_step_frame_delta = frame_index - track["last_seen_frame"]
                crossing_step_seconds = crossing_step_frame_delta / video_fps

                # Estimate boat speed over the crossing step in feet per second.
                if crossing_step_seconds > 0:
                    crossing_step_speed_ft_per_s = crossing_step_distance_ft / crossing_step_seconds
                else:
                    crossing_step_speed_ft_per_s = 0

                # Draw the full ROI box on the saved crossing image.
                cv2.rectangle(
                    annotated_frame,
                    (x1, y1),
                    (x2, y2),
                    (0, 0, 0),
                    3,
                )

                # Draw the vertical crossing x-line inside the ROI.
                cv2.line(
                    annotated_frame,
                    (int(round(full_line_x)), y1),
                    (int(round(full_line_x)), y2),
                    (0, 0, 0),
                    3,
                )

                cv2.rectangle(
                    annotated_frame,
                    (prev_box_full[0], prev_box_full[1]),
                    (prev_box_full[2], prev_box_full[3]),
                    (0, 0, 255),
                    3,
                )
                cv2.rectangle(
                    annotated_frame,
                    (current_box_full[0], current_box_full[1]),
                    (current_box_full[2], current_box_full[3]),
                    (0, 255, 0),
                    3,
                )
                # Draw only the boat path segment between the previous and current x-center points.
                cv2.line(
                    annotated_frame,
                    prev_bottom_center,
                    current_bottom_center,
                    (0, 0, 0),
                    3,
                )

                # Mark the previous and current x-center endpoints in blue.
                cv2.circle(
                    annotated_frame,
                    prev_bottom_center,
                    7,
                    (255, 0, 0),
                    -1,
                )
                cv2.circle(
                    annotated_frame,
                    current_bottom_center,
                    7,
                    (255, 0, 0),
                    -1,
                )

                # Mark the crossing-line intersection point in yellow.
                cv2.circle(
                    annotated_frame,
                    crossing_intersection_pixel,
                    7,
                    (0, 255, 255),
                    -1,
                )

                # OCR the timestamp from the original frame so annotations do not interfere.
                timestamp_crop = frame[timestamp_y1:timestamp_y2, timestamp_x1:timestamp_x2]
                timestamp_crop = cv2.cvtColor(timestamp_crop, cv2.COLOR_BGR2GRAY)
                timestamp_crop = cv2.resize(timestamp_crop, None, fx=2, fy=2)
                _, timestamp_crop = cv2.threshold(
                    timestamp_crop,
                    0,
                    255,
                    cv2.THRESH_BINARY + cv2.THRESH_OTSU,
                )
                camera_time = pytesseract.image_to_string(
                    timestamp_crop,
                    config="--psm 7",
                ).strip()

                cv2.imwrite(image_path, annotated_frame)
                crossing_events.append({
                    "frame_index": frame_index,
                    "direction": direction,
                    "camera_time": camera_time,
                    "previous_bottom_center_pixel": prev_bottom_center,
                    "current_bottom_center_pixel": current_bottom_center,
                    "average_box_width_pixels": average_box_width_pixels,
                    "average_box_height_pixels": average_box_height_pixels,
                    "previous_box_boat_length_ft": previous_box_boat_length_ft,
                    "current_box_boat_length_ft": current_box_boat_length_ft,
                    "average_boat_length_ft": average_boat_length_ft,
                    "crossing_intersection_east": crossing_intersection_east,
                    "crossing_intersection_north": crossing_intersection_north,
                    "distance_to_pressure_sensor_ft": distance_to_pressure_sensor_ft,
                    "crossing_step_distance_ft": crossing_step_distance_ft,
                    "crossing_step_speed_ft_per_s": crossing_step_speed_ft_per_s,
                })

            active_tracks[track_id] = {
                "x": detection["x"],
                "y": detection["y"],
                "box": detection["box"],
                "velocity_x": detection["x"] - track["x"],
                "velocity_y": detection["y"] - track["y"],
                "last_seen_frame": frame_index,
                "missing_count": 0,
            }

        # Start new custom tracks for detections that did not match an existing track.
        for detection_index in unmatched_detection_indices:
            detection = current_detections[detection_index]
            active_tracks[next_track_id] = {
                "x": detection["x"],
                "y": detection["y"],
                "box": detection["box"],
                "velocity_x": 0,
                "velocity_y": 0,
                "last_seen_frame": frame_index,
                "missing_count": 0,
            }
            print(f"new custom track - frame {frame_index}, track_id {next_track_id}")
            next_track_id += 1

        # Mark unmatched existing tracks as missing, then remove stale tracks.
        for track_id in list(unmatched_track_ids):
            active_tracks[track_id]["missing_count"] += 1
            print(
                f"missing custom track - frame {frame_index}, track_id {track_id}, "
                f"missing_count {active_tracks[track_id]['missing_count']}"
            )

            if active_tracks[track_id]["missing_count"] > max_missing_frames:
                print(f"removing custom track - frame {frame_index}, track_id {track_id}")
                del active_tracks[track_id]

        # If a boat was seen, reset the no-boat timeout.
        if current_detections:
            last_boat_frame = frame_index


# Release the video after processing is complete.
cap.release()
cv2.destroyAllWindows()

# Print the final number of detected crossings.
print(crossing_count)
print(crossing_events)

elapsed_time = time.perf_counter() - start_time
print(f"Elapsed time: {elapsed_time:.2f} seconds")


