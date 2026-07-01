# READ ME

The main file that should be used in the AIBoatDetector folder is BoatFinder.ipynb. BoatFinder was created in google colab as a python notebook. It was made to do all
of its computing virtually using googles computers. Make sure to switch the computing type to GPU for faster computing (should be free). All other files in this folder are
are not finished nor important. The code is designed to create its own virtual environment everytime it is run. The only things that are designed to be changed is below in the "Modify these" section. Since this is computing in the cloud, all imported files and exported files must be stored in your google drive. The filepaths for these files must be changed. The main filepaths include data.mat, stored in this repo folder, and the export filepath for the CSV that is created. When running the code, a pop-up should occur that will want you to sign into your 
google drive.

## What does BoatFinder.ipynb do?

BoatFinder takes a specified start and end date of interest and analyzes all available CMS south-view camera videos within those dates.  A region of interest (ROI) is created around the center portion of the frame, and YOLO is used to detect boats inside that ROI. The code records an event when a tracked boat crosses the selected vertical x-line. A GCP survey is used to create a homography matrix that converts full-frame pixel coordinates into NC State Plane feet/UTM-style Easting and Northing coordinates; tidal influence has not been adjusted for. For each crossing, the frame before the boat crosses the x-line and the frame after it crosses the x-line are captured. The bottom-center point of each YOLO box is converted into real-world coordinates, and the Pythagorean theorem is used to estimate the distance traveled between those two detections. The elapsed frame count is divided by the video FPS to calculate elapsed time, and speed is calculated as distance traveled over elapsed time. A line is also created between the before and after boat positions. Where that line intersects the x-line, the pixel coordinate is saved, converted into real-world coordinates, and compared to the pressure sensor location with the Pythagorean theorem to estimate boat distance from the sensor. The camera time is calculated from the video filename plus the crossing frame offset; the filename timestamp is slightly behind the true displayed camera time, so a manual time correction is applied with `camera_time_offset_seconds`. In the current state of the code, the corrected `time` should be accurate to approximately +/- 1 second.

<img width="1920" height="1080" alt="Image" src="https://github.com/user-attachments/assets/da2e53d1-e7d7-40a8-8706-dee53652a2a2" /> 

This is an image of what the yolo is looking at. The box is made fairly small to avoid capturing unwanted boats.

Checkpointing: the final processing loop writes results to the CSV after each full day is processed. This means completed days are already saved in `results_csv_path` before the next day begins, so a crash or disconnect should only risk the current unfinished day. The loop prints a checkpoint message after each day showing how many rows were appended and how many videos were skipped.

## Main Variables (can be modified in "Configuration, Model, and Calibration" block):

# Modify these:

- `start_date`: first archive date to process, formatted as `YYYY-MM-DD`.
- `end_date`: final archive date to process, formatted as `YYYY-MM-DD`.
- `results_csv_path`: Google Drive path where the final CSV is written.
- `mat_path`: path to `data.mat`, which stores the GCP image and real-world coordinate pairs.
- `pressure_sensor_east_ft`: pressure sensor Easting coordinate in NC State Plane feet.
- `pressure_sensor_north_ft`: pressure sensor Northing coordinate in NC State Plane feet.

# Should not need modification:

- `VIDEO_FPS`: fixed video frame rate used to convert frame counts into seconds. Should consistently be 5. 
- `search_seconds`: how often search mode checks for boats.
- `rewind_seconds`: how far back track mode rewinds after search mode finds a boat.
- `track_stride_frames`: frame spacing used while tracking an active boat event.
- `no_boat_timeout_seconds`: time with no boat detections before returning to search mode.
- `duplicate_crossing_cooldown_seconds`: suppresses extra crossing events detected within this many seconds of the last accepted crossing.
- `camera_time_offset_seconds`: manual timestamp correction added after filename/frame timing; negative values move the output time earlier.
- `boat_class_id`: YOLO class ID used for boats.
- `x1, y1, x2, y2`: ROI coordinates in the original full-frame image.
- `max_match_distance`: maximum pixel distance allowed when matching detections between frames.
- `max_missing_frames`: number of missed track updates allowed before removing a custom track.
- `direction_penalty`: extra matching cost applied when a detection moves opposite the expected direction.
- `duplicate_x_center_threshold`: removes duplicate detections with nearly identical x-centers.
- `minimum_track_box_area`: minimum YOLO box area accepted during track mode.
- `maximum_valid_boat_length_ft`: skips a crossing when both before/after boat-length estimates are greater than this value.
- `slow_boat_speed_threshold_ft_per_s`: speed threshold used to identify slow-boat duplicate crossing risk.
- `slow_boat_cooldown_seconds`: cooldown time used to suppress repeated slow-boat crossings.
- `search_confidence`: YOLO confidence threshold used during search mode.
- `track_confidence`: YOLO confidence threshold used during track mode.
- `search_model`: YOLO model used for sparse search.
- `track_model`: YOLO model used for active tracking.
- `pixel_to_utm_matrix`: homography matrix that converts image pixels to real-world coordinates.

## CSV Outputs

- `video_url`: source video URL processed for this row. Every URL within the date range is pasted, even if there is no captured boat.
- `skip_reason`: reason a video was skipped; blank for successful crossing events.
- `time`: crossing time calculated from the video filename timestamp minus 4 hours, plus `frame_index / VIDEO_FPS`, plus `camera_time_offset_seconds`.
- `frame_index`: video frame index where the crossing was detected.
- `direction`: crossing direction, either `left_to_right` or `right_to_left`.
- `previous_bottom_center_pixel`: full-frame pixel coordinate of the boat's bottom-center point before crossing.
- `current_bottom_center_pixel`: full-frame pixel coordinate of the boat's bottom-center point after crossing.
- `crossing_intersection_pixel`: full-frame pixel coordinate where the boat path line intersects the x-line.
- `crossing_step_elapsed_seconds`: elapsed time between the previous and current boat detections used for the crossing.
- `average_box_width_pixels`: average YOLO box width from the before and after crossing detections.
- `average_box_height_pixels`: average YOLO box height from the before and after crossing detections.
- `previous_box_boat_length_ft`: estimated real-world length of the previous YOLO box bottom edge (UTM conversion/pythagorean theorem).
- `current_box_boat_length_ft`: estimated real-world length of the current YOLO box bottom edge (UTM conversion/pythagorean theorem).
- `average_boat_length_ft`: average of the previous and current box-based boat length estimates (UTM conversion/pythagorean theorem).
- `crossing_intersection_east`: real-world Easting coordinate of the path/x-line intersection (UTM conversion).
- `crossing_intersection_north`: real-world Northing coordinate of the path/x-line intersection (UTM conversion).
- `distance_to_pressure_sensor_ft`: distance from the path/x-line intersection to the pressure sensor, in feet (UTM conversion/pythagorean theorem).
- `crossing_step_distance_ft`: real-world distance traveled between the previous and current bottom-center detections, in feet (UTM conversion/pythagorean theorem).
- `crossing_step_speed_ft_per_s`: estimated speed between the previous and current detections, in feet per second.

  ## Issues

  This is not perfect code. The boat detection works roughly 90% of the time. When heavy traffic occurs, the boat detection is less accurate. It does not pick up false positives
  but it will tend to miss boats occasionally. For standard isolated boats, this code works very well.
