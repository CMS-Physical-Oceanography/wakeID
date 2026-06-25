# Shows and saves the first video frame with a chosen ROI box, plus the cropped ROI,
# so the x1, y1, x2, y2 coordinates can be checked before using them in the detector.

# Put the ROI coordinates you want to test here.
x1, y1, x2, y2 = 270, 359, 1400, 675

### This code is meant to be run before the main BoatDectector.py, to verify the ROI coordinates.
#------------------------------------------------------------------------------------------------
import cv2


video_path = r"D:\Boatvideos\October\01\cms_dock_south-2022-10-01-110405Z.mp4"

preview_width = 1200
frame_output_path = "first_frame_with_roi.jpg"
roi_output_path = "first_frame_roi_crop.jpg"

cap = cv2.VideoCapture(video_path)
ret, frame = cap.read()
cap.release()

if not ret:
    raise RuntimeError("Could not read the first frame from the video.")

frame_height, frame_width = frame.shape[:2]

x1 = max(0, min(x1, frame_width))
x2 = max(0, min(x2, frame_width))
y1 = max(0, min(y1, frame_height))
y2 = max(0, min(y2, frame_height))

if x2 <= x1 or y2 <= y1:
    raise ValueError("Invalid ROI coordinates. Make sure x2 > x1 and y2 > y1.")

roi = frame[y1:y2, x1:x2]
frame_with_roi = frame.copy()
cv2.rectangle(frame_with_roi, (x1, y1), (x2, y2), (0, 255, 0), 3)

cv2.imwrite(frame_output_path, frame_with_roi)
cv2.imwrite(roi_output_path, roi)

print(f"x1, y1, x2, y2 = {x1}, {y1}, {x2}, {y2}")
print(f"Full frame size: width={frame_width}, height={frame_height}")
print(f"ROI size: width={x2 - x1}, height={y2 - y1}")
print(f"Saved marked frame: {frame_output_path}")
print(f"Saved ROI crop: {roi_output_path}")

scale = min(1, preview_width / frame_width)
preview_height = int(frame_height * scale)
preview = cv2.resize(frame_with_roi, (preview_width, preview_height))

cv2.imshow("First Frame With ROI", preview)
cv2.imshow("ROI Crop", roi)
cv2.waitKey(0)
cv2.destroyAllWindows()
