# Tests YOLO boat detection on a chosen video frame and ROI, then displays the detections.
import cv2
from ultralytics import YOLO


model = YOLO("yolov8n.pt")

video_path = r"D:\Boatvideos\October\02\cms_dock_south-2022-10-02-160338Z.mp4"
x1, y1, x2, y2 =  315, 340, 750, 675
boat_class_id = 8
preview_width = 1200
frame_index_to_test = 190

cap = cv2.VideoCapture(video_path)
cap.set(cv2.CAP_PROP_POS_FRAMES, frame_index_to_test)
ret, frame = cap.read()
cap.release()

if not ret:
    raise RuntimeError(f"Could not read frame {frame_index_to_test} from the video.")

roi = frame[y1:y2, x1:x2]
results = model.predict(
    roi,
    classes=[boat_class_id],
    verbose=False,
)

roi_with_boxes = roi.copy()
full_frame_with_roi = frame.copy()
cv2.rectangle(full_frame_with_roi, (x1, y1), (x2, y2), (0, 255, 0), 3)

for result in results:
    for box in result.boxes:
        x1b, y1b, x2b, y2b = box.xyxy[0]
        confidence = float(box.conf[0])

        x1b = int(x1b)
        y1b = int(y1b)
        x2b = int(x2b)
        y2b = int(y2b)
        box_area = (x2b - x1b) * (y2b - y1b)

        print(f"box area: {box_area} pixels, confidence: {confidence:.2f}")

        cv2.rectangle(roi_with_boxes, (x1b, y1b), (x2b, y2b), (0, 0, 255), 2)
        cv2.putText(
            roi_with_boxes,
            f"area {box_area}px",
            (x1b, max(20, y1b - 10)),
            cv2.FONT_HERSHEY_SIMPLEX,
            0.7,
            (0, 0, 255),
            2,
        )

scale = min(1, preview_width / full_frame_with_roi.shape[1])
preview_height = int(full_frame_with_roi.shape[0] * scale)
full_frame_preview = cv2.resize(full_frame_with_roi, (preview_width, preview_height))

print(f"Detected {sum(len(result.boxes) for result in results)} boat boxes on frame {frame_index_to_test}.")

cv2.imshow("Full Frame With ROI", full_frame_preview)
cv2.imshow("ROI YOLO Boat Detections", roi_with_boxes)
cv2.waitKey(0)
cv2.destroyAllWindows()
