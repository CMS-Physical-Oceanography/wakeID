import cv2
import scipy.io


video_path = r"D:\Boatvideos\October\02\cms_dock_south-2022-10-02-163337Z.mp4"
mat_path = r"C:\Users\colin\OneDrive - UNC-Wilmington\WakeProjData\Summer2022\ICW_survey_20220922\mats\data.mat"
preview_width = 1200

mat_data = scipy.io.loadmat(mat_path, squeeze_me=True, struct_as_record=False)
both_coord = mat_data["bothCoord"]

cap = cv2.VideoCapture(video_path)
ret, frame = cap.read()
cap.release()

if not ret:
    raise RuntimeError("Could not read the first frame from the video.")

for gcp_number in range(1, 23):
    gcp = getattr(both_coord, f"gcp{gcp_number}")
    x_img = int(round(float(gcp.Ximg)))
    y_img = int(round(float(gcp.Yimg)))
    east = float(gcp.East)
    north = float(gcp.North)

    cv2.circle(frame, (x_img, y_img), 6, (0, 0, 255), -1)
    cv2.putText(
        frame,
        f"gcp{gcp_number}",
        (x_img + 8, y_img - 8),
        cv2.FONT_HERSHEY_SIMPLEX,
        0.6,
        (0, 0, 255),
        2,
    )

    print(f"gcp{gcp_number}: Ximg={x_img}, Yimg={y_img}, East={east:.3f}, North={north:.3f}")

frame_height, frame_width = frame.shape[:2]
scale = min(1, preview_width / frame_width)
preview_height = int(frame_height * scale)
preview = cv2.resize(frame, (int(frame_width * scale), preview_height))

cv2.imshow("GCP Pixel Coordinates", preview)
cv2.waitKey(0)
cv2.destroyAllWindows()
