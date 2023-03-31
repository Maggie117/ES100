# ES100

## Description
This project uses EL images to make solar cell power predictions.

## Usage
The EL images, extracted maximum power point (MPP) values, and the MATLAB scripts should be downloaded. The image processing and power prediction algorithms can be found under the Code directory. The EL images can be found in the Training Images or Validation directory. The following files contain the MPP values for the training set and validation set, respectively: "training_mpp.xlsx" and "measured_mpp_val.xlsx." Once everything is in the working directory in MATLAB, run the final image processing script followed by the final power prediction script. The first figure represents the lines that the algorithm was able to detect. The second figure shows the approximation error of the predictive algorithm, and the third figure shows the corresponding residual plot. The predicted MPP values can be viewed and extracted by opening the predicted_pw array in the Workspace in MATLAB.
