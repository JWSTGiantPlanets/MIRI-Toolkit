# Important

For MIRI MRS stripe correction: When spec3.cube_build.coord_system = 'ifualign': the pipeline aligns in the FOV of the IFU in the y direction. This may result in the 'horizontal' stripes still appearing diagonal. To make the stripes actually be horizontal (required for calibration) the IFU should be aligned in the x direction. To do this, replace ifu_cube.py under the directory path of form:

/opt/anaconda3/envs/jwst/lib/python3.10/site-packages/jwst/cube_build

with the script calibration_pipeline/ifu_align.py

# MIRI-Toolkit
Repository for MIRI visualisation and analysis tools

## Available Codes
!!!! Note that all codes here are experimental, use at your own risk !!!!

* `MIRIClean_FFT.ipynb` - Attempt to clean IFU-aligned MIRI frames by inspecting and filtering their 2D FFTs - this doesn't currently work as implemented, we're exploring other techniques.


## Known Challenges

* MIRI/MRS Wavelength solution currently insufficient, leading to striping in images associated with narrow emission lines.  The wavelengths must be shifted spaxel by spaxel to a uniform solution.
* MRS Saturation for Jupiter beyond 11 µm - this was expected, but workarounds are needed to access the unsaturated frames in the 4/5-group integrations.
