- README -

For further anaylsis, the stiffness map is in the variable 'stiffness_map' and the r-square in the variable 'rsquare_map'

1.) Open a spectroscopy file with XEI software
2.) Choose 'Young's modulus' panel
3.) If you want to export only point data files, look for a list of point data under a reference image. 
4.) Currently, only the 1st data point is chosen, you have to move a scroll bar of the point data panel to the buttom and press SHIFT while checking the final data point. (see figure 1.0)
5.) Move your mouse cursor to the center of the point data panel, then right click and choose export -> Point data 
6.) Specify your data directory (Right click + New Folder). This is the directory that you have to fill in the variable 'DataLocation' inside the MATLAB script file 'data_reader.m'.
7.) Follow the steps as explained in the "data_reader.m"

## If you want to show the reference image you need to follow the following steps
1.) Follow the steps 1-7 above
2.) Move your mouse cursor to the center of the poin data panel, then right click and choose export -> image (You have to put the image file (*.png) in the same directory as the point data)
3.) right click again and choose export -> text to the same directory as the image file and the point data file
4.) Follow the steps as explained in the "data_reader.m" 




