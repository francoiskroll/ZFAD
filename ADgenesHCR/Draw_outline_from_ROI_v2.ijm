roiManager("reset");

//full brain
// roiManager("open", "/Users/francoiskroll/Dropbox/ZFAD/ADgenesHCR/outline/RoiSet_fullbrainoutline_new+5.zip");
//max projection
// roiManager("open", "/Users/francoiskroll/Dropbox/ZFAD/ADgenesHCR/outline/Roi_AVG_outline_+5.roi");

//sagittal view max projection
roiManager("open", "/Users/francoiskroll/Dropbox/ZFAD/ADgenesHCR/outline/Roi_MAX_outline_+5_sagittal.roi");

getDimensions(width, height, channels, slices, frames);


n = roiManager('count');
for (i = 0; i < n; i++) {
    roiManager('select', i);
    
    if (Roi.contains(0, 0)) {
		run("Select None");
	    } else {
		Overlay.addSelection
		Overlay.addSelection("white", "2");
		Overlay.show
		run("Select None");


    }
}

	
	roiManager("reset");



