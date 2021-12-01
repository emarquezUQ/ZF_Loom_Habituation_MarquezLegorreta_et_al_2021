

output = getDirectory("Choose Destination Directory to save results");
dir = getDirectory("Choose a Directory ");
setBatchMode(true);
find(dir)

function find(dir) { 
	list = getFileList(dir); 
	for (i=0; i<list.length; i++) {           
		if (endsWith(list[i], ")/")) 
		find(""+dir+list[i]);			
		else if (endsWith(list[i], ".tif")) {
			open(""+dir+list[i]);
			title1=getTitle();  
			run("Properties...", "channels=1 slices=241 frames=1 unit=micron pixel_width=0.3226 pixel_height=0.3226 voxel_depth=1 global");
			title2=replace(title1,".tif",".nrrd");
			run("Nrrd ... ", "nrrd="+output+toString(title2));
			selectWindow(title1);			
			close();
		}

	}

}

//run("Properties...", "channels=1 slices=50 frames=1 unit=micron pixel_width=1.5 pixel_height=1.5 voxel_depth=5 global");


//run("Nrrd ... ", "nrrd=D:/temp/NRRD/hs_as_1022_fish1.nrrd");

