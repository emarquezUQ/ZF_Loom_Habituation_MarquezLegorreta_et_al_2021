

dir = getDirectory("Choose a Directory from which to load Data");
output = getDirectory("Destination Directory to save results");
setBatchMode(true);
find(dir)

function find(dir) { 
	list = getFileList(dir); 
	for (i=0; i<list.length; i++) {           
		if (endsWith(list[i], ")/")) 
			find(""+dir+list[i]);			
		else if (endsWith(list[i], "/")) {
			open(""+dir+list[i]);			
			title1=getTitle();    
			getDimensions(w, h, channels, slices, frames);
			run("Stack to Hyperstack...", "order=xyczt(default) channels=1 slices="+toString(slices/3)+" frames=3 display=Grayscale");	/// here you need to change the frames to your video								
			selectWindow(title1);
			run("Re-order Hyperstack ...", "channels=[Channels (c)] slices=[Frames (t)] frames=[Slices (z)]");
			run("Z Project...", "projection=[Average Intensity] all");
					
				//title2="AVG_"+substring(title1,0,lengthOf(title1))+".tif";
				title2=getTitle();
				run("Z Project...", "projection=[Average Intensity]");
				wait(200);
				selectWindow(title2);
				save(output + title2); 	
				if (isOpen("Exception")) {
					selectWindow("Exception");
					run("Close");
				}
				selectWindow(title1);		
				close("\\Others");
				call("java.lang.System.gc");
			 
		}
		close();
		call("java.lang.System.gc"); 
	} 
} 

setBatchMode(false);


