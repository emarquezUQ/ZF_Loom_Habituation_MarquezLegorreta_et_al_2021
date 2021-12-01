dir = getDirectory("Choose a Directory ");
output = getDirectory("Choose Destination Directory to save results");
setBatchMode(true);

list = getFileList(dir); 
for (i=0; i<list.length; i++) { 
	open(""+dir+list[i]);
	title1=getTitle();    
	//save(output + title1);
	for (j=1; j<=50; j++) {			
		Stack.setSlice(j);
		run("Reduce Dimensionality...", "frames keep");			
		run("Align slices in stack...", "method=5 windowsizex=450 windowsizey=450 x0=100 y0=100 swindow=0 subpixel=true itpmethod=0 ref.slice=1 show=true");
		saveAs("Tiff", output+"Slice"+toString(j)+"_"+substring(title1,0,10)+"_"+substring(title1,20,lengthOf(title1)));
		selectWindow(title1);
		close("\\Others");
		call("java.lang.System.gc");
	}  
}

setBatchMode(false);
