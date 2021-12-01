getDimensions(w, h, channels, slices, frames);
run("Stack to Hyperstack...", "order=xyczt(default) channels=1 slices=25 frames=420 display=Grayscale");
output = getDirectory("Choose Destination Directory to save results");
setBatchMode(true);
                title1=getTitle();    
                //save(output + title1);
                getDimensions(w, h, channels, slices, frames);
                for (j=1; j<=slices; j++) {                                 
                	Stack.setSlice(j);
                	run("Reduce Dimensionality...", "frames keep");                                
                	run("Align slices in stack...", "method=5 windowsizex=400 windowsizey=400 x0=50 y0=50 swindow=0 subpixel=true itpmethod=0 ref.slice=1 show=true");
                	saveAs("Tiff", output+"Slice"+toString(j)+"_"+substring(title1,0,lengthOf(title1))+".tif");
                	//run("Z Project...", "projection=[Average Intensity]");
                	//saveAs("Tiff", output+"AVGSlice"+toString(j)+"_"+substring(title1,0,10)+"_"+substring(title1,20,lengthOf(title1)));
               		selectWindow(title1);
                	close("\\Others");
                	call("java.lang.System.gc");
                }  


setBatchMode(false);