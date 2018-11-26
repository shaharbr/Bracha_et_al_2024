// close all image windows
while (nImages>0) { 
	selectImage(nImages); 
	close(); 
	} 


input_dic = "path\\";

output_dic = input_dic + "splitchannel\\"; 
File.makeDirectory(output_dic); 

filelist = getFileList(input_dic);

print(input_dic);

for (i = 0; i < filelist.length; i++) {
	while (nImages>0) { 
	selectImage(nImages); 
	close(); 
	} 
     	if (endsWith(filelist[i], ".dv")) {
     		print(filelist[i]);
     		run("Bio-Formats Importer", "open="+input_dic+filelist[i]+" autoscale color_mode=Grayscale rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");
     		title = getTitle();
     		print(title);
     		run("Split Channels");
     		title = getTitle();
     		print(title);
     		saveAs("Tiff", output_dic+title);
     		close();
     		title = getTitle();
     		print(title);
     		saveAs("Tiff", output_dic+title);
     		close();
     		title = getTitle();
     		print(title);
     		saveAs("Tiff", output_dic+title);
     		close();
     		title = getTitle();
     		print(title);
     		saveAs("Tiff", output_dic+title);
     		close();
     		}
     		// close all image windows
      		while (nImages>0) { 
          		selectImage(nImages); 
          		close(); 
          		} 
			}


input_dic = output_dic;
output_dic = output_dic + "deconvolved\\"; 
File.makeDirectory(output_dic); 

filelist = getFileList(input_dic);

print(input_dic);


for (i = 0; i < filelist.length; i++) {
	while (nImages>0) { 
	selectImage(nImages); 
	close(); 
	} 
     	if (endsWith(filelist[i], ".tif")) {
     		print(filelist[i]);
     		open(filelist[i]);
     		title = getTitle();
     		print(title);
     		if (startsWith(filelist[i], "C1")) {
     			run("Diffraction PSF 3D", "index=1.470 numerical=1.36 wavelength=519 longitudinal=0 image=64.03 slice=400 width,=1024 height,=1024 depth,=9 normalization=[Sum of pixel values = 1] title=PSF");
     			run("Iterative Deconvolve 3D", "image="+title+" point=PSF output=deconvolved normalize show log perform detect wiener=0.000 low=0 z_direction=0 maximum=10 terminate=0.010");
	     		saveAs("Tiff", output_dic + filelist[i] + "_deconvolved.tif");
     		}
     		else if (startsWith(filelist[i], "C2")) {
     			run("Diffraction PSF 3D", "index=1.470 numerical=1.36 wavelength=617 longitudinal=0 image=64.03 slice=400 width,=1024 height,=1024 depth,=9 normalization=[Sum of pixel values = 1] title=PSF");
     			run("Iterative Deconvolve 3D", "image="+title+" point=PSF output=deconvolved normalize show log perform detect wiener=0.000 low=0 z_direction=0 maximum=10 terminate=0.010");
	     		saveAs("Tiff", output_dic + filelist[i] + "_deconvolved.tif");
     		}
     		else if (startsWith(filelist[i], "C3")) {
     			run("Diffraction PSF 3D", "index=1.470 numerical=1.36 wavelength=461 longitudinal=0 image=64.03 slice=400 width,=1024 height,=1024 depth,=9 normalization=[Sum of pixel values = 1] title=PSF");
     			run("Iterative Deconvolve 3D", "image="+title+" point=PSF output=deconvolved normalize show log perform detect wiener=0.000 low=0 z_direction=0 maximum=10 terminate=0.010");
	     		saveAs("Tiff", output_dic + filelist[i] + "_deconvolved.tif");
     		}
     		// close all image windows
      		while (nImages>0) { 
          		selectImage(nImages); 
          		close(); 
          		} 
			}}




// close all image windows
while (nImages>0) { 
	selectImage(nImages); 
	close(); 
	} 
