//Automated image analysis tools for trypanosomes
//=============================================================================================================================
//A set of macros for use in ImageJ
//  ImageJ is a free and open source cross platform piece of scientific image analysis software
//  For more information and to download ImageJ please visit: http://rsbweb.nih.gov/ij/
//  For more information about macros in ImageJ please visit: http://rsbweb.nih.gov/ij/developer/macro/macros.html
//Install these macros via "Plugins>Macros>Install" or "Ctrl+Shift+M"
//
//Copyright 2011 Richard J Wheeler (www.richardwheeler.net)
//This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as
//published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.
//This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
//of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
//You should have received a copy of the GNU General Public License along with this program; if not, see
//http://www.gnu.org/licenses.

var slice_phase=3;
var slice_dapi=1;
var slice_pi=2;
var image_scale=6.236;
var distort_scale=1;
var distort_xori=0;
var distort_yori=0;
var dapi_nucleus=1000;
var dapi_kinetoplast=1000;
var dapi_background=200;
var pi_nucleus=1000;
var pi_kinetoplast=2000;
var pi_background=200;

var analyseall=true;

//Macro - Image Setup Tool
//=============================================================================================================================
//Global variables changed
//Slice numbers
//	slice_phase
//	slice_dapi
//	slice_pi
//Image scale
//	image_scale

macro "Image Setup Action Tool -C888L88fcLfc8fL8f1cL1c88C00fL84f8Lf88bL8b18L1884Cf00L80f4Lf488L8814L1480" {
	setBatchMode(true);

	Dialog.create("Image Setup");
		Dialog.addNumber("Phase slice number:", slice_phase, 0, 5, "");
		Dialog.addNumber("DAPI slice number:", slice_dapi, 0, 5, "");
		Dialog.addNumber("PI slice number:", slice_pi, 0, 5, "");
		Dialog.addNumber("Image scale:", image_scale, 3, 5, "px/um");
	Dialog.show();
		slice_phase=Dialog.getNumber();
		slice_dapi=Dialog.getNumber();
		slice_pi=Dialog.getNumber();
		image_scale=Dialog.getNumber();

	setBatchMode(false);
}

//Macro - Measure Chromatic Aberation Tool
//=============================================================================================================================
//Global variables used
//Slice numbers
//	slice_dapi
//	slice_pi
//Global variables modified
//PI channel distort variables
//	distort_scale
//	distort_xori
//	distort_yori

var mcattolerance=200;
var mcatmaxdist=5;

macro "Measure Chromatic Aberation Action Tool -Cf88o1022Cf00o8177C00fo2322C88fo9577C000L0fffL0c0fL4d4fL8c8fLcdcfLfcff" {
	setBatchMode(true);

	//Display options
	printresults=false;
	drawplot=false;
	scalefactor=20;
	arrowwidth=8;
	Dialog.create("Measure Chromatic Aberation");
		Dialog.addCheckbox("Analyse all images currently open", analyseall);
		Dialog.addMessage("Analysis options:");
		Dialog.addNumber("Noise tolerance:", mcattolerance, 0, 5, "");
		Dialog.addNumber("Maximum pairing distance:", mcatmaxdist, 0, 5, "px");
		Dialog.addMessage("Output options:");
		Dialog.addCheckbox("Print data to log", printresults);
		Dialog.addCheckbox("Make image plot of distortion", drawplot);
		Dialog.addNumber("     Scale factor:", scalefactor, 0, 5, "px/px");
		Dialog.addNumber("     Arrow width:", arrowwidth, 0, 5, "px");
	Dialog.show();
		analyseall=Dialog.getCheckbox();
		mcattolerance=Dialog.getNumber();
		mcatmaxdist=Dialog.getNumber();
		printresults=Dialog.getCheckbox();
		drawplot=Dialog.getCheckbox();
		scalefactor=Dialog.getNumber();
		arrowwidth=Dialog.getNumber();

	if (analyseall==true) {
		images=newArray(nImages());
		for (n=0; n<nImages; n++) {
			selectImage(n+1);
			images[n]=getImageID();
		}
	} else {
		images=newArray(1);
		images[0]=getImageID();
	}

	finalxd=newArray(0);
	finalyd=newArray(0);
	finalxp=newArray(0);
	finalyp=newArray(0);
	finald=newArray(0);
	for (m=0; m<lengthOf(images); m++) {
		//Batch mode and select image
		setBatchMode(true);
		selectImage(images[m]);

		//Get points
		setSlice(slice_dapi);
		run("Find Maxima...", "noise="+mcattolerance+" output=[Point Selection]");
		getSelectionCoordinates(dapix, dapiy);
		setSlice(slice_pi);
		run("Find Maxima...", "noise="+mcattolerance+" output=[Point Selection]");
		getSelectionCoordinates(pix, piy);
		run("Select None");

		//Assign points to pairs
		dapip=newArray(lengthOf(dapix));
		dapid=newArray(lengthOf(dapix));
		for (i=0; i<lengthOf(dapix); i++) {
			mindist=getWidth()+getHeight();
			for (j=0; j<lengthOf(pix); j++) {
				xd=dapix[i]-pix[j];
				yd=dapiy[i]-piy[j];
				dist=pow(xd*xd+yd*yd, 0.5);
				if (dist<mindist) {
					mindist=dist;
					dapip[i]=j+1;
					dapid[i]=dist;
				}
			}
		}

		//Tidy pairs arrays
		p=0;
		for (i=0; i<lengthOf(dapix); i++) {
			if (dapip[i]!=0 && dapid[i]<=mcatmaxdist) {
				p++;
			}
		}
		tfinalxd=newArray(p);
		tfinalyd=newArray(p);
		tfinalxp=newArray(p);
		tfinalyp=newArray(p);
		tfinald=newArray(p);
		p=0;
		for (i=0; i<lengthOf(dapix); i++) {
			if (dapip[i]!=0 && dapid[i]<=mcatmaxdist) {
				partner=dapip[i]-1;
				tfinalxd[p]=dapix[i];
				tfinalyd[p]=dapiy[i];
				tfinalxp[p]=pix[partner];
				tfinalyp[p]=piy[partner];
				tfinald[p]=dapid[i];
				p++;
			}
		}

		//append t[final] arrays to [final] arrays
		finalxd=appendArrayToArray(finalxd, tfinalxd);
		finalyd=appendArrayToArray(finalyd, tfinalyd);
		finalxp=appendArrayToArray(finalxp, tfinalxp);
		finalyp=appendArrayToArray(finalyp, tfinalyp);
		finald=appendArrayToArray(finald, tfinald);
	}

	if (printresults==true) {
		//Print results excluding those over the max dist
		print("n dapix dapiy pix piy dist");
		n=0;
		for (i=0; i<lengthOf(finalxd); i++) {
			print(n, finalxd[i], finalyd[i], finalxp[i], finalyp[i], finald[i]);
			n++;
		}
	}

	if (drawplot==true) {
		setBatchMode(false);
		//Draws an exaggerated plot of the distortion
		newImage("Distortion", "32-bit Black", getWidth(), getHeight(), 1);
		setBatchMode(true);
		for (i=0; i<lengthOf(finalxd); i++) {
			xd=finalxd[i]-finalxp[i];
			yd=finalyd[i]-finalyp[i];
			x2=finalxd[i]+scalefactor*xd;
			y2=finalyd[i]+scalefactor*yd;
			if (xd!=0 || yd!=0) {
				drawArrow(finalxd[i], finalyd[i], x2, y2, arrowwidth, 1);
			} else {
				makeRectangle(finalxd[i]-arrowwidth/2, finalyd[i]-arrowwidth/2, arrowwidth, arrowwidth);
				run("Add...", "value=1");
			}
		}
		run("Enhance Contrast", "saturated=0.0");
		updateDisplay();
	}

	//Calculate linear regression parameters for result
	//Calculate xdif and ydif arrays
	xdif=newArray(lengthOf(finalxd));
	ydif=newArray(lengthOf(finalxd));
	for (i=0; i<lengthOf(finalxd); i++) {
		xdif[i]=finalxd[i]-finalxp[i];
		ydif[i]=finalyd[i]-finalyp[i];
	}
	dxsum=0;
	dysum=0;
	dixsum=0;
	diysum=0;
	for (i=0; i<lengthOf(finalxd); i++) {
		dxsum=dxsum+finalxd[i];
		dysum=dysum+finalyd[i];
		dixsum=dixsum+xdif[i];
		diysum=diysum+ydif[i];
	}
	dxmean=dxsum/lengthOf(finalxd);
	dymean=dysum/lengthOf(finalxd);
	dixmean=dixsum/lengthOf(finalxd);
	diymean=diysum/lengthOf(finalxd);
	xsxx=0;
	xsyy=0;
	xsxy=0;
	ysxx=0;
	ysyy=0;
	ysxy=0;
	for (i=0; i<lengthOf(finalxd); i++) {
		xsxx=xsxx+pow(finalxd[i]-dxmean, 2);
		xsyy=xsyy+pow(xdif[i]-dixmean, 2);
		xsxy=xsxy+(finalxd[i]-dxmean)*(xdif[i]-dixmean);
		ysxx=ysxx+pow(finalyd[i]-dymean, 2);
		ysyy=ysyy+pow(ydif[i]-diymean, 2);
		ysxy=ysxy+(finalyd[i]-dymean)*(ydif[i]-diymean);
	}
	xb=xsxy/xsxx;
	yb=ysxy/ysxx;
	xa=dixmean-dxmean*xb;
	ya=diymean-dymean*yb;
	if (printresults==true) {
		print("Transformation required to make PI signal overlay DAPI:");
		print("Scale: "+1+(xb+yb)/2);
		print("Origin: ("+-xa/xb+", "+-ya/yb+")");
	}

	Dialog.create("Distortion Analysis Result");
	Dialog.addMessage(""+lengthOf(images)+" images analysed.");
//****I HACKED THIS!!!!!!*****
	Dialog.addMessage("Result:\nScale: "+1+((xb+yb)/2)*1.2+"\nOrigin: ("+-xa/xb+", "+-ya/yb+")");
	Dialog.show();
//****I HACKED THIS!!!!!!*****
	distort_scale=1+((xb+yb)/2)*1.2;
	distort_xori=-xa/xb;
	distort_yori=-ya/yb;

	setBatchMode(false);
}

//Distortion test function - appendArrayToArray
//-----------------------------------------------------------------------------------------------------------------------------
//Joins two arrays together and returns the result
function appendArrayToArray(array1, array2) {
	result=newArray(lengthOf(array1)+lengthOf(array2));
	for (i=0; i<lengthOf(array1); i++) {
		result[i]=array1[i];
	}
	for (i=0; i<lengthOf(array2); i++) {
		result[i+lengthOf(array1)]=array2[i];
	}
	return result;
}

//Distortion test function - drawArrow
//-----------------------------------------------------------------------------------------------------------------------------
//Draws an arrow between two points with a width w and colour

function drawArrow(x1, y1, x2, y2, w, color) {
	xd=x1-x2;
	yd=y1-y2;
	arrowlength=pow(xd*xd+yd*yd, 0.5);
	headwidth=1;
	headlength=2;
	arrowshapex=newArray(-w/2, -w/2, -w*headwidth, 0, w*headwidth, w/2, w/2);
	arrowshapey=newArray(0, arrowlength-headlength*w, arrowlength-headlength*w, arrowlength, arrowlength-headlength*w, arrowlength-headlength*w, 0);
	angle=-3.1415926/2+atan2(yd, xd);
	arrowtransformx=newArray(lengthOf(arrowshapex));
	arrowtransformy=newArray(lengthOf(arrowshapey));
	for (i=0; i<lengthOf(arrowshapex); i++) {
		x=arrowshapex[i];
		y=arrowshapey[i];
		xt=x*cos(angle)-y*sin(angle);
		yt=x*sin(angle)+y*cos(angle);
	arrowtransformx[i]=xt+x1;
		arrowtransformy[i]=yt+y1;
	}
	makeSelection("polygon", arrowtransformx, arrowtransformy);
	run("Add...", "value="+color);
}

//Macro - Correct Chromatic Aberation Tool
//=============================================================================================================================
//Global variables used
//Slice numbers
//	slice_dapi
//	slice_pi
//PI channel distort variables
//	distort_scale
//	distort_xori
//	distort_yori

macro "Correct Chromatic Aberation Action Tool -Cf88o1422Cf00o8577C00fo2722C88fo9977C000L08f8Lf8c5Lf8cb" {
	setBatchMode(true);

	//Display options
	Dialog.create("Correct Chromatic Aberation");
		Dialog.addCheckbox("Modify all images currently open", analyseall);
		Dialog.addNumber("Scale factor:", distort_scale, 5, 9, "");
		Dialog.addNumber("Origin (X):", distort_xori, 5, 9, "px");
		Dialog.addNumber("Origin (Y):", distort_yori, 5, 9, "px");
	Dialog.show();
		analyseall=Dialog.getCheckbox();
		distort_scale=Dialog.getNumber();
		distort_xori=Dialog.getNumber();
		distort_yori=Dialog.getNumber();

	//Grab images IDs for analysis
	if (analyseall==true) {
		images=newArray(nImages());
		for (n=0; n<nImages; n++) {
			selectImage(n+1);
			images[n]=getImageID();
		}
	} else {
		images=newArray(1);
		images[0]=getImageID();
	}

	//Loop through images performing modification to pi channel
	for (n=0; n<lengthOf(images); n++) {
		run("Select None");
		selectImage(images[n]);
		setSlice(slice_pi);	
		modimage=originScale(distort_xori, distort_yori, distort_scale);
		run("Select All");
		run("Copy");
		selectImage(modimage);
		close();
		selectImage(images[n]);
		setSlice(slice_pi);
		setPasteMode("Copy");
		run("Paste");
		run("Select None");
	}

	exit(""+lengthOf(images)+" images modified.");
	setBatchMode(false);
}

//Distortion Adjust Function - originScale
//-------------------------------------------------------------------------------------------------------------------------------
//requires a scalefactor of greater than 1
//x and y represent the origin of scaling
//returns the imageID of the resulting image
function originScale(x, y, scalefactor) {

	w=getWidth();
	h=getHeight();
	run("Duplicate...", "title=originScaleTemp");
	originScaleTemp=getImageID();

	//Determine anchor position and canvas size
	if (x<=w/2) {
		xanchor="Right";
		wnew=round((w-x)*2);
	} else {
		xanchor="Left";
		wnew=round(x*2);
	}
	if (y<=h/2) {
		yanchor="Bottom";
		hnew=round((h-y)*2);
	} else {
		yanchor="Top";
		hnew=round(y*2);
	}

	//Adjust canvas size to centre image
	run("Canvas Size...", "width="+wnew+" height="+hnew+" position="+yanchor+"-"+xanchor);

	//Scale the image and crop to original size
	run("Size...", "width="+wnew*scalefactor+" height="+hnew*scalefactor+" constrain interpolation=Bicubic");
	xcentre=round((wnew*scalefactor)/2);
	ycentre=round((hnew*scalefactor)/2);
	makeRectangle(xcentre-x, ycentre-y, w, h);
	run("Crop");

	//return imageID
	return getImageID();
}

//Macro - Measure K/N Signal Tool
//=============================================================================================================================
//Global variables used
//Slice numbers
//	slice_dapi
//	slice_pi
//Global variables modified
//PI channel distort variables
//	dapi_nucleus
//	dapi_kinetoplast
//	dapi_background
//	pi_nucleus
//	pi_kinetoplast
//	pi_background

var mkntolerance=200;
var excludeoutlierstdevs=3;

macro "Measure K/N Signal Action Tool -C80fo1322Cf08o8477C000L0fffL0c0fL4d4fL8c8fLcdcfLfcff" {
	setBatchMode(true);

	//Display options
	printresults=false;
	Dialog.create("Measure K/N Signal");
		Dialog.addCheckbox("Analyse all images currently open", analyseall);
		Dialog.addMessage("Analysis options:");
		Dialog.addNumber("Noise tolerance:", mkntolerance, 0, 5, "");
		Dialog.addNumber("Exclude outliers:", excludeoutlierstdevs, 0, 5, "standard deviations");
		Dialog.addCheckbox("Print data to log", printresults);
	Dialog.show();
		analyseall=Dialog.getCheckbox();
		mkntolerance=Dialog.getNumber();
		excludeoutlierstdevs=Dialog.getNumber();
		printresults=Dialog.getCheckbox();

	//Get image IDs
	setBatchMode(true);
	if (analyseall==true) {
		images=newArray(nImages());
		for (n=0; n<nImages; n++) {
			selectImage(n+1);
			images[n]=getImageID();
		}
	} else {
		images=newArray(1);
		images[0]=getImageID();
	}

	sumdmode=0;
	sumpmode=0;
	for (i=0; i<lengthOf(images); i++) {
		//Get slice modes
		run("Select None");
		setSlice(slice_dapi);
		getRawStatistics(area, mean, min, max, stdev, histogram);
		hmax=0;
		modepos=0;
		for (a=0; a<lengthOf(histogram); a++) {
			if(histogram[a]>hmax) {
				mode=a;
				hmax=histogram[a];
			}
		}
		sumdmode=sumdmode+mode+1;
		setSlice(slice_pi);
		getRawStatistics(area, mean, min, max, stdev, histogram);
		hmax=0;
		modepos=0;
		for (a=0; a<lengthOf(histogram); a++) {
			if(histogram[a]>hmax) {
				mode=a;
				hmax=histogram[a];
			}
		}
		sumpmode=sumpmode+mode+1;
	}

	//Calculate average backgrounds
	tdapi_background=sumdmode/lengthOf(images);
	tpi_background=sumpmode/lengthOf(images);

	ratioso=newArray(0);
	dvalueso=newArray(0);
	pvalueso=newArray(0);
	for (i=0; i<lengthOf(images); i++) {
		selectImage(images[i]);

		//Get values for maxima and record log2 ratios to an array
		setSlice(slice_dapi);
		run("Find Maxima...", "noise="+mkntolerance+" output=[Point Selection]");
		getSelectionCoordinates(x1, y1);
		setSlice(slice_pi);
		run("Find Maxima...", "noise="+mkntolerance+" output=[Point Selection]");
		getSelectionCoordinates(x2, y2);
		v1=newArray(lengthOf(x1)+lengthOf(x2));
		v2=newArray(lengthOf(x1)+lengthOf(x2));
		for (n=0; n<lengthOf(x1); n++) {
			setZCoordinate(slice_dapi-1);
			v1[n]=getPixel(x1[n], y1[n])-tdapi_background;
			setZCoordinate(slice_pi-1);
			v2[n]=getPixel(x1[n], y1[n])-tpi_background;
		}
		for (n=0; n<lengthOf(x2); n++) {
			setZCoordinate(slice_dapi-1);
			v1[n+lengthOf(x1)]=getPixel(x2[n], y2[n])-tdapi_background;
			setZCoordinate(slice_pi-1);
			v2[n+lengthOf(x1)]=getPixel(x2[n], y2[n])-tpi_background;
		}
		r=newArray(lengthOf(v1));
		for(n=0; n<lengthOf(v1); n++) {
			if (v1[n]<=0) {
				v1[n]=1;
			}
			if (v2[n]<=0) {
				v2[n]=1;
			}				
			r[n]=log(v1[n]/v2[n])/log(2);
		}

		//Append data to final array
		ratioso=appendArrayToArray2(ratioso, r);
		dvalueso=appendArrayToArray2(dvalueso, v1);
		pvalueso=appendArrayToArray2(pvalueso, v2);

	}

	//Iterate through ratios array and exclude outliers (due to shot noise, dead pixels, fluorescent debris, etc)
	nstdev=excludeoutlierstdevs;
	ratiosmean=arrayAverage(ratioso);
	ratiosstdev=arrayStdDev(ratioso);
	includearr=newArray(lengthOf(ratioso));
	includecount=0;
	for (n=0; n<lengthOf(ratioso); n++) {
		if (dvalueso[n]!=0 && pvalueso[n]!=0) {
			if (ratioso[n]>ratiosmean-nstdev*ratiosstdev && ratioso[n]<ratiosmean+nstdev*ratiosstdev) {
				includearr[n]=1;
				includecount++;
			}
		}
	}
	ratios=newArray(includecount);
	dvalues=newArray(includecount);
	pvalues=newArray(includecount);
	index=0;
	for (n=0; n<lengthOf(ratioso); n++) {
		if (includearr[n]==1) {
			ratios[index]=ratioso[n];
			dvalues[index]=dvalueso[n];
			pvalues[index]=pvalueso[n];
			index++;
		}
	}

	//1D K means clustering of ratio data
	//Calculate mean and stdev of data
	sum=0;
	for (n=0; n<lengthOf(ratios); n++) {
		sum=sum+ratios[n];
	}
	mean=sum/lengthOf(ratios);
	sumsq=0;
	for (n=0; n<lengthOf(ratios); n++) {
		dif=ratios[n]-mean;
		difsq=dif*dif;
		sumsq=sumsq+difsq;
	}
	stdevsq=sumsq/lengthOf(ratios);
	stdev=pow(stdevsq, 0.5);

	//Starting cluster values
	cluster1=mean-stdev;
	cluster2=mean+stdev;

	//K means clustering
	nchanges=lengthOf(ratios);
	group=newArray(lengthOf(ratios));
	while (nchanges>0) {
		nchanges=0;
		for (n=0; n<lengthOf(ratios); n++) {
			dist1=abs(cluster1-ratios[n]);
			dist2=abs(cluster2-ratios[n]);
			oldgroup=group[n];
			if (dist1<dist2) {
				group[n]=1;
			} else {
				group[n]=2;
			}
			if (group[n]!=oldgroup) {
				nchanges++;
			}
		}
		sum1=0;
		sum2=0;
		n1=0;
		n2=0;
		for (n=0; n<lengthOf(ratios); n++) {
			if (group[n]==1) {
				sum1=sum1+ratios[n];
				n1++;
			} else {
				sum2=sum2+ratios[n];
				n2++;
			}
		}
		if (n1!=0) {
			cluster1=sum1/n1;
		}
		if (n2!=0) {
			cluster2=sum2/n2;
		}
	}

	//Calculate mean dapi and pi intensities for all points within each cluster
	d1sum=0;
	p1sum=0;
	d2sum=0;
	p2sum=0;
	for (n=0; n<lengthOf(ratios); n++) {
		if (group[n]==1) {
			d1sum=d1sum+dvalues[n];
			p1sum=p1sum+pvalues[n];
		} else if (group[n]==2) {
			d2sum=d2sum+dvalues[n];
			p2sum=p2sum+pvalues[n];
		}
	}
	d1mean=d1sum/n1;
	p1mean=p1sum/n1;
	d2mean=d2sum/n1;
	p2mean=p2sum/n1;

	if (printresults==true) {
		print("DAPI value,PI value,Log Ratio,Cluster");
		for (n=0; n<lengthOf(ratios); n++) {
			print(dvalues[n]+","+pvalues[n]+","+ratios[n]+","+group[n]);
		}
	}

	//Assign the two groups to kin and nuc according to ratio value
	if (cluster1<cluster2) {
		tdapi_nucleus=d1mean;
		tpi_nucleus=p1mean;
		tdapi_kinetoplast=d2mean;
		tpi_kinetoplast=p2mean;
	} else {
		tdapi_nucleus=d2mean;
		tpi_nucleus=p2mean;
		tdapi_kinetoplast=d1mean;
		tpi_kinetoplast=p1mean;
	}
	ratioratio=(tdapi_nucleus/tpi_nucleus)/(tdapi_kinetoplast/tpi_kinetoplast);

	if (printresults==true) {
		print("DAPI Background", "PI Background");
		print(tdapi_background, tpi_background);
		print("DAPI Nucleus", "PI Nucleus");
		print(tdapi_nucleus, tpi_nucleus);
		print("DAPI Kinetoplast", "PI Kinetoplast");
		print(tdapi_kinetoplast, tpi_kinetoplast);
		print("Points in cluster 1", "Points in cluster 2");
		print(n1, n2);
		print("Ratio ratio");
		print(ratioratio);
	}

	Dialog.create("Automatic Intensity Sampling Tool");
		Dialog.addMessage("Nucleus DAPI: "+tdapi_nucleus);
		Dialog.addMessage("Nucleus PI: "+tpi_nucleus);
		Dialog.addMessage("Kinetoplast DAPI: "+tdapi_kinetoplast);
		Dialog.addMessage("Kinetoplast PI: "+tpi_kinetoplast);
		Dialog.addMessage("Background DAPI: "+tdapi_background);
		Dialog.addMessage("Background PI: "+tpi_background);
		Dialog.addMessage("Points in each cluster: ("+n1+", "+n2+")");
		Dialog.addMessage("Ratio ratio: "+ratioratio);
	Dialog.show();
	dapi_nucleus=tdapi_nucleus;
	dapi_kinetoplast=tdapi_kinetoplast;
	dapi_background=tdapi_background;
	pi_nucleus=tpi_nucleus;
	pi_kinetoplast=tpi_kinetoplast;
	pi_background=tpi_background;
}

//Intensity test function - appendArrayToArray2
//-----------------------------------------------------------------------------------------------------------------------------
//Joins two arrays together and returns the result
function appendArrayToArray2(array1, array2) {
	result=newArray(lengthOf(array1)+lengthOf(array2));
	for (i=0; i<lengthOf(array1); i++) {
		result[i]=array1[i];
	}
	for (i=0; i<lengthOf(array2); i++) {
		result[i+lengthOf(array1)]=array2[i];
	}
	return result;
}

//Intensity test function - arrayAverage
//-----------------------------------------------------------------------------------------------------------------------------
//Returns the mean of an array
function arrayAverage(array) {
	sum=0;
	for (a=0; a<lengthOf(array); a++) {
		sum+=array[a];
	}
	return sum/lengthOf(array);
}

//Intensity test function - arrayStdDev
//-----------------------------------------------------------------------------------------------------------------------------
//Returns the standard deviation of an array
function arrayStdDev(array) {
	sumsqd=0;
	mean=arrayAverage(array);
	for (a=0; a<lengthOf(array); a++) {
		sumsqd+=pow(array[a]-mean, 2);
	}
	return pow(sumsqd/lengthOf(array), 0.5);
}

//Macro - Colour Deconvolution Tool
//=============================================================================================================================
//Global variables used
//Slice numbers
//	slice_dapi
//	slice_pi
//PI channel distort variables
//	dapi_nucleus
//	dapi_kinetoplast
//	dapi_background
//	pi_nucleus
//	pi_kinetoplast
//	pi_background

var cdsubtractbackground=false;
var cdrollingballradius=15;

macro "Colour Deconvolution Action Tool -C0ffo1622Cf80o8877C000L08f8Lf8c5Lf8cb" {
	setBatchMode(true);

	//Display options
	Dialog.create("Colour Deconvolution");
		Dialog.addCheckbox("Modify all images currently open", analyseall);
		Dialog.addNumber("Nucleus DAPI:", dapi_nucleus, 3, 9, "");
		Dialog.addNumber("Nucleus PI:", pi_nucleus, 3, 9, "");
		Dialog.addNumber("Kinetoplast DAPI:", dapi_kinetoplast, 3, 9, "");
		Dialog.addNumber("Kinetoplast PI:", pi_kinetoplast, 3, 9, "");
		Dialog.addNumber("Background DAPI:", dapi_background, 3, 9, "");
		Dialog.addNumber("Background PI:", pi_background, 3, 9, "");
		Dialog.addCheckbox("Subtract background following processing", cdsubtractbackground);
		Dialog.addNumber("Rolling ball radius:", cdrollingballradius, 0, 9, "px");
	Dialog.show();
		analyseall=Dialog.getCheckbox();
		dapi_nucleus=Dialog.getNumber();
		pi_nucleus=Dialog.getNumber();
		dapi_kinetoplast=Dialog.getNumber();
		pi_kinetoplast=Dialog.getNumber();
		dapi_background=Dialog.getNumber();
		pi_background=Dialog.getNumber();
		cdsubtractbackground=Dialog.getCheckbox();
		cdrollingballradius=Dialog.getNumber();

	//Grab images IDs for analysis
	if (analyseall==true) {
		images=newArray(nImages());
		for (n=0; n<nImages; n++) {
			selectImage(n+1);
			images[n]=getImageID();
			run("Select None");
		}
	} else {
		images=newArray(1);
		images[0]=getImageID();
		run("Select None");
	}

	//Normalise and subtract bg from kinetoplast and nucleus vectors
	kinsum=dapi_kinetoplast+pi_kinetoplast;
	nucsum=dapi_nucleus+pi_nucleus;
	ndapi_kinetoplast=(dapi_kinetoplast)/kinsum;
	npi_kinetoplast=(pi_kinetoplast)/kinsum;
	ndapi_nucleus=(dapi_nucleus)/nucsum;
	npi_nucleus=(pi_nucleus)/nucsum;
	mdetinv=1/(ndapi_kinetoplast*npi_nucleus-ndapi_nucleus*npi_kinetoplast);

	//Loop through all images and perform deconvolution
	for (i=0; i<lengthOf(images); i++) {

		selectImage(images[i]);
		run("Select None");
		setSlice(slice_dapi);
		run("Duplicate...", "dapi1");
		rename("dapi1");
		dapi1=getImageID();
		run("Duplicate...", "dapi2");
		rename("dapi2");
		dapi2=getImageID();
		selectImage(images[i]);
		setSlice(slice_pi);
		run("Duplicate...", "pi1");
		rename("pi1");
		pi1=getImageID();
		run("Duplicate...", "pi2");
		rename("pi2");
		pi2=getImageID();

		selectImage(dapi1);
		run("32-bit");
		run("Multiply...", "value="+npi_nucleus);
		selectImage(dapi2);
		run("32-bit");
		run("Multiply...", "value="+-npi_kinetoplast);
		selectImage(pi1);
		run("32-bit");
		run("Multiply...", "value="+-ndapi_nucleus);
		selectImage(pi2);
		run("32-bit");
		run("Multiply...", "value="+ndapi_kinetoplast);

		imageCalculator("Add create 16-bit", "dapi1","pi1");
		rename("pi");
		pi=getImageID();
		run("Multiply...", "value="+mdetinv);
		selectImage(dapi1);
		close();
		selectImage(pi1);
		close();
		imageCalculator("Add create 16-bit", "dapi2","pi2");
		rename("dapi");
		dapi=getImageID();
		run("Multiply...", "value="+mdetinv);
		selectImage(dapi2);
		close();
		selectImage(pi2);
		close();

		run("Options...", "iterations=1 black edm=Overwrite count=1");
		selectImage(pi);
		run("Select None");
		if (cdsubtractbackground==true) {
			run("Subtract Background...", "rolling="+cdrollingballradius);
		}
		getRawStatistics(area, mean, min, max, stdev, histogram);
		maxh=0;
		for (m=0; m<lengthOf(histogram); m++) {
			if (histogram[m]>maxh) {
				medv=m;
				maxh=histogram[m];
			}
		}
		median=(max-min)*(medv/255)+min;
		run("Subtract...", "value="+median);
		selectImage(dapi);
		run("Select None");
		if (cdsubtractbackground==true) {
			run("Subtract Background...", "rolling="+cdrollingballradius);
		}
		getRawStatistics(area, mean, min, max, stdev, histogram);
		maxh=0;
		for (m=0; m<lengthOf(histogram); m++) {
			if (histogram[m]>maxh) {
				medv=m;
				maxh=histogram[m];
			}
		}
		median=(max-min)*(medv/255)+min;
		run("Subtract...", "value="+median);

		//Copy images into pi and dapi channels of stack
		selectImage(dapi);
		run("Select All");
		run("Copy");
		close();
		selectImage(images[i]);
		setSlice(slice_dapi);
		setPasteMode("Copy");
		run("Paste");
		selectImage(pi);
		run("Select All");
		run("Copy");
		close();
		selectImage(images[i]);
		setSlice(slice_pi);
		run("Paste");
	}

	exit(""+lengthOf(images)+" images modified.");
	setBatchMode(false);

}

//Macro - Cell Analysis Tool
//=============================================================================================================================
//Global variables used
//Slice numbers
//	slice_phase
//	slice_dapi
//	slice_pi
//Image scale
//	image_scale

var displaymasks=true;
var rollingballradius=10;
var cellautomanual="Automatic";
var cellthresholdvalue=6000;
var cellautothreshtype="Triangle";
var cellthresholdvalue=6000;
var cellminarea=200;
var trimbranches=true;
var minbranch=10;
var threshtype="Manual";
var autothreshtype="MaxEntropy";
var nthresholdvalue=600;
var kthresholdvalue=500;
var nucminarea=50;
var kinminarea=5;

macro "Cell Analysis Action Tool -C0ffo1322Cf80o8477C000T0f081T4f082T8f083" {
	setBatchMode(true);
	getDateAndTime(year, month, dayOfWeek, dayOfMonth, hour, minute, second, msec);
	starttime=second+100*minute+100*100*hour+100*100*100*dayOfMonth+100*100*100*100*month+100*100*100*100*100*year;
	starttime=""+year+"_"+month+"_"+dayOfMonth+"_"+hour+"_"+minute+"_"+second;

	//User interface
	threshchoice=newArray("Default", "Huang", "Intermodes", "IsoData", "Li", "MaxEntropy", "Minimum", "Moments", "Otsu", "Percentile", "RenyiEntropy", "Shanbhag", "Triangle", "Yen");
	automanual=newArray("Manual", "Automatic");
	Dialog.create("K/N Count")
		Dialog.addCheckbox("Process all images currently open", analyseall);
		//Dialog.addCheckbox("Give a visual display of results", displaymasks);
		Dialog.addMessage("Cell thresholding options:");
		Dialog.addNumber("Rolling ball radius:", rollingballradius, 0, 5, "px");
		Dialog.addChoice("Thresholding type:", automanual, automanual[getArrayIndex(automanual, cellautomanual)]);
		Dialog.addNumber("Manual cell threshold value:", cellthresholdvalue, 0, 5, "");
		Dialog.addChoice("Automatic threshold type:", threshchoice, threshchoice[getArrayIndex(threshchoice, cellautothreshtype)]);
		Dialog.addNumber("Minimum cell area:", cellminarea, 2, 5, "px^2");
		Dialog.addCheckbox("Trim short skeleton branches", trimbranches);
		Dialog.addNumber("Minimum skeleton branch length:", minbranch, 2, 5, "px");
		Dialog.addMessage("Nucleus and Kintoplast thresholding options:");
		Dialog.addChoice("Thresholding type:", automanual, automanual[getArrayIndex(automanual, threshtype)]);
		Dialog.addNumber("Manual nucleus threshold value:", nthresholdvalue, 0, 5, "");
		Dialog.addNumber("Manual kinetoplast threshold value:", kthresholdvalue, 0, 5, "");
		Dialog.addChoice("Automatic threshold type:", threshchoice, threshchoice[getArrayIndex(threshchoice, autothreshtype)]);
		Dialog.addNumber("Minimum nucleus area:", nucminarea, 2, 5, "px^2");
		Dialog.addNumber("Minimum kinetoplast area:", kinminarea, 2, 5, "px^2");
	Dialog.show();
		analyseall=Dialog.getCheckbox();
		//displaymasks=Dialog.getCheckbox();
		rollingballradius=Dialog.getNumber();
		cellautomanual=Dialog.getChoice();
		cellthresholdvalue=Dialog.getNumber();
		cellautothreshtype=Dialog.getChoice();
		cellminarea=Dialog.getNumber();
		trimbranches=Dialog.getCheckbox();
		minbranch=Dialog.getNumber();
		threshtype=Dialog.getChoice();
		nthresholdvalue=Dialog.getNumber();
		kthresholdvalue=Dialog.getNumber();
		autothreshtype=Dialog.getChoice();
		nucminarea=Dialog.getNumber();
		kinminarea=Dialog.getNumber();


	//Grab images IDs for analysis
	if (analyseall==true) {
		images=newArray(nImages());
		for (n=0; n<nImages; n++) {
			selectImage(n+1);
			images[n]=getImageID();
		}
	} else {
		images=newArray(1);
		images[0]=getImageID();
	}

	run("Options...", "iterations=1 black edm=Overwrite count=1");
	cellno=0;
	for (i=0; i<lengthOf(images); i++) {
		selectImage(images[i]);
		rename("Analysed_"+starttime+"_"+i);
		run("Select None");
		w=getWidth();
		h=getHeight();

		//Detect and size filter cell masks
		selectImage(images[i]);
		setSlice(slice_phase);
		run("Duplicate...", "title=phase");
		phase=getImageID();
		run("Subtract Background...", "rolling="+rollingballradius+" light");
		setAutoThreshold(cellautothreshtype);
		run("Convert to Mask");
		run("Dilate");
		//run("Dilate");
		run("Erode");
		run("Erode");
		run("Dilate");
		for (x=0; x<w; x++) {
			for (y=0; y<h; y++) {
				if (getPixel(x, y)==255) {
					doWand(x, y);
					getRawStatistics(area);
					if (area<cellminarea) {
						setColor(0);
						fill();
					} else {
						setColor(128);
						fill();
					}
				}
			}
		}
		run("Select None");

		//Generate skeletonised cell shapes
		selectImage(phase);
		run("Duplicate...", "title=skeleton");
		skeleton=getImageID();
		run("Make Binary");
		run("Skeletonize");
		//Trim branches
		if (trimbranches==true) {
			//Tidy 4-connected branch points which have ambiguity when trimming branches
			for (x=1; x<w-1; x++) {
				for (y=1; y<h-1; y++) {
					if (getPixel(x, y)!=0) {
						count=0;
						xcarr=newArray(0, 1, -1, 0);
						ycarr=newArray(1, 0, 0, -1);
						for (ic=0; ic<lengthOf(xcarr); ic++) {
							if (getPixel(x+xcarr[ic], y+ycarr[ic])!=0) {
								count++;
							}
						}
						if (count>2) {
							setPixel(x, y, 0);
						}
					}
				}
			}
			//Find all white pixels in the skeleton image with only 1 neighbour
			//These pixels are the branch ends
			nterm=0;
			for (x=1; x<w-1; x++) {
				for (y=1; y<h-1; y++) {
					if (getPixel(x, y)==255) {
						count=0;
						for (a=-1; a<2; a++) {
							for (b=-1; b<2; b++) {
								if ((a!=0 || b!=0) && getPixel(x+a, y+b)==255) {
									count++;
								}
							}
						}
						if (count==1) {
							nterm++;
						}
					}
				}
			}
			termx=newArray(nterm);
			termy=newArray(nterm);
			tnum=0;
			nterm=0;
			for (x=1; x<w-1; x++) {
				for (y=1; y<h-1; y++) {
					if (getPixel(x, y)==255) {
						count=0;
						for (a=-1; a<2; a++) {
							for (b=-1; b<2; b++) {
								if ((a!=0 || b!=0) && getPixel(x+a, y+b)==255) {
									count++;
								}
							}
						}
						if (count==1) {
							termx[tnum]=x;
							termy[tnum]=y;
							tnum++;
						}
					}
				}
			}
			//Start from each branch end
			//Walk along the branch counting its length
			for (n=0; n<lengthOf(termx); n++) {
				cx=termx[n];
				cy=termy[n];
				visitx=newArray(minbranch+2);
				visitx[0]=cx;
				visity=newArray(minbranch+2);
				visity[0]=cy;
				finished=0;
				leng=0;
				while (finished==0) {
					finished=1;
					count=0;
					for (a=-1; a<2; a++) {
						for (b=-1; b<2; b++) {
							if ((a!=0 || b!=0) && getPixel(cx+a, cy+b)==255) {
								visited=0;
								for (m=0; m<minbranch; m++) {
									if (cx+a==visitx[m] && cy+b==visity[m]) {
										visited=1;
									}
								}
								if (visited==0) {
									count++;
									ctx=cx+a;
									cty=cy+b;
								}
							}
						}
					}
					leng++;
					if (count!=1 || leng>minbranch+1) {
						finished=1;
					} else {
						finished=0;
						visitx[leng]=cx;
						visity[leng]=cy;
						cx=ctx;
						cy=cty;
					}
				}
				//Erase the branch if under the minimum length
				if (leng<=minbranch) {
					for (m=0; m<lengthOf(visitx); m++) {
						setPixel(visitx[m], visity[m], 0);
					}
				}
			}
			//Tidy branch stubs which can form following pruning
			for (x=0; x<w; x++) {
				for (y=0; y<w; y++) {
					if (getPixel(x, y)!=0) {
						cvalid=0;
						for (ori=0; ori<3; ori++) {
							//0=must be blank, 1=must be filled, -1=can be either
							if (ori==0) {
								tav=newArray(
									0, 0, 0,
									-1, 1, 0,
									-1, 1 ,1
								);
							} else if (ori==1) {
								tav=newArray(
									-1, -1, 0,
									1, 1, 0,
									1, 0 ,0
								);
							} else if (ori==2) {
								tav=newArray(
									0, 1, -1,
									0, 1, 1,
									0, 0 ,0
								);
							}
							for (od=0; od<4; od++) {
								if (od==0) {
									c=1;
									d=1;
								} else if (od==1) {
									c=-1;
									d=1;
								} else if (od==2) {
									c=1;
									d=-1;
								} else if (od==3) {
									c=-1;
									d=-1;
								}
								valid=true;
								for (a=-1; a<2; a++) {
									for (b=-1; b<2; b++) {
										if (tav[(a*c+1)+3*(b*d+1)]==0 && getPixel(x+a, y+b)!=0) {
											valid=false;
										} else if (tav[(a*c+1)+3*(b*d+1)]==1 && getPixel(x+a, y+b)==0) {
											valid=false;
										}
									}
								}
								if (valid==true) {
									cvalid++;
								}
							}
						}
						if (cvalid!=0) {
							setPixel(x, y, 0);
						}
					}
				}
			}
		}

		//Create distance maps and medial axis transform
		selectImage(phase);
		run("Duplicate...", "title=mat");
		mat=getImageID();
		run("Make Binary");
		run("Distance Map");
		selectImage(skeleton);
		run("Select All");
		run("Copy");
		run("Select None");
		selectImage(mat);
		setPasteMode("Transparent-white");
		run("Paste");
		setPasteMode("Copy");
		run("Select None");

		//Detect and size filter kin masks
		selectImage(images[i]);
		setSlice(slice_pi);
		run("Duplicate...", "title=kin");
		kin=getImageID();
		if (threshtype=="Automatic") {
			setAutoThreshold(autothreshtype+" dark");
		} else {
			setThreshold(kthresholdvalue, pow(2, 16));
		}
		run("Convert to Mask");
		//run("Watershed");
		for (x=0; x<w; x++) {
			for (y=0; y<h; y++) {
				if (getPixel(x, y)==255) {
					doWand(x, y);
					getRawStatistics(area);
					if (area<kinminarea) {
						setColor(0);
						fill();
					} else {
						setColor(128);
						fill();
					}
				}
			}
		}
		run("Select None");

		//Detect and size filter nuc masks
		selectImage(images[i]);
		setSlice(slice_dapi);
		run("Duplicate...", "title=nuc");
		nuc=getImageID();
		if (threshtype=="Automatic") {
			setAutoThreshold(autothreshtype+" dark");
		} else {
			setThreshold(nthresholdvalue, pow(2, 16));
		}
		run("Convert to Mask");
		//run("Watershed");
		for (x=0; x<w; x++) {
			for (y=0; y<h; y++) {
				if (getPixel(x, y)==255) {
					doWand(x, y);
					getRawStatistics(area);
					if (area<nucminarea) {
						setColor(0);
						fill();
					} else {
						setColor(128);
						fill();
					}
				}
			}
		}
		run("Select None");

		//Compile images into a single stack
		newImage("AnalysedMasks_"+starttime+"_"+i, "8-bit Black", w, h, 5);
		masks=getImageID();
		setPasteMode("Copy");
		selectImage(phase);
		run("Select All");
		run("Copy");
		close();
		selectImage(masks);
		setSlice(1);
		run("Paste");
		selectImage(kin);
		run("Select All");
		run("Copy");
		close();
		selectImage(masks);
		setSlice(2);
		run("Paste");
		selectImage(nuc);
		run("Select All");
		run("Copy");
		close();
		selectImage(masks);
		setSlice(3);
		run("Paste");
		selectImage(mat);
		run("Select All");
		run("Copy");
		close();
		selectImage(masks);
		setSlice(4);
		run("Paste");
		selectImage(skeleton);
		close();

		selectImage(images[i]);
		run("Duplicate...", "title=temp duplicate");
		temp=getImageID();

		//Force column order for some results
		setResult("Cell Area (um2)", 0, 0);
		setResult("Total Kinetoplast DNA", 0, 0);
		setResult("Total Nuclear DNA", 0, 0);
		setResult("Touching Edge (bool)", 0, 0);
		setResult("Kinetoplast Number", 0, 0);
		setResult("Nucleus Number", 0, 0);
		setResult("Skeleton length", 0, 0);
		setResult("Skeleton termini", 0, 0);
		setResult("Skeleton branches", 0, 0);
		setResult("Cell length (um)", 0, 0);
		setResult("Cell width (um)", 0, 0);
		setResult("Image No", 0, 0);
		setResult("Year", 0, 0);
		setResult("Month", 0, 0);
		setResult("Day", 0, 0);
		setResult("Hour", 0, 0);
		setResult("Minute", 0, 0);
		setResult("Second", 0, 0);
		setResult("Cell OX", 0, 0);
		setResult("Cell OY", 0, 0);
		setResult("Cell X", 0, 0);
		setResult("Cell Y", 0, 0);
		setResult("Cell W", 0, 0);
		setResult("Cell H", 0, 0);

		//Perform particle analysis
		selectImage(masks);
		setSlice(1);
		run("Select None");
		startcellno=cellno;
		for (x=0; x<w; x++) {
			for (y=0; y<h; y++) {
				//Find cells and measure properties
				if (getPixel(x, y)==128) {
					doWand(x, y);
					getSelectionBounds(sx, sy, sw, sh);
					setColor(196);
					fill();
					getRawStatistics(area);
					setResult("Image No", cellno, i);
					setResult("Year", cellno, year);
					setResult("Month", cellno, month);
					setResult("Day", cellno, dayOfMonth);
					setResult("Hour", cellno, hour);
					setResult("Minute", cellno, minute);
					setResult("Second", cellno, second);
					setResult("Cell OX", cellno, x);
					setResult("Cell OY", cellno, y);
					setResult("Cell X", cellno, sx);
					setResult("Cell Y", cellno, sy);
					setResult("Cell W", cellno, sw);
					setResult("Cell H", cellno, sh);
					setResult("Cell Area (um2)", cellno, area/(image_scale*image_scale));
					roiManager("Reset");
					roiManager("Add");
					selectImage(temp);
					roiManager("Select", 0);
					setSlice(slice_pi);
					getRawStatistics(area, mean);
					setResult("Total Kinetoplast DNA", cellno, area*mean);
					setSlice(slice_dapi);
					getRawStatistics(area, mean);
					setResult("Total Nuclear DNA", cellno, area*mean);
					selectImage(masks);
					borderdist=2;
					if (sx<borderdist || sy<borderdist || sx+sw>=w-borderdist || sy+sh>=h-borderdist) {
						setResult("Touching Edge (bool)", cellno, 1);
					} else {
						setResult("Touching Edge (bool)", cellno, 0);
					}
					//Analyse kinetoplasts
					particleno1=0;
					for (xa=sx; xa<sx+sw; xa++) {
						for (ya=sy; ya<sy+sh; ya++) {
							setSlice(1);
							if (getPixel(xa, ya)==196) {
								setSlice(2);
								if (getPixel(xa, ya)==128) {
									doWand(xa, ya);
									roiManager("Reset");
									roiManager("Add");
									setColor(255);
									fill();
									selectImage(temp);
									roiManager("Select", 0);
									setSlice(slice_pi);
									getRawStatistics(area, mean, min, max, 	stdev, histogram);
									setResult("Kinetoplast "+particleno1+1+" Area (um2)", cellno, area/(image_scale*image_scale));
									setResult("Kinetoplast "+particleno1+1+" Sum Intensity", cellno, area*mean);
									//setResult("Kinetoplast "+particleno1+1+" Mean", cellno, mean);
									//setResult("Kinetoplast "+particleno1+1+" Min", cellno, min);
									//setResult("Kinetoplast "+particleno1+1+" Max", cellno, max);
									//setResult("Kinetoplast "+particleno1+1+" Stdev", cellno, stdev);
									setResult("Kinetoplast "+particleno1+1+" Centroid X", cellno, getCentroidX());
									setResult("Kinetoplast "+particleno1+1+" Centroid Y", cellno, getCentroidY());
									particleno1++;
									selectImage(masks);
								}
							}
						}
					}
					setResult("Kinetoplast Number", cellno, particleno1);
					//Analyse nuclei
					particleno2=0;
					for (xa=sx; xa<sx+sw; xa++) {
						for (ya=sy; ya<sy+sh; ya++) {
							setSlice(1);
							if (getPixel(xa, ya)==196) {
								setSlice(3);
								if (getPixel(xa, ya)==128) {
									doWand(xa, ya);
									roiManager("Reset");
									roiManager("Add");
									setColor(255);
									fill();
									selectImage(temp);
									roiManager("Select", 0);
									setSlice(slice_dapi);
									getRawStatistics(area, mean, min, max, stdev, histogram);
									setResult("Nucleus "+particleno2+1+" Area (um2)", cellno, area/(image_scale*image_scale));
									setResult("Nucleus "+particleno2+1+" Sum Intensity", cellno, area*mean);
									//setResult("Nucleus "+particleno2+1+" Mean", cellno, mean);
									//setResult("Nucleus "+particleno2+1+" Min", cellno, min);
									//setResult("Nucleus "+particleno2+1+" Max", cellno, max);
									//setResult("Nucleus "+particleno2+1+" Stdev", cellno, stdev);
									setResult("Nucleus "+particleno2+1+" Centroid X", cellno, getCentroidX());
									setResult("Nucleus "+particleno2+1+" Centroid Y", cellno, getCentroidY());
									particleno2++;
									selectImage(masks);
								}
							}
						}
					}
					setResult("Nucleus Number", cellno, particleno2);
					//Analyse skeleton
					selectImage(masks);
					setSlice(4);
					nskeleton=0;
					ntermini=0;
					nbranch=0;
					sumtermini=0;
					maxwidth=0;
					for (xa=sx; xa<sx+sw; xa++) {
						for (ya=sy; ya<sy+sh; ya++) {
							setSlice(1);
							if (getPixel(xa, ya)==196) {
								setSlice(4);
								if (getPixel(xa, ya)!=0) {
									nskeleton++;
									count=0;
									for (a=-1; a<2; a++) {
										for (b=-1; b<2; b++) {
											if ((a!=0 || b!=0) && getPixel(xa+a, ya+b)!=0) {
												count++;
											}
										}
									}
									if (count==1) {
										ntermini++;
										sumtermini=sumtermini+getPixel(xa, ya);
										tx=xa;
										ty=ya;
									}
									if (count>=3) {
										nbranch++;
									}
									maxwidth=maxOf(maxwidth, getPixel(xa, ya));
								}
							}
						}
					}
					setResult("Skeleton length", cellno, nskeleton);
					setResult("Skeleton termini", cellno, ntermini);
					setResult("Skeleton branches", cellno, nbranch);
					setResult("Cell length (um)", cellno, (nskeleton+sumtermini)/image_scale);
					setResult("Cell width (um)", cellno, maxwidth/image_scale);
					//Analyse KN position relative to skeleton
					selectImage(masks);
					setSlice(4);
					if (ntermini==2 && nbranch==0) {
						skelex=newArray(nskeleton);
						skeley=newArray(nskeleton);
						skeled=newArray(nskeleton);
						xc=tx;
						yc=ty;
						skelepx=0;
						while (skelepx<nskeleton) {
							for (a=-1; a<2; a++) {
								for (b=-1; b<2; b++) {
									if (getPixel(xc+a, yc+b)!=0 && (a!=0 || b!=0)) {
										visited=0;
										for (j=0; j<skelepx; j++) {
											if (xc+a==skelex[j] && yc+b==skeley[j]) {
												visited++;
											}
										}
										if (visited==0) {
											ar=a;
											br=b;
										}
									}
								}
							}
							skelex[skelepx]=xc;
							skeley[skelepx]=yc;
							skeled[skelepx]=skelepx;
							xc=xc+ar;
							yc=yc+br;
							skelepx++;
						}
						for (j=0; j<getResult("Nucleus Number", cellno); j++) {
							termdist=0;
							mindist=1000000;
							for (k=0; k<nskeleton; k++) {
								dx=skelex[k]-getResult("Nucleus "+j+1+" Centroid X", cellno);
								dy=skeley[k]-getResult("Nucleus "+j+1+" Centroid Y", cellno);
								dist=pow(dx*dx+dy*dy, 0.5);
								if (dist<mindist) {
									termdist=skeled[k];
									mindist=dist;
								}
							}
							setResult("Nucleus "+j+1+" Terminus Distance (um)", cellno, termdist/image_scale);
							setResult("Nucleus "+j+1+" Skeleton Distance (um)", cellno, mindist/image_scale);
						}
						for (j=0; j<getResult("Kinetoplast Number", cellno); j++) {
							termdist=0;
							mindist=1000000;
							for (k=0; k<nskeleton; k++) {
								dx=skelex[k]-getResult("Kinetoplast "+j+1+" Centroid X", cellno);
								dy=skeley[k]-getResult("Kinetoplast "+j+1+" Centroid Y", cellno);
								dist=pow(dx*dx+dy*dy, 0.5);
								if (dist<mindist) {
									termdist=skeled[k];
									mindist=dist;
								}
							}
							setResult("Kinetoplast "+j+1+" Terminus Distance (um)", cellno, termdist/image_scale);
							setResult("Kinetoplast "+j+1+" Skeleton Distance (um)", cellno, mindist/image_scale);
						}
					}

					//Mark cell as analysed
					cellno++;
					selectImage(masks);
					setSlice(5);
					setColor(255);
					drawString(cellno, x, y+15);
					drawString(particleno1+"K"+particleno2+"N", x, y+30);
					setSlice(1);
					doWand(x, y);
					setColor(255);
					fill();
					updateResults();
				}
			}
		}
		selectImage(temp);
		close();
		if (displaymasks==true) {
			//Display mask image
			selectImage(masks);
			run("Select None");
			setBatchMode(false);
			run("Duplicate...", "title=AnalysedMasks_"+starttime+"_"+i+" duplicate");
			run("Make Composite", "display=Composite");
			setBatchMode(true);
		}
		selectImage(masks);
		close();
	}
	exit(""+lengthOf(images)+" images analysed.");

	setBatchMode(false);
}

//DAPI counter function - appendToArray
//-----------------------------------------------------------------------------------------------------------------------------
//Adds a value to the end of an array and returns the resulting array

function appendToArray(value, array) {
	result=newArray(lengthOf(array)+1);
	for (a=0; a<lengthOf(array); a++) {
		result[a]=array[a];
	}
	result[lengthOf(result)-1]=value;
	return result;
}

//DAPI counter function - getCentroidX and Y
//-----------------------------------------------------------------------------------------------------------------------------
//Calculates the centroid of the current selection

function getCentroidX() {
	getSelectionCoordinates(x, y);
	xsum=0;
	for (i=0; i<lengthOf(x); i++) {
		xsum=xsum+x[i];
	}
	return xsum/lengthOf(x);
}

function getCentroidY() {
	getSelectionCoordinates(x, y);
	ysum=0;
	for (i=0; i<lengthOf(y); i++) {
		ysum=ysum+y[i];
	}
	return ysum/lengthOf(y);
}

//DAPI counter function - Get array index
//-----------------------------------------------------------------------------------------------------------------------------
//Gets the index of the first occurence of a value within an array

function getArrayIndex(array, value) {
	index=0;
	for (a=lengthOf(array)-1; a>=0; a--) {
		if (array[a]==value) {
			index=a;
		}
	}
	return index;
}

//Macro - K/N Count Summary
//=============================================================================================================================
//Requires data in data table

macro "K/N Count Summary Action Tool -C000T08081T4808KT88081Tc808NT0f081T4f08KT8f082Tcf08N" {

	exctouchedge=true;
	excbranchskele=true;
	Dialog.create("K/N Count Summary");
		Dialog.addCheckbox("Exclude cells touching image edge", exctouchedge);
		Dialog.addCheckbox("Exclude cells with branched skeletons", excbranchskele);
	Dialog.show();
		exctouchedge=Dialog.getCheckbox();
		excbranchskele=Dialog.getCheckbox();

	c1k1n=0;
	c1k2n=0;
	c2k1n=0;
	c2k2n=0;
	c1k0n=0;
	notcell=0;
	other=0;
	touchedge=0;
	branched=0;
	for (n=0; n<nResults; n++) {
		nuc=getResult("Nucleus Number", n);
		kin=getResult("Kinetoplast Number", n);
		if (getResult("Touching Edge (bool)", n)==0 || exctouchedge==false) {
			if ((getResult("Skeleton termini", n)==2 && getResult("Skeleton branches", n)==0) || excbranchskele==false) {
				if (kin==1 && nuc==1) {
					c1k1n++;
				} else if (kin==1 && nuc==2) {
					c1k2n++;
				} else if (kin==2 && nuc==1) {
					c2k1n++;
				} else if (kin==2 && nuc==2) {
					c2k2n++;
				} else if (kin==1 && nuc==0) {
					c1k0n++;
				} else if (kin==0 && nuc==0) {
					notcell++;
				} else {
					other++;
				}
			} else {
				branched++;
			}
		} else {
			touchedge++;
		}
	}
	totalcells=c1k1n+c1k2n+c2k1n+c2k2n+c1k0n+other;
	Dialog.create("DAPI Count Summary");
		Dialog.addMessage("Total cells: "+c1k1n+c1k2n+c2k1n+c2k2n+c1k0n+other);
		Dialog.addMessage("1K1N: "+100*c1k1n/totalcells+"% ("+c1k1n+")");
		Dialog.addMessage("1K2N: "+100*c1k2n/totalcells+"% ("+c1k2n+")");
		Dialog.addMessage("2K1N: "+100*c2k1n/totalcells+"% ("+c2k1n+")");
		Dialog.addMessage("2K2N: "+100*c2k2n/totalcells+"% ("+c2k2n+")");
		Dialog.addMessage("1K0N: "+100*c1k0n/totalcells+"% ("+c1k0n+")");
		Dialog.addMessage("Other: "+100*other/totalcells+"% ("+other+")");
		Dialog.addMessage("Filtered cells:");
		Dialog.addMessage("Not a cell (0K0N): "+notcell);
		Dialog.addMessage("Touching edge: "+touchedge);
		Dialog.addMessage("Branched skeleton: "+branched);
	Dialog.show();
}

//Macro - Save analysis
//=============================================================================================================================
//Requires data in data table and all images from the analysis open and not renamed

macro "Save Analysis Action Tool -C115F11eeCaaaF5164CeeeF37a7" {
	run("Input/Output...", "jpeg=90 gif=-1 file=.txt copy_column copy_row save_column save_row");
	setBatchMode(true);
	if (nResults==0) {
		exit("Error: No data found in the results table");
	}
	datestring=""+getResult("Year", 1)+"_"+getResult("Month", 1)+"_"+getResult("Day", 1)+"_"+getResult("Hour", 1)+"_"+getResult("Minute", 1)+"_"+getResult("Second", 1);
	cnimages=0;
	for (i=1; i<nResults; i++) {
		if (isOpen("Analysed_"+datestring+"_"+getResult("Image No", i))==false || isOpen("AnalysedMasks_"+datestring+"_"+getResult("Image No", i))==false) {
			exit("Error: Some images from the analysis have been renamed or are not still open");
		}
		selectImage("Analysed_"+datestring+"_"+getResult("Image No", i));
		run("Make Composite", "display=Composite");
		selectImage("AnalysedMasks_"+datestring+"_"+getResult("Image No", i));
		run("Make Composite", "display=Composite");
		cnimages=maxOf(cnimages, getResult("Image No", i));
	}

	path=getDirectory("");
	File.makeDirectory(path+File.separator()+datestring);
	selectWindow("Results");
	run("Text...", "save="+path+File.separator()+datestring+File.separator()+"Results_"+datestring+".txt");
	experimentid=File.open(path+File.separator()+datestring+".txt");
		print(experimentid, "Analysis performed at: "+datestring);
		print(experimentid, cnimages+1+" images analysed.");
	File.close(experimentid);
	settings=File.open(path+File.separator()+datestring+File.separator()+"AnalysisSettings_"+datestring+".txt");
		print(settings, "slice_phase="+slice_phase);
		print(settings, "slice_dapi="+slice_dapi);
		print(settings, "slice_pi="+slice_pi);
		print(settings, "image_scale="+image_scale);
		print(settings, "distort_scale="+distort_scale);
		print(settings, "distort_xori="+distort_xori);
		print(settings, "distort_yori="+distort_yori);
		print(settings, "dapi_nucleus="+dapi_nucleus);
		print(settings, "dapi_kinetoplast="+dapi_kinetoplast);
		print(settings, "dapi_background="+dapi_background);
		print(settings, "pi_nucleus="+pi_nucleus);
		print(settings, "pi_kinetoplast="+pi_kinetoplast);
		print(settings, "pi_background="+pi_background);

		print(settings, "mcattolerance="+mcattolerance);
		print(settings, "mcatmaxdist="+mcatmaxdist);

		print(settings, "mkntolerance="+mkntolerance);

		print(settings, "cdsubtractbackground="+cdsubtractbackground);
		print(settings, "cdrollingballradius="+cdrollingballradius);

		print(settings, "displaymasks="+displaymasks);
		print(settings, "rollingballradius="+rollingballradius);
		print(settings, "cellautomanual="+cellautomanual);
		print(settings, "cellthresholdvalue="+cellthresholdvalue);
		print(settings, "cellautothreshtype="+cellautothreshtype);
		print(settings, "cellthresholdvalue="+cellthresholdvalue);
		print(settings, "cellminarea="+cellminarea);
		print(settings, "trimbranches="+trimbranches);
		print(settings, "minbranch="+minbranch);
		print(settings, "threshtype="+threshtype);
		print(settings, "autothreshtype="+autothreshtype);
		print(settings, "nthresholdvalue="+nthresholdvalue);
		print(settings, "kthresholdvalue="+kthresholdvalue);
		print(settings, "nucminarea="+nucminarea);
		print(settings, "kinminarea="+kinminarea);
	File.close(settings);

	for (i=0; i<=cnimages; i++) {
		selectImage("Analysed_"+datestring+"_"+i);
		//run("Make Composite", "display=Composite");
		//wait(100);
		slicestring="";
		for (j=0; j<nSlices(); j++) {
			if (j==slice_phase-1 || j==slice_dapi-1 || j==slice_pi-1) {
				slicestring=slicestring+"1";
			} else {
				slicestring=slicestring+"0";
			}
		}
		Stack.setActiveChannels(slicestring);
		setSlice(slice_phase);
		run("Grays");
		run("Enhance Contrast", "saturated=0.0");
		setSlice(slice_dapi);
		run("Orange Hot");
		run("Enhance Contrast", "saturated=0.0");
		setSlice(slice_pi);
		run("Cyan Hot");
		run("Enhance Contrast", "saturated=0.0");
		saveAs("ZIP", path+File.separator()+datestring+File.separator()+"Analysed_"+datestring+"_"+i+".zip");
		saveAs("Jpeg", path+File.separator()+datestring+File.separator()+"Analysed_"+datestring+"_"+i+".jpg");
		rename("Analysed_"+datestring+"_"+i);
		selectImage("AnalysedMasks_"+datestring+"_"+i);
		//run("Make Composite", "display=Composite");
		//wait(100);
		slicestring="11110";
		Stack.setActiveChannels(slicestring);
		for (j=0; j<nSlices(); j++) {
			setSlice(j+1);
			if (j==0) {
				run("Grays");
				run("Invert LUT");
				setMinAndMax(-128, 255);
			} else if (j==1) {
				run("Cyan Hot");
				setMinAndMax(0, 512);
			} else if (j==2) {
				run("Orange Hot");
				setMinAndMax(0, 512);
			} else if (j==3) {
				run("Grays");
				setMinAndMax(0, 16);
			} else if (j==4) {
				setMinAndMax(0, 255);
				run("Grays");
			}
		}
		saveAs("ZIP", path+File.separator()+datestring+File.separator()+"AnalysedMasks_"+datestring+"_"+i+".zip");
		saveAs("PNG", path+File.separator()+datestring+File.separator()+"AnalysedMasks_"+datestring+"_"+i+".png");
		rename("AnalysedMasks_"+datestring+"_"+i);
	}

	Dialog.create("DAPI Count Summary");
		Dialog.addMessage((cnimages+1)*2+" images saved.");
	Dialog.show();
}