//call("BIOP_LibInstaller.installLibrary", "BIOP"+File.separator+"BIOPLib.ijm");
function prepareTable(tableName) {
                updateResults();
                if(isOpen("Results")) { IJ.renameResults("Results","Temp"); updateResults();}
                if(isOpen(tableName)) { IJ.renameResults(tableName,"Results"); updateResults();}

}

function parseJACop(input){
open(input);
name=getTitle();
C1="C1-"+name;
C2="C2-"+name;
run("32-bit");
run("Split Channels");

selectWindow(C1);
setAutoThreshold("Moments dark stack");
run("NaN Background", "stack");
run("8-bit");

selectWindow(C2);
setAutoThreshold("Moments dark stack");
run("NaN Background", "stack");
run("8-bit");

run("JACoP ");

selectWindow(C1);
selectWindow(C2);
run("JACoP ", "imga="+C1+" imgb="+C2+" thra=110 thrb=124 pearson overlap mm");




run("Close All");

logdump=split(getInfo("log"),"\n");
for (i=0;i<logdump.length;i++){
	    if (startsWith(logdump[i],"Image A"))
	        imgA=substring(logdump[i],9,lengthOf(logdump[i]));
	    if (startsWith(logdump[i],"Image B"))
	        imgB=substring(logdump[i],9,lengthOf(logdump[i]));
	    if (startsWith(logdump[i],"Pearson's Coefficient"))
	        Pc = parseFloat(substring(logdump[i+1], 2, lengthOf(logdump[i+1])));


                if (startsWith(logdump[i], "Overlap Coefficient"))
                                Oc = parseFloat(substring(logdump[i+1], 2, lengthOf(logdump[i+1])));
                        

                if (startsWith(logdump[i], "k1=")) 
                                k1 = parseFloat(substring(logdump[i], 3, lengthOf(logdump[i])));
                     
         
                if (startsWith(logdump[i], "k2="))
                                k2 = parseFloat(substring(logdump[i], 3, lengthOf(logdump[i])));
                      
                 if (startsWith(logdump[i], "Using thresholds")) {
                        thrA = parseFloat(substring(logdump[i], indexOf(logdump[i], "=")+1, indexOf(logdump[i], "and")-1));
                        thrB = parseFloat(substring(logdump[i], lastIndexOf(logdump[i], "=")+1, lastIndexOf(logdump[i], ")")));
                        thrVals = true;
                        //OcThr = parseFloat(substring(logdump[i+3], 2, lengthOf(logdump[i+1])));
                       // K1Thr = parseFloat(substring(logdump[i+6], 3, lengthOf(logdump[i])));
                       // K2Thr = parseFloat(substring(logdump[i+7], 3, lengthOf(logdump[i])));
                }
                if (startsWith(logdump[i], "Manders' Coefficients (original):")) {
                        M1 = parseFloat(substring(logdump[i+1], 3, lengthOf(logdump[i+1])));
                        M2 = parseFloat(substring(logdump[i+2], 3, lengthOf(logdump[i+2])));

                }

                if (startsWith(logdump[i], "Manders' Coefficients (using threshold")) {
                        M1Thr = parseFloat(substring(logdump[i+1], 3, lengthOf(logdump[i+1])));
                        M2Thr = parseFloat(substring(logdump[i+2], 3, lengthOf(logdump[i+2])));

                } 
     
	}
	close("Log");
	

}

//input="/Users/";

var imgA;
var imgB;
var Pc;
var Oc;
var k1;
var k2;
var M1;
var M2;
var thrA;
var thrB;
var M1Thr;
var M2Thr;
input="H:/";
setBatchMode(true);
list=getFileList(input);
prepareTable("Coloc Results");
for (i=0;i<list.length;i++){
    parseJACop(input+list[i]);

                
                n=i;

                setResult("imageA",n,imgA);
                setResult("imageB",n,imgB);
                setResult("Pearson's", n, Pc);
                setResult("Overlap Coefficient (threshold)", n, Oc);
                setResult("k1 (threshold)", n, k1);
                setResult("k2 (threshold)", n, k2);
                //setResult("Overlap Coeffecient (threshold)",n,OcThr);
                //setResult("k1 (threshold)",n,k1Thr);
                //setResult("k2 (threshold)",n,k2Thr);
                setResult("M1 (no threshold)", n, M1);
                setResult("M2 (no threshold)", n, M2);
                setResult("ThrA", n, thrA);
                setResult("ThrB", n, thrB);
                setResult("M1", n, M1Thr);
                setResult("M2", n, M2Thr);
}
setBatchMode(false);


