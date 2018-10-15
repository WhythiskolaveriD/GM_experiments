# README #

# Contents
This repo contains all the Global Mass experiments I have done so far. It now contains the followings:

* Experiment 1a: synthetic data
* Experiment 1b: Update global GIA from GPS
* Experiment 2: Update global Steric from SSH and Mass
* Experiment 3: Update North America GIA 



# Before running
Before running you should install the following on your pc or server.

* Latest version of R

* Install the latest version of [Rstudio](https://www.rstudio.com/) (One of the best development environmet software for R)

* Latest version of the [INLA](http://www.r-inla.org/) package ([installation instructions here](http://www.r-inla.org/download))

* Source the code scripts in the [gmrcode] folder in this repo

* Install the following R packages (and any other packages required by the following or as you need)

     * sp
     
     * rgdal
    
     * GEOmap
    
     * ncdf4
    
     * knitr
  
## Migrating R packages
  
  You can upgrade your R version whenever you want, but you cannot control this on the server. Meanwhile, the R packages you have installed are usually stored under the folder related to the version name of R. On you own PC, you can change this to an permanent address but this might not be possible on a server. The following summarizes three methods for migrating installed R pacages.
  
### Method (1)

Suppose your packages from the old version are all installed under  the folder ``~/R/x86_64-redhat-linux-gnu-library/3.4``; you new version will have packages all installed under ``~/R/x86_64-redhat-linux-gnu-library/3.5``.

You can simply copy all the packages you have already installed under the old version to the new version folder and then run the following line in R ``update.packages(checkBuilt=TRUE)``

This should work for most packages. If it does not work, just re-install the packages.

### Method (2)

Similar to method(1), you can create a central library folder that keep the packages for all versions and add upgrade the central folder. To add the central folder to your R search path, use ``.libPaths(“path2central”)``.

You can add this search path permanently by adding ``R_LIBS=path2central`` to your .Rprofile on linux, or change the environment variables on windows.

### Method (3)

If still not working, the brutal way is to delete everything and re-install them. Before deleting, you can keep a copy of the names of the installed packages by ``my.packages <- installed.packages()[,1]`` and save this vector in a txt file. 

Then open the new version of R, load the “my.packages” from the text file and  install all these packages by ``install.packages(my.packages, dep = TRUE)``.

## Upgrading INLA on server

INLA requires certain built of c++ libraries which may not be avaible on the server you are running and you either don't have the right to rebuilt these things from source. The INLA team provived patches for different types of server. You can download these patches and replace the corresponding folder in the R INLA package.

For example, you can solve this by the following steps.

(1)	Install the latest version of INLA as shown here http://www.r-inla.org/download

(2)	Then download and unpack the files from here http://inla.r-inla-download.org/Linux-builds/CentOS%20Linux-7%20(Core)/Version_18.08.09/ 

(3)	Replace the folder “64bit” under where you installed INLA by the folder in (2) . For me on Atlantis, it’s under  “~/R/x86_64-redhat-linux-gnu-library/3.5/INLA/bin/linux/64bit/”



# How to use

All the experiments are documented in the format of .Rmd, which is an R markdown file of both codes and comments and can be converted into html/pdf format for presentation. 

* For initial reading or final presentation, please use the rendered html file. 

* Current .Rmd file only include skeleton of the experiment. You can play with the Rmd file to add more plots or features and render an html file for presentation.

* For more details of code and preambles, please read the R script in each section.

    * The .Rmd file are set to be skipping all the code chunk for fast rendering the html file.
    
    * Some of the codes and results are hidden in the html presentation but can be found in the .Rmd file.
    
    * You can open the .Rmd file in R studio and run the code chunks interactively
    
    * Some of the longer code chunks are written in an R script file; the path is given explicitly in the .Rmd file
 
 * For running the entire experiment, please run the R script on the servers. The next setion explains how to run the R scripts on a server. 

# Running on servers

To run an experiment on servers, you should use the R script corresponding to the experiment, rather than the Rmd file. Running a R script on a server can be done in two ways.

#### 1 Running in an R session
You can run an R session on a server by typing `R` and return in an command line environment. Then you can source the script by

```r
source("mypath/myscript.R")
```

#### 2 Run Rscript in a command line

You can also run the script from the command line. This is useful when the job need to to incorporate in an command line script, e.g. submitting a job to the clusters.

The simplest one one is
```
yourname@server $ R CMD BATCH --no-save --no-restore --slave mypath/myscript.R &
```

Or you can add some arguements and output the log
```
yourname@server $ R CMD BATCH --no-save --no-restore --slave '--args input="abc" ' mypath/myscript.R output.txt &
```

### Load all required modules 

On Atlantis you usually will not need this except for geospatial related operations(gdal etc). On Bluecrystal you need to specify all related modules including R and C++ libraries. 

An example of submitting jobs to Bluecrystal can be found in the ``Exp2.Rmd``.


