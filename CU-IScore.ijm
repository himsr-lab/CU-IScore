/*  Copyright 2021 Regents of the University of Colorado
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <https://www.gnu.org/licenses/>.
 *
 *  Author:       Christian Rickert <christian.rickert@cuanschutz.edu>
 *  Group:        Human Immune Monitoring Shared Resource (HIMSR)
 *                University of Colorado, Anschutz Medical Campus
 *
 *  Title:        CU-IScore
 *  Summary:      ImageJ2 macro for the IHC scoring of multi-channel images
 *
 *  DOI:          https://doi.org/10.5281/zenodo.4599591
 *  URL:          https://github.com/christianrickert/CU-IScore/
 *
 *  Description:
 *
 *  CU-IScore uses the Fiji image processing package exclusively – there is no need
 *  to install additional plugins: You can simply edit or run this code by loading and
 *  executing it with the Macro plugin.
 *  The scoring of immunofluorescence images is performed analog to the approach
 *  described in the publication for the IHC Profiler:
 *
 *  Varghese F., Bukhari A.B., Malhotra R., De A. (2014) IHC Profiler: An Open Source
 *  Plugin for the Quantitative Evaluation and Automated Scoring of Immunohistochemistry
 *  Images of Human Tissue Samples. PLoS ONE 9(5): e96801.
 *  https://doi.org/10.1371/journal.pone.0096801
 *
 *  A reference implementation (ImageJ plugin) by Dmitry Brant can be found here:
 *
 *  https://github.com/dbrant/ihc-profiler
 *
 *  CU-IScore calculates the I-Score of images or slices of image stacks by separating
 *  the corresponding grayscale histograms into four intervals with weighted groups.
 *  Our scoring method returns higher score values for images with higher average pixel
 *  values by default – which is an applicable for the scoring of immunofluorescence data.
 *  However, the IHC Profiler scoring method for use with DAB/hematoxylin color spectra
 *  images can be activated in the code for images with separate DAB and hematoxylin
 *  channels. The grayscale histogram groups are weighted given the following table:
 *
 *  I-Score groups        |  Weight  |   Interval   |  Weight  |  IHC Profiler groups
 *  ----------------------|--------------------------------------------------------------
 *  Background        0   |    0     |  [min-24%]   |    3     |  High positive  3+
 *  Low intensity     1+  |    1     |  [25%-49%]   |    2     |  Positive       2+
 *  Medium intensity  2+  |    2     |  [50%-74%]   |    1     |  Low positive   1+
 *  High intensity    3+  |    3     |  [75%-max]   |    0     |  Negative       0
 *  Outliers          o   |    0     |  [outside]   |    0     |  -
 *
 *  The calculation key takes into account the percentages of pixels in each of the
 *  weighted groups described by these two equations:
 *
 *  n(pixel) = n('3+') + n('2+') + n('1+')
 *  Score    = (100 / n(pixel)) * (3 * n('3+') + 2 * ('2+') + ('1+'))
 *
 *  When comparing scores between datasets, the same pixel value range must be applied:
 *  The intervals can be determined automatically by CU-IScore, which will perform a
 *  pre-run to determine the pixel value ranges on a per-slice basis. Alternatively,
 *  users can enter fixed pixel value ranges to exclude outliers from the scoring.
 */

/*
 *  Imports
 */

run("Bio-Formats Macro Extensions");

/*
 * Variables
 */

// file matching pattern
suffixes = newArray(".tif", ".tiff");  // match file name endings (not only extensions)

// scoring settings
globalRanges = true;  // determine global slice histograms
classicMode = false;  // report IHC-Score, instead of I-Score

// advanced user settings
batchMode = true;  // speed up processing by limiting visual output
fixedRanges = false;  // use user-specified ranges for outlier removal
var maxima = newArray(0);  // maximum pixel values (global variable)
var minima = newArray(0);  // minimum pixel values (global variable)
versionString = "CU-IScore v0.9 (2021-04-01)";

/*
 *  Start
 */

print("\\Clear");
requires("1.52t");  // minimum ImageJ version
tableName = "CU-IScore";
Table.create(tableName);  // creates and resets a table
file = File.openDialog("Select the first TIFF of your dataset");
if ( batchMode )
  setBatchMode(true);  // enter batch mode and don't display newly opened images
processFolder(file, suffixes, tableName);
print("\n*** Saving table to file ***");
path = File.getDirectory(file);
csvFile = path + File.separator + tableName + ".csv";
waitForFileDeletion(csvFile);
Table.save(csvFile);
printDateTimeStamp();
if ( batchMode )
  setBatchMode(false);  // exit batch mode and dispose of all but the most recent image

/*
 *  Loop
 */

// Function to score files with matching suffixes from a folder
function processFolder(file, suffixes, tableName)
{
  folder = File.getParent(file);
  objects = getFileList(folder);  // files and folders
  objectCount = objects.length;

  // check files for global extrema
  if ( !fixedRanges && globalRanges )
  {
    initializeRun();
    checkFile(file);  // file selected by user

    for (i = 0; i < objectCount; ++i)  // all other files in folder
    {
      objectPath = folder + File.separator + objects[i];
      if(endsWithEither(objectPath, suffixes) && objectPath != file)
        checkFile(objectPath);
    }

  }

  // calcucate the I-Scores
  rowIndex = 0;  // table row index
  scoreFile(file, tableName, rowIndex);

  for (i = 0; i < objectCount; ++i)
  {
    objectPath = folder + File.separator + objects[i];
    if(endsWithEither(objectPath, suffixes) && objectPath != file)
    {
      rowIndex +=1;
      scoreFile(objectPath, tableName, rowIndex);
    }
  }

}

// Function to check a single file for its intervals
function checkFile(file)
{
  // For the global range setting, we check all files with their respective
  // channels for the minimum and maximum values. The values are stored in
  // two arrays, so that the value at the first indices corresponds to the
  // minimum and maximum values of the first channel, respectively.
  print("\n*** Checking file ***");
  print("\t" + file);

  // read image file
  fileName = File.getName(file);
  fileSlices = readImage(file);

  // get image metadata
  slices = nSlices();  // channel count

  // match array sizes with slice number, works with global variables only
  maxima = extendArray(maxima, slices, -1e30);
  minima = extendArray(minima, slices, 1e30);

  // update global max and min values per slice
  for ( slice = 0; slice < slices; ++slice )
  {
    setSlice(slice + 1);
    getRawStatistics(pixels, mean, min, max);
    if ( max > maxima[slice] )
      maxima[slice] = max;
    if ( min < minima[slice] )
      minima[slice] = min;
  }

  // close datasets and images
  Ext.close();
  close("*");
}

// Function to score a single file
function scoreFile(file, tableName, rowIndex)
{
  // The score calculation is centered around the getHistogram function, which
  // groups pixel values into bins and returns the bin counts. Most of the code
  // deals with input checks and optimizes the histogram parameters for use
  // with each of the channels. Monochromatic images are not scored.

  // prepare next run
  initializeRun();
  print("\n*** Scoring file ***");
  print("\t" + file);
  histoBins = 0;  // number of bins for histogram clustering
  scoreAverage = 0.0;  // stack score average
  scoreCount = 0.0;  // stack score count
  scoreSum = 0.0;  // stack score sum

  // read image file
  filePath = File.getDirectory(file);
  fileName = File.getName(file);
  fileSlices = readImage(file);
  fileTitle = getTitle();

  // get image metadata
  bits = bitDepth();
  print("\tBit depth: " + bits);
  getRawStatistics(pixels);
  slices = nSlices();

  // match array sizes with slice number
  maxima = extendArray(maxima, slices, 0.0);
  minima = extendArray(minima, slices, 0.0);
  print("\tFixed ranges: " + fixedRanges);
  print("\tGlobal ranges: " + globalRanges);

  // create table entries with incremental indexing
  Table.set("File", rowIndex, fileName, tableName);  // first column
  if ( classicMode )
  {
    Table.set("IHC-Score", rowIndex, "-", tableName);  // second column
    print("\tScoring method: IHC-Score");
  }
  else
  {
    Table.set("I-Score", rowIndex, "-", tableName);
    print("\tScoring method: I-Score");
  }

  // score slices in image stack
  for (slice = 0; slice < slices; ++slice)
  {
    setSlice(slice + 1);  // current channel
    sliceLabel = getMetadata("Label");
    print("\n\tChannel: " + sliceLabel);

    // check slice range
    if ( fixedRanges || globalRanges )
      getRawStatistics(pixels);  // get pixel count only
    else
    {
      getRawStatistics(pixels, mean, min, max);  // extend range from slice
      if ( max > maxima[slice] )
        maxima[slice] = max;
      if ( min < minima[slice] )
        minima[slice] = min;
    }
    valueRange = maxima[slice] - minima[slice];  // pixel value range
    print("\t\tValue range:     " + valueRange);

    // histogram-based score calculation
    if ( valueRange > 0 ) // image slice with dynamic range
    {
      groupCounts = initializeArray(5, 0);  // pixel counts of groups
      groupLimits = initializeArray(3, 0);  // right bounds of groups
      sliceScore = 0;  // channel score
      scorePixels = 0;  // pixel count for channel score, n(pixel)
      if ( bits == 32 )
      {
        histoBins = Math.round(valueRange / 0.001);  // bin number needs to be greater than zero
        if ( histoBins % 2 != 0 )  // avoid binning artifacts with odd number of bins
          histoBins += 1;
        if ( histoBins < 1024 )  // ensure minimum histogram resolution
          histoBins = 1024;
      }
      else if ( bits == 16 )
        histoBins = 65534;  // histoBounds calculation fails with 65536
      else // 24-bit, 8-bit,
        histoBins = 256;  // must be 256
      print("\t\tHistogram bins:  " + histoBins);
      histoBounds = initializeArray(histoBins, 0);  // right bounds of histogram bins
      histoPixels = initializeArray(histoBins, 0);  // pixel counts of histogram bins

      // create image histogram
      getHistogram(histoBounds, histoPixels, histoBins);  // count pixel intensities by bins

      // set up interval limits
      groupLimits[0] = minima[slice] + 0.25 * valueRange;  // background (0) or high positive (3+)
      groupLimits[1] = minima[slice] + 0.50 * valueRange;  // low intensity (1+) or positive (2+)
      groupLimits[2] = minima[slice] + 0.75 * valueRange;  // medium intensity (2+) or low positive (1+)

      // assign pixels in histogram bins to intervals
      groupCounts = countPixels(slice, groupLimits, histoBounds, histoPixels, histoBins);

      // calculate pixel count for channel score
      if ( classicMode )
        scorePixels += groupCounts[0];
      scorePixels += groupCounts[1];
      scorePixels += groupCounts[2];
      if ( !classicMode )
        scorePixels += groupCounts[3];

      // calculate score
      if ( scorePixels > 0 )  // scoring successful
        sliceScore = calculateScore(groupCounts, scorePixels);

      // print values to log file
      print("\t\tPixels (total):  " + (groupCounts[0] + groupCounts[1] + groupCounts[2]
                                    +  groupCounts[3] + groupCounts[4])
                                    + " (" + pixels + ")");
      reportScore(slice, minima, maxima, groupLimits, groupCounts, sliceScore, classicMode);

      // add channel value to table
      Table.set(sliceLabel, rowIndex, sliceScore, tableName);

      // keep track of scores
      scoreCount += 1.0;
      scoreSum += sliceScore;
    }
    else
      print("\t\t(No scoring)");
  }

  if ( scoreCount > 0 )  // all slices monochromatic
  {
    // report stack score
    scoreAverage = scoreSum / scoreCount;
    if ( classicMode )
    {
      print("\n\tIHC-Score: " + scoreAverage + " (channel average)");
      Table.set("IHC-Score", rowIndex, scoreAverage, tableName);  // overwrite existing value
    }
    else
    {
      print("\n\tI-Score: " + scoreAverage + " (channel average)");
      Table.set("I-Score", rowIndex, scoreAverage, tableName);
    }
  }
  // save and clear run, free memory
  finalizeRun(filePath, fileName);
  call("java.lang.System.gc");  // trigger JVM garbage collection
}

/*
 *  Functions
 */

// Function to calculate the score
function calculateScore(counts, pixels)
{
  // Keep in mind that the scores are calculated with relative pixel counts -
  // using percentage values of the relevant pixel count totals.
  output = 0;

  if ( classicMode )  // IHC-score
    output = (100.0 / pixels) * (3 * counts[0] + 2 * counts[1] + counts[2]);
  else  // I-Score
    output = (100.0 / pixels) * (3 * counts[3] + 2 * counts[2] + counts[1]);
  return output;
}

// Function to count pixels in groups
function countPixels(slice, limits, bounds, pixels, bins)
{
  // For the scoring, we need to concatenate the histogram group counts according
  // to the definition of the scoring method chosen. For performance reasons, we
  // avoid unnecessary indexing of the limits array.
  output = initializeArray(5,0);
  bckMax = limits[0];
  lowMax = limits[1];
  medMax = limits[2];

  for ( b = 0; b < bins; ++b )
  {
    if ( bounds[b] < minima[slice] || bounds[b] > maxima[slice] )  // outside specified range
      output[4] += pixels[b];
    else  // within specified range
    {
      if ( bounds[b] < bckMax )
        output[0] += pixels[b];
      else if ( bounds[b] < lowMax )
        output[1] += pixels[b];
      else if ( bounds[b] < medMax )
        output[2] += pixels[b];
      else
        output[3] += pixels[b];
    }
  }

  return output;
}

// Function to test whether a string ends with suffixes from a list
function endsWithEither(string, suffixes)
{
  found = false;

  for (i = 0; i < suffixes.length; ++i)
  {
    if (endsWith(string, suffixes[i]))
      found = true;
  }

  return found;
}

// Function extends and returns an array with default values
function extendArray(array, size, value)
{
  while ( array.length < size )
  {
    array = Array.concat(array, value);
  }

  return array;
}

// Function to finalize a completed segmentation run
function finalizeRun(path, name)
{
  print("\n*** Saving result to file ***");

  label = File.getNameWithoutExtension(name);
  directory = path + File.separator + label;  // create subfolder for result files
  File.makeDirectory(directory);
  result = directory + File.separator + "CU-IScore";  // full file path without extension
  Ext.close();  // close active Bio-Formats dataset
  close("*");  // closes all image windows
  txtFile = result + ".txt";
  print("\tWriting: " + txtFile);
  waitForFileDeletion(txtFile);
  printDateTimeStamp();
  saveLogFile(txtFile);
}

// Function initializes and returns an array with default values
function initializeArray(counts, value)
{
  return Array.fill(newArray(counts), value);
}

// Function to initialize a new segmentation run
function initializeRun()
{
  print("\\Clear");  // clear Log window
  printDateTimeStamp();
  print("ImageJ2 v" + IJ.getFullVersion);
  print(versionString);
  run("Close All");  // close all image windows
}

// Function to print date and time
function printDateTimeStamp()
{
  getDateAndTime(year, month, dayOfWeek, dayOfMonth, hour, minute, second, msec);

  print("\n[" + year + "-" + (month + 1)  + "-" + dayOfMonth + ", "
              + hour + ":" + minute + ":" + second + "]");
}

// Function to read an image file with Bio-Formats
function readImage(file)
{
  // Since we don't know all possible metadata labels for different imaging platforms,
  // we've decided to implement a generic pattern matching system that allows us to
  // simply add new labels to the list of known labels by extending the prefix list.
  // As counting of those metadata labels can be done with or without leading zeros,
  // we're using the same mechansim to try different combinations of labels with counts,
  // until we get a valid response from the Ext.getMetadataValue function.
  counting = "";  // metadata key count subsequential to key
  metadataPrefixes = newArray("ChannelName #", "Information|Image|Channel|Name #",
                              "Name #", "PageName #");
  metadataCounting = newArray("0", "01", "1");
  countingLength = metadataCounting.length;
  prefixesLength = metadataPrefixes.length;
  offset = 0;  // offset between key and slide numbering
  prefix = "";  // metadata key prefix preceeding count
  slice = 0;  // default "not found" return value from Ext.getMetadataValue()
  slices = newArray(0);  // labels
  sliceCount = 0;

  run("Bio-Formats Windowless Importer", "open=" + v2p(file));
  Ext.setId(file);  // initializes the given id (file)
  Ext.getImageCount(sliceCount);

  for (i = 0; slice == 0 && i < prefixesLength; ++i)  // try different metadata prefixes
  {

    for (j = 0; slice == 0 && j < countingLength; ++j)  // try different metadata counts
    {
      Ext.getMetadataValue(metadataPrefixes[i] + metadataCounting[j], slice);
      if ( slice != 0 )  // matching prefix/count pair
      {
        prefix = metadataPrefixes[i];
        counting = metadataCounting[j];
        if ( counting == "0" )
          offset = -1;
      }
    }

  }

  for (k = 1; k <= sliceCount; ++k)  // extract metadata for each slice
  {
    setSlice(k);
    if ( counting == "01" && k <= 9 )  // Polaris
      Ext.getMetadataValue(prefix + "0" + k, slice);
    else  // MIBI, Vectra, Zeiss
      Ext.getMetadataValue(prefix + (k + offset), slice);
    if ( slice == 0 )  // no compatible metadata label found
      slice = k;  // use number instead
    setMetadata("Label", slice);  // add label to slice
    slices = Array.concat(slices, toString(slice));
    print("\t" + k + ".) " + slice);
  }

  run("Maximize");  // window pane
  return slices;
}

// Function to report the score
function reportScore(slice, minima, maxima, limits, counts, score, classic)
{
  if ( classic )
  {
    print("\t\tHigh positive    ["  + minima[slice] + "-("  + limits[0] + ")]: " + counts[0]);
    print("\t\tPostive          ["  + limits[0] + "-("  + limits[1] + ")]: " + counts[1]);
    print("\t\tLow positive     ["  + limits[1] + "-("  + limits[2] + ")]: " + counts[2]);
    print("\t\tNegative         ["  + limits[2] + "-"   + maxima[slice] +  "]: " + counts[3]);
    print("\t\t-                [<" + minima[slice] + ", >" + maxima[slice] +  "]: " + counts[4]);
    print("\t\tIHC-Score:       " + score);
  }
  else
  {
    print("\t\tBackground       ["  + minima[slice] + "-("  + limits[0] + ")]: " + counts[0]);
    print("\t\tLow intensity    ["  + limits[0] + "-("  + limits[1] + ")]: " + counts[1]);
    print("\t\tNormal intensity ["  + limits[1] + "-("  + limits[2] + ")]: " + counts[2]);
    print("\t\tHigh intensity   ["  + limits[2] + "-"   + maxima[slice] +  "]: " + counts[3]);
    print("\t\tOutliers         [<" + minima[slice] + ", >" + maxima[slice] +  "]: " + counts[4]);
    print("\t\tI-Score:         " + score);
  }
}

// Function to save the Log window content
function saveLogFile(file)
{
  selectWindow("Log"); //select Log window
  saveAs("Text", file);
}

// Function to convert a value to a parameter string
function v2p(value)
{
  // This is a workaround for "undefined variable" errors when using
  // the address operator (&) in functions' parameter assignments.
  // By converting the numeric value to a string, we also avoid
  // string concatenation problems during parameter parsing.
  return "[" + toString(value) + "]";
}

// Function to wait for a certain file to be deleted
function waitForFileDeletion(file)
{
  // Network drives might be slow or files might be locked, i.e. opened in
  // another application and not accessible for deletion. We're therefore
  // notifying the user and waiting until the lock has been removed.
  notified = false;

  while ( File.exists(file) )
  {
    deleted = File.delete(file);
    if ( deleted != 1 && !notified )  // deleting failed, notify user
    {
      print("\tCan't delete file. Close other programs accessing the file.");
      notified = true;
    }
    wait(1000);  // ms
  }

}
