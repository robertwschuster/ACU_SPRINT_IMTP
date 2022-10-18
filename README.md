# ACU_SPRINT_IMTP
R code and shiny app for isometric mid-thigh pull analysis

For the IMTP analysis shiny app:
1.  Open the files you want to analyse
  a.  You can select multiple files containing multiple trials from multiple participants
2.  The plots below the table will allow you to identify poor trials
  a.  If a file contains multiple trials, these will be plotted in separate tabs
  b.  If an indicator of poor trial quality is detected, it will be printed below the plot
3.  Select the trials for which you want to export the performance variables by clicking the corresponding rows in the data table
4.  Enter the name of the CSV file you want to save these results to
5. Click 'Download selected rows' and navigate to where you want to save the results

For IMTP_analysis.R:
1.	Open the script in RStudio.
2.	Click “Source”. [CAVE: if “Run” is used the script won’t perform as intended.]
3.	A browser window to select the files to analyse will automatically open. [CAVE: this window is often in the background, so make sure to check your task bar.]
4.	Navigate to where your files are located, select the files you want to analyse and click “Open”.
5.	If certain criteria for questionable IMTP trial quality are met, you will be prompted to choose whether or not to continue with the analysis. In the console, type either “y” (for yes) or “n” (for no) and hit enter.
  a.	You may get multiple prompts for the same file/trial.
  b.	The file/trial in question will be plotted in the plot window so you can visually inspect its quality.
6.	Once all files are processed, you will be asked whether or not to export the results. Type “y” or “n” in the console and hit enter.
  a.	If you decide to export the results, you will be asked to enter the filename under which you want to save the results. Type the filename in the console and hit enter.
  b.	A csv file will be saved in the same directory the files you selected for analysis are located.
