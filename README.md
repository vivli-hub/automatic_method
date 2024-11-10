
# Method
An automatic method for the identification of cycles in Covid-19 time series data

# Authors
Miaotian Li, Ciprian Doru Giurcaneanu, Jiamou Liu

## Data
- [WHO_data.csv](#WHO_data)
- country_label.xlsx
- Latitude_label.xlsx

### WHO Data
In our experiments, we use the daily counts of COVID-19 cases that are publicly available on the website of WHO: [Data](https://covid19.who.int/data).

## Code
- [read_countries_WHO.m](#read_countries_WHO)
- [run_filt_WHO.m](#run_filt_WHO)
- [run_filt_fullsearch.m](#run_filt_fullsearch)

### read_countries_WHO
Data can be downloaded from: [Data](https://covid19.who.int/WHO-COVID-19-global-data.csv). This file contains the data that we have used in our experiments. Function read_countries_WHO.m reads the data from WHO_data.csv and creates a Mat file for each country in the folder WHO/WHO_Data. We have already created these files for all the countries; there are 236 countries in our experiments. We have also done the following data processing steps: (i) converted all NA and negative values to 0; (ii) Padded each time series with zeros.

### run_filt_WHO
Once you have these Mat files, you can run the main function run_filt_WHO.m.
For example, you can call the function like this:
run_filt_WHO(wreg, flag_transf, flag_plot, sp, ftype, fdomain, window)

**Input**

| **Variable**   | **Description**                                                                                         |
|----------------|---------------------------------------------------------------------------------------------------------|
| **`%wreg`**    | WHO region for which the data is analysed                                                               |
| **`%flag_transf`** | Flag for data transform:<br> 0 = No transform applied to the data<br> 1 = Transform from GP Nason, Scientific Reports, 2020<br> 2 = Transform usually applied to daily stock market indices at closing time<br> 3 = First order differences |
| **`%flag_plot`**   | 0 = No plots<br> 1 = Plots for each country                                                        |
| **`%sp`**      | 'day' = Daily data<br> 'week' = weekly data                                                             |
| **`%ftype`**   | 'stop' = Stop-band filter<br> 'high' = High-pass filter                                                |
| **`%fdomain`** | 'time' = Filtering in time domain<br> 'freq' = Filtering in frequency domain                           |
| **`%window`**  | 'hann' = Hanning window<br> 'halfhann' = Half Hanning window<br> 'rect' = Rectangular window           |

**Output**

run_filt_WHO produces two different files. Both of them have the same name, which is given by fina = strcat('Results_', wreg, '_tr', num2str(flag_transf), '_', sp, '_', ftype, '_', fdomain, ‘_’, window). However, one file has the extension .txt and the other file has the extension .mat.

**`%'Results_', wreg, '_tr', num2str(flag_transf), '_', sp, '_', ftype, '_', fdomain, ‘_’, window.txt`**

The .txt files record the region each country belongs to, along with the detection results for the fundamental period and harmonics. Note that the IT criterion used for detection is indicated in the results. If the result is marked with 'B', it means that BIC was able to detect the corresponding fundamental period or harmonics; if marked with 'S', it indicates that SC was employed in detection.

**`%'Results_', wreg, '_tr', num2str(flag_transf), '_', sp, '_', ftype, '_', fdomain, ‘_’, window.mat`**

| **Variable**   | **Description**                                                                                         |
|----------------|---------------------------------------------------------------------------------------------------------|
| **`seven`**    | `T0./(1:5)`                                                                                             |
| **`BIC_o`** or **`SC_o`** | The number of CTA's within a specific WHO region where fundamental period and harmonics can be detected (before filtering) by BIC and SC, respectively   |
| **`BIC_f`** or **`SC_f`** | The number of CTA's within a specific WHO Region where fundamental period and harmonics can be detected (after filtering) by BIC and SC, respectively |
| **`noc`**      | The total number of CTA’s in a specific WHO region                                                          |
| **`t_before`** | A matrix of size `n*5`, where `n` represents the number of CTA's in a particular region that can be detected with at least one periodogram or harmonics before applying the filter. The 5 columns correspond to the detection results from `T0` to `T0/5`, where the numbers represent the periods associated with frequencies detected by BIC or SC. Note: `n` does not represent the number of countries within a particular region, as some countries cannot be detected with any periodogram or harmonics. |
| **`t_after`**  | Similar to `t_before`, but stores the data conclusions after applying the filter. Only countries that were detected with a fundamental periodogram before applying the filter will have the filter applied. |

### run_filt_fullsearch

Once you have these Mat files, you can run the main function run_filt_fullsearch.m. You need to give as input a flag which indicates how the data are transformed (have a look at the comments in the function to understand how to select this flag)

**Input**

| **Variable**   | **Description**                                                                                         |
|----------------|---------------------------------------------------------------------------------------------------------|
| **`%wreg`**    | WHO region for which the data is analyzed                                                               |
| **`%flag_transf`** | Flag for data transformation:<br> 0 = No transformation applied to the data<br> 1 = Transformation from GP Nason, Scientific Reports, 2020<br> 2 = Transformation usually applied to daily stock market indices at closing time<br> 3 = First order differences (computed in time domain) |
| **`%sp`**      | 'day' = Daily data<br> 'week' = Iekly data                                                             |
| **`%ftype`**   | 'stop' = Stop band filter<br> 'high' = High pass filter                                                |
| **`%fdomain`** | 'time' = Filtering in time domain<br> 'freq' = Filtering in frequency domain                           |


**Output**

The output generated by this function is very similar to that of run_filt_WHO.m, but the .mat file only includes seven, t_before, and t_after. Both are algorithms for detecting the fundamental period and harmonics, with slight differences in their processes. The specific distinctions are explained in the PDF document.

## Data management
- [mean_sd.m](#mean_sd)
- [data_tidy.Rmd](#data_tidy)
- [generate_plot.Rmd](#generate_plot)
- [patition_comparison.Rmd](#partition_index)

### mean_sd
The mean and standard deviation of the estimation errors (before and after filtering) are computed for each WHO region. Calculate the mean and standard deviation of the bias caused by different WHO Regions, using different filters and domains.

**Input**
| **Variable**   | **Description**                                                                                         |
|----------------|---------------------------------------------------------------------------------------------------------|
| **`%wreg`**    | WHO region for which the data is analyzed                                                               |
| **`%flag_transf`** | Flag for data transformation:<br> 0 = No transformation is applied to the data<br> 1 = Transformation from GP Nason, Scientific Reports, 2020<br> 2 = Transformation usually applied to daily stock market indices at closing time<br> 3 = First order differences (computed in time domain) |
| **`%sp`**      | 'day' = Daily data<br> 'week' = Iekly data                                                             |
| **`%ftype`**   | 'stop' = Stop band filter<br> 'high' = High pass filter                                                |
| **`%fdomain`** | 'time' = Filtering in time domain<br> 'freq' = Filtering in frequency domain                           |

**Output**

**`%results_mean_before`** 

A 5x1 matrix, with the first column through the fifth column corresponding to T0, ..., T0/5. Compute the mean of the estimated error before filtering.

**`%results_sd_before`** 

A 5x1 matrix, with the first column through the fifth column corresponding to T0, ..., T0/5. Compute the standard deviation of the estimated error before filtering.

**`%results_mean_after`** 

A 5x1 matrix, with the first column through the fifth column corresponding to T0, ..., T0/5. Compute the mean of the estimated error after filtering.

**`%results_sd_after`** 

A 5x1 matrix, with the first column through the fifth column corresponding to T0, ..., T0/5. Compute the standard deviation of the estimated error after filtering.

### data_tidy
We will have daily data and weekly data, and the fundamental period and harmonics detected for each type of data are different. Therefore, we will have four functions in this Rmd file. 

`process_lines_daily(lines)`：Extract the data related to 'before filtering' from the TXT file and convert it into wide format. If this CTA corresponds to the relevant fundamental period or harmonics, mark it as 1; otherwise, mark it as 0. This function applies only to daily data. In this function, the input is coming from TXT files.

`process_lines_weekly(lines)`: This function is as same as the function `process_lines_daily`. The only difference is that this function is for weekly data.

`process_lines_daily_after(lines)` : This function is as same as the function `process_lines_daily`. The only difference is that this function analyzes the results after filtering.

`process_lines_week_after(lines)`: This function is as same as the function `process_lines_daily_after`. The only difference is that this function is for weekly data.

**Input**

**`'Results_', wreg, '_tr', num2str(flag_transf), '_', sp, '_', ftype, '_', fdomain, ‘_’, window.txt`**

**Output**

**`before_filter_daily.csv`**: 

| Variable Name | Levels | Explanation |
|---------------|--------|-------------|
| country       | --     | The name of CTA’s |
| WHO_Region    | AFRO; AMRO; EMRO; EURO; SEARO; WPRO | The name of six WHO regions |
| T0            | 1 (can be detected); 0 (cannot be detected) | 7-day cycle |
| T0/2          | 1 (can be detected); 0 (cannot be detected) | 7/2-day cycle |
| T0/3          | 1 (can be detected); 0 (cannot be detected) | 7/3-day cycle |
| T0/4          | 1 (can be detected); 0 (cannot be detected) | 7/4-day cycle |
| T0/5          | 1 (can be detected); 0 (cannot be detected) | 7/5-day cycle |
| criterion     | BIC; SC | Two IT criteria |
| transform     | 0; 1; 2; 3 | Three different kinds of transforms |
| filter        | none   | Don’t use any filter |
| domain        | none   | -- |



**`before_filter_week.csv`**: 

It is exactly the same as before_filter_daily.csv, with the only difference being that it focuses on confirming the annual cycle (52-week cycle). In addition, there is a new column called `latitude_region`.

**`after_filter_daily.csv`**: 

It is exactly the same as before_filter_daily.csv. But the levels of filter and domain are different. There are two kinds of filters: high-pass and stop-band, and also two domains: time and frequency.

**`after_filter_daily.csv`**: 

It is exactly the same as before_filter_week.csv.

### generate_plot
Represent the results of the above experiments. Use ggplot to aggregate and summarize the daily data and weekly data by WHO Region group; summarize the weekly data by Latitude group.

| **Type**   | **Files**                                                                                             |
|------------|-------------------------------------------------------------------------------------------------------|
| **Input**  | `before_filter_daily.csv`<br>`after_filter_daily.csv`<br>`before_filter_week.csv`<br>`after_filter_week.csv` |
| **Output** | `figure_2.pdf`<br>`figure_3.pdf`<br>`figure_4.pdf`<br>`figure_5.pdf`<br>`figure_6.pdf`<br>`figure_7.pdf`<br>`figure_B1.pdf`<br>`figure_B2.pdf`<br>`figure_B3.pdf`<br>`figure_B4.pdf` |


### partition_index
We use the ‘partitionComparison’ package in R environment for computing the following similarity indexes for the two partitions: Rand, Jaccard and Fowlkes-Mallows ). 

**Input**

**`country_label.xlsx`**

**`Latitude_label.xlsx`**

**Output**

Rand, Jaccard and Fowlkes-Mallows Index


The full report can be found at [Link]().
