
# Method
Automatic method for the identification of cycles in Covid-19 time-series data

# Authors
Miaotian Li, Ciprian Doru Giurcaneanu, Jiamou Liu

## Data
- [WHO_data.csv](#WHO_data)
- [country_label.xlsx](#partition_index)
- [Latitude_label.xlsx](#partition_index)

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
| **`%wreg`**    | WHO region for which the data is analyzed                                                               |
| **`%flag_transf`** | Flag for data transform:<br> 0 = No transform applied to the data<br> 1 = Transform from GP Nason, Scientific Reports, 2020<br> 2 = Transform usually applied to daily stock market indices at closing time<br> 3 = First order differences |
| **`%flag_plot`**   | 0 = No plots<br> 1 = Plots for each country                                                        |
| **`%sp`**      | 'day' = Daily data<br> 'week' = weekly data                                                             |
| **`%ftype`**   | 'stop' = Stop-band filter<br> 'high' = High-pass filter                                                |
| **`%fdomain`** | 'time' = Filtering in time domain<br> 'freq' = Filtering in frequency domain                           |
| **`%window`**  | 'hann' = Hanning window<br> 'halfhann' = Half Hanning window<br> 'rect' = Rectangular window           |

**Output**

run_filt_WHO produces two different files. Both of them have the same name, which is given by `fina = strcat(Results_', wreg, '_tr', num2str(flag_transf), '_', sp, '_', ftype, '_', fdomain, '_', window)`. However, one file has the extension .txt and the other file has the extension .mat.

**`%'Results_', wreg, '_tr', num2str(flag_transf), '_', sp, '_', ftype, '_', fdomain, '_', window.txt`**

The .txt files record the region each country belongs to, along with the detection results for the fundamental period and harmonics. Note that the IT criterion used for detection is indicated in the results. If the result is marked with 'B', it means that BIC was able to detect the corresponding fundamental period or harmonics; if marked with 'S', it indicates that SC was employed in detection.

**`%'Results_', wreg, '_tr', num2str(flag_transf), '_', sp, '_', ftype, '_', fdomain, ‘_’, window.mat`**

| **Variable**   | **Description**                                                                                         |
|----------------|---------------------------------------------------------------------------------------------------------|
| **`seven`**    | `T0./(1:5)`                                                                                             |
| **`BIC_o`** or **`SC_o`** | The number of CTA's within a specific WHO region where fundamental period and harmonics have been detected (before filtering) by BIC and SC, respectively   |
| **`BIC_f`** or **`SC_f`** | The number of CTA's within a specific WHO Region where fundamental period and harmonics have been be detected (after filtering) by BIC and SC, respectively |
| **`noc`**      | The total number of CTA’s in a specific WHO region                                                          |
| **`t_before`** | An array of size `nx6`, where `n` represents the number of CTA's in a particular region for which at least one oscillation with period `T0,...,T0/5` was detected (before filtering). Note that `n` is not necessarily equal to the total number of CTA's in that region because it is possible that no oscilattions are detected for some of the CTA's. In the first column, it is saved the name of the CTA. In the `j`-th column (`j=2,..,6`) are saved the estimated periods that correspond to the oscillation with period `T0/(j-1)`. The estimated periods are computed based on the frequencies detected by BIC and SC. |
| **`t_after`**  | Similar to `t_before`, except that it stores the results obtained after filtering. Note that the filter is applied only to the time series data from the CTA's where the fundamental period was detected. |

### run_filt_fullsearch

For the Full-Search method, you should run the function run_filt_fullsearch.m. The input arguments are explained below.

**Input**

| **Variable**   | **Description**                                                                                         |
|----------------|---------------------------------------------------------------------------------------------------------|
| **`%wreg`**    | WHO region for which the data is analyzed                                                               |
| **`%flag_transf`** | Flag for data transform:<br> 0 = No transform applied to the data<br> 1 = Transform from GP Nason, Scientific Reports, 2020<br> 2 = Transform usually applied to daily stock market indices at closing time<br> 3 = First order differences  |
| **`%sp`**      | 'day' = Daily data<br> 'week' = weekly data                                                             |
| **`%ftype`**   | 'stop' = Stop-band filter<br> 'high' = High-pass filter                                                |
| **`%fdomain`** | 'time' = Filtering in time domain<br> 'freq' = Filtering in frequency domain                           |


**Output**

The output generated by this function is very similar to that of run_filt_WHO.m, except that the .mat file contains only the variables seven, t_before, and t_after.

## Data management
- [mean_sd.m](#mean_sd)
- [data_tidy.Rmd](#data_tidy)
- [generate_plot.Rmd](#generate_plot)
- [patition_comparison.Rmd](#partition_index)

### mean_sd
The means and the standard deviations of the estimation errors (before and after filtering) for each WHO region. The input arguments for the function mean_sd are explained below.

**Input**
| **Variable**   | **Description**                                                                                         |
|----------------|---------------------------------------------------------------------------------------------------------|
| **`%wreg`**    | WHO region for which the data is analyzed                                                               |
| **`%flag_transf`** | Flag for data transform:<br> 0 = No transform is applied to the data<br> 1 = Transform from GP Nason, Scientific Reports, 2020<br> 2 = Transform usually applied to daily stock market indices at closing time<br> 3 = First order differences |
| **`%sp`**      | 'day' = Daily data<br> 'week' = weekly data                                                             |
| **`%ftype`**   | 'stop' = Stop-band filter<br> 'high' = High-pass filter                                                |
| **`%fdomain`** | 'time' = Filtering in time domain<br> 'freq' = Filtering in frequency domain                           |

**Output**

**`%results_mean_before`** 

A vector of size `5x1`; the `j`-th entry (`j=1,...,5`) is the mean of the estimation errors for the oscillation with period `T0/j` (before filtering)

**`%results_sd_before`** 

A vector of size `5x1`; the `j`-th entry (`j=1,...,5`) is the standard deviation of the estimation errors for the oscillation with period `T0/j` (before filtering)

**`%results_mean_after`** 

A vector of size `5x1`; the `j`-th entry (`j=1,...,5`) is the mean of the estimation errors for the oscillation with period `T0/j` (after filtering)

**`%results_sd_after`** 

A vector of size `5x1`; the `j`-th entry (`j=1,...,5`) is the standard deviation of the estimation errors for the oscillation with period `T0/j` (after filtering)

### data_tidy

`process_lines_daily(lines)`：For daily data - extracts the results obtained before filtering from the TXT files and convert them into wide format. The significance of the binary value written down for each CTA and each fundamental frequency/harmonic: 1 - detected; 0 - not detected 

`process_lines_weekly(lines)`: Similar to the function `process_lines_daily`, except that it is designed for weekly data. 

`process_lines_daily_after(lines)`: Similar to the function `process_lines_daily`, except that it extracts the results obtained after filtering.

`process_lines_week_after(lines)`: Similar to the function `process_lines_daily`, except that it extracts the results obtained after filtering, for weekly data.

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

It is exactly the same as before_filter_daily.csv, with the only difference being that it focuses on confirming the annual cycle (52-week cycle). In addition, there is a column called `latitude_region`.

**`after_filter_daily.csv`**: 

It is exactly the same as before_filter_daily.csv. There are two kinds of filters: high-pass and stop-band, and also two domains: time and frequency.

**`after_filter_daily.csv`**: 

It is exactly the same as before_filter_week.csv.

### generate_plot
Plot the results of the experiments. Use ggplot to aggregate and summarize the daily data and weekly data by WHO region; summarize the weekly data by Latitude group.

| **Type**   | **Files**                                                                                             |
|------------|-------------------------------------------------------------------------------------------------------|
| **Input**  | `before_filter_daily.csv`<br>`after_filter_daily.csv`<br>`before_filter_week.csv`<br>`after_filter_week.csv` |
| **Output** | `figure_3.pdf`<br>`figure_4.pdf`<br>`figure_5.pdf`<br>`figure_6.pdf`<br>`figure_7.pdf`<br>`figure_8.pdf`<br>`figure_B1.pdf`<br>`figure_B2.pdf`<br>`figure_B3.pdf`<br>`figure_B4.pdf` |


### partition_index
We use the ‘partitionComparison’ package in R environment for computing the following similarity indexes for the two partitions: Rand, Jaccard and Fowlkes-Mallows ). 

**Input**

**`country_label.xlsx`**

**`Latitude_label.xlsx`**

**Output**

Rand, Jaccard and Fowlkes-Mallows Index


The full report can be found at [Link]().
